## Modified calculate_morphological_distances() function from 'Claddis' package
## that silences line [original on GitHub] 416 to avoid triggering
## 'character_dependencies' error (which in the morphological data set case is a
## false error because the dependencies are not always straightforward when
## coding phylum-wide characters.) Code by Graeme Lloyd and taken from
## https://github.com/graemetlloyd/Claddis/blob/master/R/calculate_morphological_distances.R
## , committed 4/1/2021. Thanks, Graeme!

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 8/5/2021, v.
# 0.6.3. Users familiar with earlier versions will need to either download the
# archived version of the package from GitHub or alter the code accordingly to
# use appropriate argument and function names.

calculate_morphological_distances2 <- function(cladistic_matrix, distance_metric = "mord", ged_type = "wills", distance_transformation = "arcsine_sqrt", polymorphism_behaviour = "min_difference", uncertainty_behaviour = "min_difference", inapplicable_behaviour = "missing", character_dependencies = NULL, alpha = 0.5) {
  
  # ADD HOPKINS SUGGESTION (VIA EMAIL) FOR FOURTH GEDTYPE WHERE MEAN DISTANCE FOR CHARACTER REPLACES MISSING VALUES.
  # CHECK POLYMORPHISM UNCERTAINTY IN GENERAL AS NOT CLEAR IT IS DOING WHAT IT SHOULD DO.
  # CHECK TRANSFORM IS APPROPRIATE AND WARN USER IF NOT
  # MAYBE ALLOW MANHATTAN TYPE DISTANCES TOO.
  # ADD LEHMANN REFERENCE!
  # ALLOW MANHATTAN DISTANCES
  # CONSIDER DISTANCES FOR POLYMORPHISMS IN SAME WAY PHYTOOLS DOES WITH POLYMK
  # ADD TRANSORMATION USED TO OUTPUT (AS MAY CHANGE IF OPTIONS COLLIDE)
  # ALLOW WMPD SOMEHOW? MAYBE A SEPARATE FUNCTION WITH GROUPS (WHICH NEEDS TO BE IMPLEMENTED ACROSS THE PACKAGE FOR DISPARITY PLOTS ETC.)
  # RETOOL AROUND STOCHASTIC CHARACTER MAPS IF DOING PHYLOGENY TOO
  # CHECK HSJ MAKES SENSE IN CONTEXT OF COMPARABLE WEIGHTS - MIGHT NEED A THEORETICAL SOLUTION BEFORE AN IMPLEMENTATION ONE.
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Subfunction to find comparable characters for a pairwise taxon comparison:
  find_comparable <- function(taxon_pair, cladistic_matrix) {
    
    # Get intersection of characters that are coded for both taxa in a pair:
    output <- intersect(intersect(which(x = !is.na(cladistic_matrix[taxon_pair[[1]], ])), which(x = cladistic_matrix[taxon_pair[[1]], ] != "")), intersect(which(x = !is.na(cladistic_matrix[taxon_pair[[2]], ])), which(x = cladistic_matrix[taxon_pair[[2]], ] != "")))
    
    # Return output:
    list(output)
  }
  
  # Subfunction to get character strings for each pair of taxa:
  get_pairwise_strings <- function(taxon_pair, cladistic_matrix) {
    
    # Get character states for first taxon in pair:
    row1 <- cladistic_matrix[rownames(x = cladistic_matrix)[taxon_pair[[1]]], ]
    
    # Get character states for second taxon in pair:
    row2 <- cladistic_matrix[rownames(x = cladistic_matrix)[taxon_pair[[2]]], ]
    
    # Return output as a list:
    list(row1, row2)
  }
  
  # Subfunction to subset pairwise comparisons by just comparable characters:
  subset_by_comparable <- function(row_pair, comparable_characters) {
    
    # Collapse first row to just comparable characters:
    row_pair[[1]] <- row_pair[[1]][comparable_characters]
    
    # Collapse second row to just comparable characters:
    row_pair[[2]] <- row_pair[[2]][comparable_characters]
    
    # Output colapsed row pair:
    row_pair
  }
  
  # Subfunction to edit polymorphic characters down to a single value:
  edit_polymorphisms <- function(comparisons, comparable_characters, ordering, polymorphism_behaviour, uncertainty_behaviour) {
    
    # Set first taxon values:
    first_row <- comparisons[[1]]
    
    # Set second taxon values:
    second_row <- comparisons[[2]]
    
    # If there are any inapplicables:
    if (any(c(first_row, second_row) == "")) {
      
      # Find inapplicable positions:
      inapplicable_positions <- sort(x = unique(x = c(which(x = first_row == ""), which(x = second_row == ""))))
      
      # Find polymorphism and uncertainty positions:
      polymorphism_and_uncertainty_positions <- sort(x = unique(x = c(grep("/|&", first_row), grep("/|&", second_row))))
      
      # If there are polymorphisms or uncertianties that match up with inapplicables:
      if (length(x = intersect(inapplicable_positions, polymorphism_and_uncertainty_positions)) > 0) {
        
        # Find positions where collapsing to a single value is required:
        collapse_positions <- intersect(inapplicable_positions, polymorphism_and_uncertainty_positions)
        
        # Collapse any polymorphisms or uncertianties in first row to just first value:
        first_row[collapse_positions] <- unlist(x = lapply(X = strsplit(first_row[collapse_positions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))
        
        # Collapse any polymorphisms or uncertianties in second row to just first value:
        second_row[collapse_positions] <- unlist(x = lapply(X = strsplit(second_row[collapse_positions], split = "/|&"), function(x) ifelse(length(x = x) == 0, "", x[1])))
      }
    }
    
    # Set ordering for comparable characters:
    character_ordering <- ordering[comparable_characters]
    
    # Only if there are polymorphisms or uncertainties:
    if (length(x = c(grep("&", unique(x = c(first_row, second_row))), grep("/", unique(x = c(first_row, second_row))))) > 0) {
      
      # Find ampersands (polymorphisms):
      ampersand_elements <- sort(x = c(grep("&", first_row), grep("&", second_row)))
      
      # Find slashes (uncertianties):
      slash_elements <- sort(x = c(grep("/", first_row), grep("/", second_row)))
      
      # Combine to find all characters to check:
      characters_to_check <- sort(x = unique(x = c(ampersand_elements, slash_elements)))
      
      # Set behaviours as either the shared version or minimum difference if they contradict (may need to modify this later for more complex options):
      behaviour <- unlist(x = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(apply(apply(rbind(first_row[characters_to_check], second_row[characters_to_check]), 2, gsub, pattern = "[:0-9:]", replacement = ""), 2, list), unlist), function(x) x[nchar(x = x) > 0]), function(x) ifelse(nchar(x = x) > 0, strsplit(x, split = "")[[1]][1], x)), function(x) gsub(pattern = "&", replacement = polymorphism_behaviour, x = x)), function(x) gsub(pattern = "/", replacement = uncertainty_behaviour, x = x)), unique), function(x) ifelse(length(x) > 1, "min_difference", x)))
      
      # If behaviour is to find minimum differences:
      if (any(behaviour == "min_difference")) {
        
        # Set up minimum difference characters to check:
        min_characters_to_check <- characters_to_check[behaviour == "min_difference"]
        
        # Find intersecting character states for each character:
        intersection_character <- lapply(X = lapply(X = lapply(X = lapply(X = apply(rbind(first_row[min_characters_to_check], second_row[min_characters_to_check]), 2, strsplit, split = "&|/"), unlist), sort), rle), function(x) x$values[x$lengths > 1][1])
        
        # If at least one intersecting character state was found:
        if (any(!is.na(unlist(x = intersection_character)))) {
          
          # Record rows to update:
          rows_to_update <- which(x = !is.na(unlist(x = intersection_character)))
          
          # Store (first) shared state for both taxa:
          first_row[min_characters_to_check[rows_to_update]] <- second_row[min_characters_to_check[rows_to_update]] <- unlist(x = intersection_character)[rows_to_update]
          
          # Update minimum characters to check:
          min_characters_to_check <- min_characters_to_check[-rows_to_update]
        }
        
        # Only continue if there are still characters that need to be fixed:
        if (length(x = min_characters_to_check) > 0) {
          
          # Build two option matrices for every comparison:
          two_option_matrices <- lapply(X = apply(rbind(first_row[min_characters_to_check], second_row[min_characters_to_check]), 2, strsplit, split = "&|/"), function(x) rbind(c(min(as.numeric(x[[1]])), max(as.numeric(x[[2]]))), c(max(as.numeric(x[[1]])), min(as.numeric(x[[2]])))))
          
          # Pick smallest difference as minimum and maximum states:
          min_max_states <- lapply(X = lapply(X = lapply(X = two_option_matrices, function(x) x[which(x = abs(apply(x, 1, diff)) == min(abs(apply(x, 1, diff)))), ]), sort), as.character)
          
          # Set first row values(s):
          first_row[min_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 1))
          
          # Set second row values(s):
          second_row[min_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 2))
        }
      }
      
      # If any behaviour is to find mean differences:
      if (any(behaviour == "mean_difference")) {
        
        # Set up minimum difference characters to check:
        mean_characters_to_check <- characters_to_check[behaviour == "mean_difference"]
        
        # Build initial state matrices with column and row names as states for first and second rows:
        state_matrices <- lapply(X = lapply(X = apply(rbind(first_row[mean_characters_to_check], second_row[mean_characters_to_check]), 2, list), lapply, strsplit, split = "&|/"), function(x) matrix(nrow = length(x = x[[1]][[1]]), ncol = length(x = x[[1]][[2]]), dimnames = list(x[[1]][[1]], x[[1]][[2]])))
        
        # Fill state matrices with raw differences between each state:
        state_matrices <- lapply(X = state_matrices, function(x) {
          for (i in 1:ncol(x)) for (j in 1:nrow(x)) x[j, i] <- abs(as.numeric(colnames(x = x)[i]) - as.numeric(rownames(x = x)[j]))
          return(x)
        })
        
        # If there are unordered characters present convert maximum distances to one:
        if (any(character_ordering[mean_characters_to_check] == "unordered")) {
          state_matrices[which(x = character_ordering[mean_characters_to_check] == "unordered")] <- lapply(X = state_matrices[which(x = character_ordering[mean_characters_to_check] == "unordered")], function(x) {
            x[x > 1] <- 1
            return(x)
          })
        }
        
        # Extract minimum and maximum states from each matrix with maximum being the mean distance:
        min_max_states <- lapply(X = lapply(X = lapply(X = state_matrices, as.vector), mean), function(x) c(0, x))
        
        # Set first row values(s):
        first_row[mean_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 1))
        
        # Set second row values(s):
        second_row[mean_characters_to_check] <- unlist(x = lapply(X = min_max_states, "[[", 2))
      }
    }
    
    # Return the first and second rows either without polymorphisms or with them removed:
    return(list(first_row, second_row))
  }
  
  # Subfunction to get the absolute difference between the two rows:
  calculate_absolute_difference <- function(column) {
    
    # Isolate first row values:
    first_row <- column[[1]]
    
    # Isolate second row values:
    second_row <- column[[2]]
    
    # Get absolute differences between each pair of characters:
    return(list(abs(as.numeric(first_row) - as.numeric(second_row))))
  }
  
  # Subfunction to correct unordered distances to one:
  fix_unordered <- function(differences, comparable_characters, ordering) {
    
    # If unordered and distance greater than one replace with one:
    if (length(x = which(x = differences > 1)) > 0) differences[which(x = differences > 1)[which(x = ordering[comparable_characters[which(x = differences > 1)]] == "unordered")]] <- 1
    
    # Return corrected unordered distances:
    return(list(differences))
  }
  
  # Subfunction to find incomparable characters:
  find_incomparable <- function(comparable_characters, cladistic_matrix) {
    setdiff(x = 1:ncol(cladistic_matrix), y = comparable_characters)
  }
  
  # Subfunction to get weighted differences:
  weigh_differences <- function(differences, comparable_characters, character_weights) {
    list(as.numeric(character_weights[comparable_characters]) * differences)
  }
  
  # Subfunction to get raw Euclidean distance:
  calculate_red <- function(differences) {
    dist(rbind(differences, rep(0, length(x = differences))), method = "euclidean")
  }
  
  # Subfunction to find maximum possible differences for the comparable characters:
  find_maximum_difference <- function(comparable_characters, maximum_values, minimum_values) {
    as.numeric(maximum_values[comparable_characters]) - as.numeric(minimum_values[comparable_characters])
  }
  
  # Subfunction to transform list of distances into an actual distance matrix:
  convert_list_to_matrix <- function(list, cladistic_matrix, diag = NULL) {
    
    # Set the number of rows:
    k <- nrow(cladistic_matrix)
    
    # Create the empty matrix:
    matrix_output <- matrix(ncol = k, nrow = k)
    
    # Fill up the lower triangle:
    matrix_output[lower.tri(matrix_output)] <- unlist(x = list)
    
    # Make the matrix a distance matrix (both triangles have the same values):
    matrix_output <- as.matrix(as.dist(matrix_output))
    
    # If no diagonal is supplied:
    if (is.null(diag)) {
      
      # Set diagonal as zero:
      diag(x = matrix_output) <- 0
      
      # If a diagonal is supplied:
    } else {
      
      # Add supplied diagonal as diagonal:
      diag(x = matrix_output) <- diag
    }
    
    # Return matrix:
    matrix_output
  }
  
  # Subfunction to get count of complete characters for each taxon (diagonal in comparable characters matrix:
  count_complete <- function(column) {
    length(x = column) - sum(is.na(column))
  }
  
  # Subfunction to calculate the Gower Coefficient:
  calculate_gc <- function(differences, comparable_characters, character_weights) {
    sum(differences) / sum(character_weights[comparable_characters])
  }
  
  # Subfunction to calculate MORD:
  calculate_mord <- function(differences, maximum_differences) {
    sum(differences) / sum(maximum_differences)
  }
  
  # Subfunction for building starting GED data:
  build_ged_data <- function(differences, comparable_characters, cladistic_matrix, character_weights) {
    rbind(c(differences, rep(NA, length(x = find_incomparable(comparable_characters, cladistic_matrix)))), c(character_weights[comparable_characters], character_weights[find_incomparable(comparable_characters, cladistic_matrix)]))
  }
  
  # Subfunction to apply Hopkins and St John (2018) Alpha weighting of inapplicables:
  weigh_inapplicable_alpha <- function(differences, comparable_characters, ordering, character_weights, character_dependencies, characters_by_level, alpha) {
    
    # Set ordering for comparable characters:
    character_ordering <- ordering[comparable_characters]
    
    # Set weights for comparable characters:
    character_weights <- character_weights[comparable_characters]
    
    # Fof each character level (from most to least nested):
    for (i in length(x = characters_by_level):2) {
      
      # Get independent characters for current levels dependent characters:
      independent_characters <- unique(x = unlist(x = lapply(X = as.list(x = characters_by_level[[i]]), function(x) unname(character_dependencies[character_dependencies[, "dependent_character"] == x, "independent_character"]))))
      
      # For each independent character:
      for (j in independent_characters) {
        
        # Find dependent characters:
        dependent_characters <- unname(character_dependencies[character_dependencies[, "independent_character"] == j, "dependent_character"])
        
        # Check characters are present in current distance:
        present_characters <- intersect(comparable_characters, dependent_characters)
        
        # If characters are present:
        if (length(x = present_characters) > 0) {
          
          # Set positions of dependent characters in current differences vector:
          dependent_positions <- match(present_characters, comparable_characters)
          
          # Get position of independent character in current differences vector:
          independent_position <- which(x = comparable_characters == j)
          
          # Stop and warn user if matrix contains an impossible coding (i.e., dependent character coded when independent character is missing):
          # *** THIS IS THE ERROR THAT GETS TRIGGERED: NOW SILENCED ***
          # if (length(x = independent_position) == 0) stop("Found a dependent character coded when character it depends on is missing. Check matrix codings.")
          
          # Overwrite independent position with alpha-weighted value:
          differences[independent_position] <- 1 - (alpha * (1 - (sum(differences[dependent_positions] * character_weights[dependent_positions]) / sum(character_weights[dependent_positions]))) + (1 - alpha))
          
          # Overwrite dependent positions with NAs:
          differences[dependent_positions] <- NA
        }
      }
    }
    
    # Return modified character comparisons:
    differences
  }
  
  # Check for step matrices and stop and warn user if found:
  if (is.list(cladistic_matrix$topper$step_matrices)) stop("Function cannot currently deal with step matrices.")
  
  # Check input of distance_transformation is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_transformation, y = c("arcsine_sqrt", "none", "sqrt"))) > 0) stop("distance_transformation must be one of \"none\", \"sqrt\", or \"arcsine_sqrt\".")
  
  # Check input of distance is valid and stop and warn if not:
  if (length(x = setdiff(x = distance_metric, y = c("red", "ged", "gc", "mord"))) > 0) stop("distance_metric must be one or more of \"red\", \"ged\", \"gc\", or \"mord\".")
  
  # Check input of GED type is valid and stop and warn if not:
  if (length(x = setdiff(x = ged_type, y = c("legacy", "hybrid", "wills"))) > 0) stop("ged_type must be one or more of \"legacy\", \"hybrid\", or \"wills\".")
  
  # Check input for polymorphism_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = polymorphism_behaviour, y = c("mean_difference", "min_difference", "random"))) > 0) stop("polymorphism_behaviour must be one or more of \"mean_difference\", \"min_difference\", or \"random\".")
  
  # Check input for uncertainty_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = uncertainty_behaviour, y = c("mean_difference", "min_difference", "random"))) > 0) stop("uncertainty_behaviour must be one or more of \"mean_difference\", \"min_difference\", or \"random\".")
  
  # Check input for inapplicable_behaviour is valid and stop and warn if not:
  if (length(x = setdiff(x = inapplicable_behaviour, y = c("missing", "hsj"))) > 0) stop("inapplicable_behaviour must be one or more of \"missing\", or \"hsj\".")
  
  # Check that if using HSJ character dependencies have been specified:
  if (inapplicable_behaviour == "hsj" && is.null(character_dependencies)) stop("If using the \"hsj\" inapplicable_behaviour then character_dependencies must be specified.")
  
  # If using HSJ and character_dependencies is set (will check data are formatted correctly):
  if (inapplicable_behaviour == "hsj" && !is.null(character_dependencies)) {
    
    # Check character_dependencies is a matrix and stop and warn user if not:
    if (!is.matrix(character_dependencies)) stop("character_dependencies must be in the form of a two-column matrix.")
    
    # Check character_dependencies has two columns and stop and warn user if not:
    if (ncol(character_dependencies) != 2) stop("character_dependencies must be in the form of a two-column matrix.")
    
    # Check character_dependencies column names are correct and stop and warn user if not:
    if (length(x = setdiff(x = c("dependent_character", "independent_character"), y = colnames(x = character_dependencies))) > 0) stop("character_dependencies column names must be exactly \"dependent_character\" and \"independent_character\".")
    
    # Check character_dependencies are numeric values and stop and warn user if not:
    if (!is.numeric(character_dependencies)) stop("character_dependencies values must be numeric.")
    
    # Check character_dependencies values are within range of matrix dimensions and stop and warn user if not:
    if (length(x = setdiff(x = as.vector(character_dependencies), y = 1:sum(unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) ncol(x$matrix))))))) > 0) stop("character_dependencies can only contain character numbers within the dimensions of the cladistic_matrix specified.")
    
    # Check character_dependencies values do not lead to duplicated parent characters and stop and warn user if not:
    if (any(duplicated(character_dependencies[, "dependent_character"]))) stop("character_dependencies characters can not be dependent on two or more different independent characters.")
    
    # Find any characters that are both dependent and independent (and hence may lead to circularity issues):
    potentially_circular_dependencies <- intersect(character_dependencies[, "dependent_character"], character_dependencies[, "independent_character"])
    
    # If there is the possibility for circularity:
    if (length(x = potentially_circular_dependencies) > 0) {
      
      # For the ith independent character:
      for (i in unique(x = character_dependencies[, "independent_character"])) {
        
        # Set current character as ith character:
        current_character <- i
        
        # Ste starting found character as ith character:
        found_characters <- i
        
        # Keep going until the current character is not an independent character:
        while (sum(unlist(x = lapply(X = as.list(x = current_character), function(x) sum(character_dependencies[, "independent_character"] == x)))) > 0) {
          
          # Find any dependent character(s):
          dependent_character <- unlist(x = lapply(X = as.list(x = current_character), function(x) unname(character_dependencies[character_dependencies[, "independent_character"] == x, "dependent_character"])))
          
          # Check character was not already found (creating a circularity) and stop and wanr user if true:
          if (length(x = intersect(dependent_character, found_characters)) > 0) stop("Circularity found in character_dependencies. Fix and try again.")
          
          # Update found characters:
          found_characters <- c(found_characters, dependent_character)
          
          # Update current character(s):
          current_character <- dependent_character
        }
      }
    }
    
    # Check alpha is a value between zero and one and stop and warn user if not:
    if (alpha > 1 || alpha < 0) stop("alpha must be a value between zero and one")
  }
  
  # Isolate ordering element:
  ordering <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))
  
  # Isolate minimum values:
  minimum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))
  
  # Isolate maximum values:
  maximum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))
  
  # Isolate weights:
  character_weights <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))
  
  # Combine matrix blocks into a single matrix:
  cladistic_matrix <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))
  
  # If polymorphism_behaviour is to randomly sample one state:
  if (polymorphism_behaviour == "random") {
    
    # Find cells with polymorphisms:
    polymorphism_cells <- grep("&", cladistic_matrix)
    
    # If there are polymorphisms randomly sample one value and store:
    if (length(x = polymorphism_cells) > 0) cladistic_matrix[polymorphism_cells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[polymorphism_cells]), function(x) sample(strsplit(x, split = "&")[[1]], size = 1)))
    
    # Reset behaviour as mean difference to allow it to interact correctly with uncertainty_behaviour later:
    polymorphism_behaviour <- "mean_difference"
  }
  
  # If uncertainty_behaviour is to randomly sample one state:
  if (uncertainty_behaviour == "random") {
    
    # Find cells with uncertainties:
    uncertainty_cells <- grep("/", cladistic_matrix)
    
    # If there are uncertainties randomly sample one value and store:
    if (length(x = uncertainty_cells) > 0) cladistic_matrix[uncertainty_cells] <- unlist(x = lapply(X = as.list(x = cladistic_matrix[uncertainty_cells]), function(x) sample(strsplit(x, split = "/")[[1]], size = 1)))
    
    # Reset behaviour as mean difference to allow it to interact correctly with polymorphism_behaviour later:
    uncertainty_behaviour <- "mean_difference"
  }
  
  # If there are inapplicables and using the missing option then convert these to NAs:
  if (any(sort(x = cladistic_matrix == "")) && inapplicable_behaviour == "missing") cladistic_matrix[cladistic_matrix == ""] <- NA
  
  # Find all possible (symmetric) pairwise comparisons for the N taxa in the matrix (excluding self-comparisons):
  comparisons <- combn(1:nrow(cladistic_matrix), 2)
  
  # Find all comparable characters for each pair of taxa:
  comparable_character_list <- unlist(x = apply(comparisons, 2, find_comparable, cladistic_matrix), recursive = FALSE)
  
  # Get character states for each pairwise comparison:
  rows_pairs <- apply(comparisons, 2, get_pairwise_strings, cladistic_matrix)
  
  # Subset each pairwise comparison by just the comparable characters:
  matrix_of_comparable_characters <- mapply(subset_by_comparable, rows_pairs, comparable_character_list)
  
  # Deal with any polymorphisms found and collapse appropriately:
  matrix_of_comparable_characters <- mapply(edit_polymorphisms, unlist(x = apply(matrix_of_comparable_characters, 2, list), recursive = FALSE), comparable_character_list, MoreArgs = list(ordering, polymorphism_behaviour, uncertainty_behaviour))
  
  # Get the absolute differences between each comparable character for each pairwise comparison:
  absolute_differences <- unlist(x = apply(matrix_of_comparable_characters, 2, calculate_absolute_difference), recursive = FALSE)
  
  # Correct distances for unordered characters where distance is greater than one:
  absolute_differences <- mapply(fix_unordered, absolute_differences, comparable_character_list, MoreArgs = list(ordering))
  
  # If applying the Hopkins and St John alpha approach:
  if (inapplicable_behaviour == "hsj") {
    
    # Set primary-level characters in a list (where secondary etc. level characters will be added in turn):
    characters_by_level <- list(unname(setdiff(x = unique(x = character_dependencies[, "independent_character"]), y = unique(x = character_dependencies[, "dependent_character"]))))
    
    # Set starting more nested characters:
    higher_level_characters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = characters_by_level))
    
    # Whilst there are still more nested levels of characters:
    while (length(x = higher_level_characters) > 0) {
      
      # Add next level characters to characters by level list at next level:
      characters_by_level[[(length(x = characters_by_level) + 1)]] <- unname(character_dependencies[unlist(x = lapply(X = as.list(x = characters_by_level[[length(x = characters_by_level)]]), function(x) which(x = character_dependencies[, "independent_character"] == x))), "dependent_character"])
      
      # Set new higher level characters:
      higher_level_characters <- setdiff(x = unique(x = c(character_dependencies)), y = unlist(x = characters_by_level))
    }
    
    # Update differences with HSJ alpha weights:
    absolute_differences <- mapply(weigh_inapplicable_alpha, absolute_differences, comparable_character_list, MoreArgs = list(ordering, character_weights, character_dependencies, characters_by_level, alpha))
    
    # Reweight dependent characters zero:
    character_weights[unlist(x = characters_by_level[2:length(x = characters_by_level)])] <- 0
    
    # Update comparable characters by pruning out NAs:
    comparable_character_list <- mapply(function(x, y) y[!is.na(x)], x = absolute_differences, y = comparable_character_list, SIMPLIFY = FALSE)
    
    # Update differences by pruning out NAs:
    absolute_differences <- lapply(X = absolute_differences, function(x) x[!is.na(x)])
  }
  
  # Weight differences:
  absolute_differences <- mapply(weigh_differences, absolute_differences, comparable_character_list, MoreArgs = list(character_weights))
  
  # Get raw Euclidean distance (if using it):
  if (distance_metric == "red") raw_distances <- lapply(X = absolute_differences, calculate_red)
  
  # Only calculate the max differences for "ged" or "mord" matrices:
  if (distance_metric == "ged" || distance_metric == "mord") {
    
    # Find maximum possible differences for the comparable characters:
    maximum_possible_differences <- lapply(X = comparable_character_list, find_maximum_difference, maximum_values, minimum_values)
    
    # Correct maximum differences for unordered characters:
    maximum_possible_differences <- mapply(weigh_differences, mapply(fix_unordered, maximum_possible_differences, comparable_character_list, MoreArgs = list(ordering)), comparable_character_list, MoreArgs = list(character_weights))
  }
  
  # If calculating Raw Euclidean Distances build the distance matrix:
  if (distance_metric == "red") distance_matrix <- convert_list_to_matrix(raw_distances, cladistic_matrix)
  
  # If calculating the Gower Coefficient build the distance matrix:
  if (distance_metric == "gc") distance_matrix <- convert_list_to_matrix(as.list(x = mapply(calculate_gc, absolute_differences, comparable_character_list, MoreArgs = list(character_weights))), cladistic_matrix)
  
  # If calculating the MORD build the distance matrix:
  if (distance_metric == "mord") distance_matrix <- convert_list_to_matrix(mapply(calculate_mord, absolute_differences, maximum_possible_differences), cladistic_matrix)
  
  # If calculating the GED:
  if (distance_metric == "ged") {
    
    # Build starting GED data:
    ged_data <- mapply(build_ged_data, absolute_differences, comparable_character_list, MoreArgs = list(cladistic_matrix, character_weights), SIMPLIFY = FALSE)
    
    # Transpose matrices:
    ged_data <- lapply(X = ged_data, t)
    
    # Now build into matrix of pairwise comparisons (odds to be compared with adjacent evens):
    ged_data <- matrix(data = (unlist(x = ged_data)), ncol = ncol(cladistic_matrix), byrow = TRUE)
    
    # Calculate single weighted mean univariate distance for calculating GED Legacy or Hybrid (after equation 2 in Wills 2001):
    if (ged_type != "wills") nonwills_s_ijk_bar <- rep(sum(unlist(x = absolute_differences)) / sum(unlist(x = maximum_possible_differences)), length.out = length(x = absolute_differences))
    
    # Calculate individual pairwise weighted mean univariate distance for calculating GED Hybrid or Wills (after equation 2 in Wills 2001):
    if (ged_type != "legacy") {
      
      # Generate individual mean pairwise distance for each comparison:
      nonlegacy_s_ijk_bar <- unlist(x = lapply(X = absolute_differences, sum)) / unlist(x = lapply(X = maximum_possible_differences, sum))
      
      # Find NaNs (divide by zero errors for when there are no characters in common in a pairwise comparison):
      not_a_number <- which(x = is.nan(nonlegacy_s_ijk_bar))
      
      # If using Wills version replace NaNs with NA:
      if (ged_type == "wills" && length(x = not_a_number) > 0) nonlegacy_s_ijk_bar[not_a_number] <- NA
      
      # If using Hybrid replace NaNs with single global mean distance value:
      if (ged_type == "hybrid" && length(x = not_a_number) > 0) nonlegacy_s_ijk_bar[not_a_number] <- nonwills_s_ijk_bar[not_a_number]
      
      # Set modified non-legacy s_ijk_bar as main s_ijk_bar:
      s_ijk_bar <- nonlegacy_s_ijk_bar
    }
    
    # If using Legacy set nonwills_s_ijk_bar as main s_ijk_bar:
    if (ged_type == "legacy") s_ijk_bar <- nonwills_s_ijk_bar
    
    # For each set of differences:
    for (i in seq(from = 1, to = nrow(ged_data) - 1, length.out = length(x = absolute_differences))) {
      
      # Find missing distances (if any):
      missing_distances <- which(x = is.na(ged_data[i, ]))
      
      # Replace missing distances with s_ijk_bar (i.e., results of equation 2 in Wills 2001 into equation 1 of Wills 2001):
      if (length(x = missing_distances) > 0) ged_data[i, missing_distances] <- s_ijk_bar[ceiling(i / 2)]
    }
    
    # Isolate the distances:
    s_ijk <- ged_data[which(x = (1:nrow(ged_data) %% 2) == 1), ]
    
    # Isolate the weights:
    w_ijk <- ged_data[which(x = (1:nrow(ged_data) %% 2) == 0), ]
    
    # Calculate the GED (equation 1 of Wills 2001) for each pairwise comparison (ij):
    ged_ij <- sqrt(apply(w_ijk * (s_ijk^2), 1, sum))
    
    # Create GED distance matrix:
    distance_matrix <- convert_list_to_matrix(as.list(x = ged_ij), cladistic_matrix)
  }
  
  # Build comparable characters matrix:
  comparable_character_matrix <- convert_list_to_matrix(lapply(X = comparable_character_list, length), cladistic_matrix, diag = apply(cladistic_matrix, 1, count_complete))
  
  # Build comparable weights matrix:
  comparable_weights_matrix <- convert_list_to_matrix(lapply(X = comparable_character_list, FUN = function(x) sum(x = character_weights[x])), cladistic_matrix, diag = unname(obj = apply(X = cladistic_matrix, MARGIN = 1, FUN = function(x) sum(x = character_weights[which(x = !is.na(x = x))]))))
  
  # Add row and column names (taxa) to distance matrices:
  rownames(x = distance_matrix) <- colnames(x = distance_matrix) <- rownames(x = comparable_character_matrix) <- colnames(x = comparable_character_matrix) <- rownames(x = comparable_weights_matrix) <- colnames(x = comparable_weights_matrix) <- rownames(x = cladistic_matrix)
  
  # If there are any NaNs replace with NAs:
  if (any(is.nan(distance_matrix))) distance_matrix[is.nan(distance_matrix)] <- NA
  
  # If using a proportional distance:
  if (distance_metric == "mord" || distance_metric == "gc") {
    
    # If transforming distance matrix by taking the square root - take the square root:
    if (distance_transformation == "sqrt") distance_matrix <- sqrt(distance_matrix)
    
    # If transforming distance matrix by taking the arcsine square root:
    if (distance_transformation == "arcsine_sqrt") {
      
      # Check for squared distances greater than 1:
      if (any(sort(x = sqrt(distance_matrix)) > 1)) {
        
        # Warn user that distances were rescaled:
        print("Squared distances found of greater than 1 so matrix was rescaled prior to taking arcsine.")
        
        # Take the arcsine square root of the rescaled distance matrix:
        distance_matrix <- asin(sqrt(distance_matrix) / max(sort(x = sqrt(distance_matrix))))
        
        # If squared distances are less than or equal to one:
      } else {
        
        # Take the arcsine square root directly:
        distance_matrix <- asin(sqrt(distance_matrix))
      }
    }
  }
  
  # Return compiled output:
  list(distance_metric = distance_metric, distance_matrix = distance_matrix, comparable_character_matrix = comparable_character_matrix, comparable_weights_matrix = comparable_weights_matrix)
}