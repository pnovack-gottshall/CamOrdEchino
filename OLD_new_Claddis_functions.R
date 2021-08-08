## New version of DiscreteCharacterRate() function from 'Claddis' package,
## Modified to accept pre-inferred ancestral states output from
## Claddis::AncStateEstMatrix(). Code by Graeme Lloyd and taken from
## https://github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R.
## Downloaded 8/20/2020. Thanks, Graeme!

## We did not update the previously run code (based on older function and arg
## names from Claddis v. 0.4.1), so the code below also copies all functions
## called within test_rates below, as well as related functions drawn from new
## Claddis package used in manuscript code.

## MODIFIED test_rates() #######################################################
test_rates2 <-
  function(time_tree,
           cladistic_matrix,
           time_bins,
           branch_partitions = NULL,
           character_partitions = NULL,
           clade_partitions = NULL,
           time_partitions = NULL,
           change_times = "random",
           test_type = "aic",
           alpha = 0.01,
           multiple_comparison_correction = "benjaminihochberg",
           polymorphism_state = "missing",
           uncertainty_state = "missing",
           inapplicable_state = "missing",
           time_binning_approach = "lloyd",
           all_weights_integers = FALSE,
           estimate_all_nodes = FALSE,
           estimate_tip_values = FALSE,
           inapplicables_as_missing = FALSE,
           polymorphism_behaviour = "equalp",
           uncertainty_behaviour = "equalp",
           threshold = 0.01,
           all_missing_allowed = FALSE) {
    
  # Check for step matrices and stop and warn user if found:
  if (is.list(cladistic_matrix$topper$step_matrices)) stop("Function cannot currently deal with step matrices.")
  
  # Check tree has branch lengths:
  if (is.null(time_tree$edge.length)) stop("time_tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")
  
  # Check tree has root age:
  if (is.null(time_tree$root.time)) stop("time_tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")
  
  # Check change_times is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = change_times, y = c("midpoint", "spaced", "random"))) > 0) stop("change_times must be one of \"midpoint\", \"spaced\", or \"random\".")
  
  # Check multiple_comparison_correction is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = multiple_comparison_correction, y = c("benjaminihochberg", "bonferroni"))) > 0) stop("multiple_comparison_correction must be one of \"benjaminihochberg\" or \"bonferroni\".")
  
  # Check polymorphism_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = polymorphism_state, y = c("missing", "random"))) > 0) stop("polymorphism_state must be one of \"missing\" or \"random\".")
  
  # Check uncertainty_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = uncertainty_state, y = c("missing", "random"))) > 0) stop("uncertainty_state must be one of \"missing\" or \"random\".")
  
  # Check inapplicable_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = inapplicable_state, y = c("missing"))) > 0) stop("inapplicable_state must be \"missing\".")
  
  # Check time_binning_approach is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = time_binning_approach, y = c("close", "lloyd"))) > 0) stop("time_binning_approach must be one of \"close\" or \"lloyd\".")
  
  # Check partitions are not all NULL values:
  if (is.null(branch_partitions) && is.null(character_partitions) && is.null(clade_partitions) && is.null(time_partitions)) stop("No partitions are requested. Set at least one of branch_partitions, character_partitions, clade_partitions, or time_partitions to a list of appropriate values. Type \"?test_rates\" for help.")
  
  # Get tip number:
  n_tips <- ape::Ntip(phy = time_tree)
  
  # Get node number:
  n_nodes <- ape::Nnode(phy = time_tree)
  
  # Get internal node numbers:
  node_numbers <- 1:n_nodes + n_tips
  
  # Get edge numbers:
  edge_numbers <- 1:nrow(time_tree$edge)
  
  # Get character numbers:
  character_numbers <- 1:sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))
  
  # Ensure time bins are in correct order:
  time_bins <- sort(x = unique(x = time_bins), decreasing = TRUE)
  
  # Find the Time bin midpoints:
  time_bin_midpoints <- (time_bins[2:length(x = time_bins)] + time_bins[1:(length(x = time_bins) - 1)]) / 2
  
  # Get the numbers for each time bins:
  time_bin_numbers <- 1:length(x = time_bin_midpoints)
  
  # Subfunction to ensure partitions are formatted correctly:
  format_partition <- function(partitions_to_test, valid_values, partition_name) {
    
    # Check partitions are in the form of a list of lists:
    if (!all(c(all(unlist(x = lapply(X = partitions_to_test, is.list))), is.list(partitions_to_test)))) stop(paste(partition_name, " must be in the form of a list of lists.", sep = ""))
    
    # Get a vector of any absent valid values:
    absent_values <- setdiff(x = unique(x = unlist(x = partitions_to_test)), y = valid_values)
    
    # Check valid values have been used and if not stop and warn user:
    if (length(x = absent_values) > 0) stop(paste(partition_name, "Partitions to test must be defined using the valid range of values (", paste(range(valid_values), collapse = " to "), ") only.", sep = ""))
    
    # Check partitions never overlap and stop and warn user if they do:
    overlap_check <- lapply(X = partitions_to_test, function(x) if (any(duplicated(sort(x = unlist(x = x))))) stop(paste("Each partition of ", partition_name, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))
    
    # Subfunction to add the missing partition (if exists):
    add_missing_partitions <- function(x, valid_values) {
      
      # Define any missing values:
      missing_values <- setdiff(x = valid_values, y = unlist(x = x))
      
      # If there are missing values add them to list at end:
      if (length(x = missing_values) > 0) x[[(length(x = x) + 1)]] <- missing_values
      
      # Return x:
      x
    }
    
    # Add in missing partitions (if any):
    partitions_to_test <- lapply(X = partitions_to_test, add_missing_partitions, valid_values = valid_values)
    
    # Check partitions are all at least two in size or else no comparison can be made:
    if (any(unlist(x = lapply(X = partitions_to_test, length)) == 1) && test_type == "lrt") stop("Partitions must divide the available data into at least two parts if performing likelihood ratio tests.")
    
    # Return formatted partitions to test:
    return(partitions_to_test)
  }
  
  # Subfunction to pack partitions to short format for output:
  pack_partitions <- function(formatted_partitions) {
    unlist(x = lapply(X = formatted_partitions, function(x) {
      paste(unlist(x = lapply(X = x, function(y) {
        
        # First make sure y is sorted:
        y <- sort(x = y)
        
        # Covnvert y to a list (splitting if gaps greater than 1 are found):
        y <- unname(split(y, cumsum(c(TRUE, diff(y) > 1))))
        
        # Collapse gaps of one with hyphens:
        paste0(unlist(x = lapply(X = y, function(z) {
          res <- as.character(z)
          if (length(x = z) > 1) {
            r <- rle(c(1, pmin(diff(z), 2)))
            res <- paste0(z[c(1, cumsum(r$lengths))], c("-", " ")[r$values], collapse = "")
            res <- substr(res, 1, nchar(x = res) - 1)
          }
          res
        })), collapse = " ")
      })), collapse = " | ")
    }))
  }
  
  # If performing branch partition test(s) check and reformat branch partitions:
  if (!is.null(branch_partitions)) branch_partitions <- format_partition(partitions_to_test = branch_partitions, valid_values = edge_numbers, partition_name = "branch_partitions")
  
  # If performing character partition test(s) check and reformat character partitions:
  if (!is.null(character_partitions)) character_partitions <- format_partition(partitions_to_test = character_partitions, valid_values = character_numbers, partition_name = "character_partitions")
  
  # If performing clade partition test(s)
  if (!is.null(clade_partitions)) {
    
    # Convert clade partitions to edge partitions:
    clade_partitions <- lapply(X = clade_partitions, lapply, find_descendant_edges, tree = time_tree)
    
    # Check and reformat clade partitions:
    clade_partitions <- format_partition(partitions_to_test = clade_partitions, valid_values = edge_numbers, partition_name = "clade_partitions")
  }
  
  # If performing time bin partition test(s) check and reformat time bin partitions:
  if (!is.null(time_partitions)) time_partitions <- format_partition(partitions_to_test = time_partitions, valid_values = time_bin_numbers, partition_name = "time_partitions")
  
  # Check test_type is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = test_type, y = c("aic", "lrt"))) > 0) stop("test_type must be one of \"aic\" or \"lrt\".")
  
  # Subfunction to calculate maximum likelihood p value:
  get_likelihood_p <- function(mean_rate, sampled_rates, sampled_changes, sampled_completeness, sampled_time) {
    
    # Get log numerator:
    log_numerator <- sum(log(dpois(round(sampled_changes), mean_rate * sampled_completeness * sampled_time)))
    
    # Get log denominator:
    log_denominator <- sum(log(dpois(round(sampled_changes), sampled_rates * sampled_completeness * sampled_time)))
    
    # Get test statistic:
    test_statistic <- -2 * (log_numerator - log_denominator)
    
    # Calculate position of test statistic in chi-square distribution to get probability and return:
    pchisq(test_statistic, length(x = sampled_rates) - 1, lower.tail = FALSE)
  }
  
  # Subfunction to calculate AIC:
  calculate_aic <- function(sampled_rates, sampled_changes, sampled_completeness, sampled_time) {
    
    # Get log maximum likelihood estimate:
    log_mle <- sum(log(dpois(x = round(sampled_changes), lambda = sampled_rates * sampled_completeness * sampled_time)))
    
    # Calculate and return AIC:
    (2 * length(x = sampled_rates)) - (2 * log_mle)
  }
  
  # Subfunction to calculate AIC from partition (with columns labelled partition, rate, completeness, duration):
  calculate_partition_aic <- function(partition, aicc = FALSE) {
    
    # Get log maximum likelihood estimate:
    log_mle <- sum(log(dpois(round(partition[, "changes"]), partition[, "rate"] * partition[, "completeness"] * partition[, "duration"])))
    
    # Get k (number of parameters) term:
    k <- max(partition[, "partition"])
    
    # Calculate AIC:
    aic <- (2 * k) - (2 * log_mle)
    
    # If AICc is desired then calculate this and overwrite AIC with it:
    if (aicc) aic <- aic + (((2 * (k^2)) + (2 * k)) / (nrow(partition) - k - 1))
    
    # Return AIC:
    aic
  }
  
  # Get ages for each (tip and internal) node:
  node_dates <- date_nodes(time_tree = time_tree)
  
  # Get branch ages (from and to):
  branch_ages <- unname(cbind(node_dates[as.character(time_tree$edge[, 1])], node_dates[as.character(time_tree$edge[, 2])]))
  
  # Build edge list from node numbers (from-to) for each branch:
  edge_list <- lapply(X = apply(time_tree$edge, 1, list), function(x) {
    names(x) <- "node_number_from_to"
    x
  })
  
  # Add node ages to edge list:
  for (i in 1:length(x = edge_list)) edge_list[[i]]$node_age_from_to <- branch_ages[i, ]
  
  # Add node ages (from-to) to each edge in list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$branch_duration <- x$node_age_from_to[1] - x$node_age_from_to[2]
    x
  })
  
  # Get vector of branch types:
  branch_types <- gsub(pattern = "0", replacement = "internal", x = gsub(pattern = "1", replacement = "terminal", x = as.numeric(time_tree$edge[, 2] <= n_tips)))
  
  # Add branch type to edge list:
  for (i in 1:length(x = edge_list)) edge_list[[i]]$branch_type <- branch_types[i]
  
  # Find descendant edges for each internal node:
  find_descendant_edges_for_each_internal_node <- lapply(X = as.list(x = node_numbers), find_descendant_edges, tree = time_tree)
  
  # Get ancestral character states:
  # *** MODIFIED TO ACCEPT PREVIOUSLY ESTIMATED ANCESTRAL MATRIX
  ancestral_states <- cladistic_matrix
  
  # Build single matrix of all states in tip label then node number order:
  all_states <- do.call(what = cbind, args = lapply(X = lapply(X = ancestral_states[2:length(x = ancestral_states)], "[[", "matrix"), function(x) x[c(time_tree$tip.label, 1:n_nodes + n_tips), , drop = FALSE]))
  
  # Make vector of ordering of characters:
  ordering <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))
  
  # Make vector of weights of characters:
  character_weights <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))
  
  # Make vector of minimum values:
  minimum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))
  
  # Make vector of maximum values:
  maximum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))
  
  # Find positions in matrix with polymorphisms:
  polymorphism_positions <- grep("&", all_states)
  
  # Find positions in matrix with uncertainties:
  uncertainty_positions <- grep("/", all_states)
  
  # Find positions in matrix with inapplicables:
  inapplicable_positions <- which(x = all_states == "")
  
  # If polymorphisms were found:
  if (length(x = polymorphism_positions) > 0) {
    
    # If replacing polymorphsims with missing do so:
    if (polymorphism_state == "missing") all_states[polymorphism_positions] <- NA
    
    # If replacing polymorphisms with random values draw and replace:
    if (polymorphism_state == "random") all_states[polymorphism_positions] <- unlist(x = lapply(X = strsplit(all_states[polymorphism_positions], "&"), sample, size = 1))
  }
  
  # If uncertainties were found:
  if (length(x = uncertainty_positions) > 0) {
    
    # If replacing uncertainties with missing do so:
    if (uncertainty_state == "missing") all_states[uncertainty_positions] <- NA
    
    # If replacing uncertainties with random values draw and replace:
    if (uncertainty_state == "random") all_states[uncertainty_positions] <- unlist(x = lapply(X = strsplit(all_states[uncertainty_positions], "/"), sample, size = 1))
  }
  
  # If inapplicable states were found:
  if (length(x = inapplicable_positions) > 0) {
    
    # If replacing inapplicables with missing do so:
    if (inapplicable_state == "missing") all_states[inapplicable_positions] <- NA
  }
  
  # Set default converted continuous characters to FALSE:
  continuous_characters_discretized <- FALSE
  
  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if (any(ordering == "cont")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Set default converted continuous characters to TRUE:
    continuous_characters_discretized <- TRUE
    
    # Find out which characters are continuous:
    continuous_characters_found <- which(x = ordering == "cont")
    
    # Rescale continous characters as zero to one values:
    list_of_continuous_values_rescaled_zero_to_one <- lapply(X = lapply(X = lapply(X = apply(all_states[, continuous_characters_found, drop = FALSE], 2, list), unlist), as.numeric), function(x) {
      x <- x - min(sort(x = x))
      x <- x / max(sort(x = x))
      return(x)
    })
    
    # Now discretize and store these characters (0 to 31 scale):
    all_states[, continuous_characters_found] <- do.call(what = cbind, args = lapply(X = lapply(X = lapply(X = list_of_continuous_values_rescaled_zero_to_one, function(x) as.list(x = x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x = x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    ordering[continuous_characters_found] <- "ord"
    
    # Convert weights to 1/31:
    character_weights[continuous_characters_found] <- 1 / 31
    
    # Set minimum value to zero:
    minimum_values[continuous_characters_found] <- 0
    
    # Set maximum value to 31:
    maximum_values[continuous_characters_found] <- 31
  }
  
  # If all_weights_integers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if (all_weights_integers) while (is.character(all.equal(sum(character_weights %% 1), 0))) character_weights <- (1 / (character_weights %% 1)[(character_weights %% 1) > 0])[1] * character_weights
  
  # Add from-to node states for each character to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$character_states_from_to <- matrix(all_states[x$node_number_from_to, , drop = FALSE], nrow = 2, dimnames = list(c("from", "to")))
    x
  })
  
  # Subfunction to define character changes:
  build_changes_matrix <- function(x) {
    
    # Find only comparable characters (those scored for both from and to states):
    comparable_characters <- which(x = apply(!apply(x$character_states_from_to, 2, is.na), 2, all))
    
    # Isolate comparable ordering:
    comparable_ordering <- ordering[comparable_characters]
    
    # Isolate comparable weights:
    comparable_weights <- character_weights[comparable_characters]
    
    # Isolate only characters that actually differ (change):
    character_differences <- which(x = x$character_states_from_to["from", comparable_characters] != x$character_states_from_to["to", comparable_characters])
    
    # Build character change matrix:
    character_changes <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("character", "from", "to", "steps", "weight")))
    
    # If characters change then make a matrix from them:
    if (length(x = character_differences) > 0) character_changes <- rbind(character_changes, cbind(as.numeric(comparable_characters[character_differences]), as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]), as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]), ifelse(comparable_ordering[character_differences] == "unord", 1, abs(as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]) - as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]))), comparable_weights[character_differences]))
    
    # Store character changes as new sublist for x:
    x$character_changes <- character_changes
    
    # Store comparable characters as new sublist of x:
    x$comparable_characters <- comparable_characters
    
    # Return x:
    x
  }
  
  # Get character changes and comparable characters and add to edge list:
  edge_list <- lapply(X = edge_list, build_changes_matrix)
  
  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if (!is.null(time_partitions)) {
    
    # Subfunction to add change times to character changes:
    add_change_times <- function(x, change_times) {
      
      # Isolate character changes:
      character_changes <- x$character_changes
      
      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if (any(character_changes[, "steps"] > 1)) {
        
        # Get multistep character changes:
        multistep_characters <- which(x = character_changes[, "steps"] > 1)
        
        # For each multistep character change:
        for (i in rev(multistep_characters)) {
          
          # Isolate other rows:
          other_row_numbers <- setdiff(x = 1:nrow(character_changes), y = i)
          
          # Get unpacked changes (X:Y, e.g., 0:2 would become 0 1 2):
          unpacked_changes <- character_changes[i, "from"]:character_changes[i, "to"]
          
          # Update character changes with multistep changes unpacked:
          character_changes <- rbind(character_changes[other_row_numbers, ], unname(cbind(rep(character_changes[i, "character"], length.out = length(x = unpacked_changes) - 1), unpacked_changes[1:(length(x = unpacked_changes) - 1)], unpacked_changes[2:length(x = unpacked_changes)], rep(1, length.out = length(x = unpacked_changes) - 1), rep(character_changes[i, "weight"], length.out = length(x = unpacked_changes) - 1))))
        }
        
        # Resort character changes by character number:
        character_changes <- character_changes[order(character_changes[, "character"]), ]
      }
      
      # If using midpoint option set character change times as midpoint of branch:
      if (change_times == "midpoint") character_changes <- cbind(character_changes, rep(x$node_age_from_to[1] - (x$branch_duration / 2), length.out = nrow(character_changes)))
      
      # If using spaced then set character change times as equally spaced along branch:
      if (change_times == "spaced") character_changes <- cbind(character_changes, x$node_age_from_to[1] - (seq(from = 0, to = x$branch_duration, length.out = nrow(character_changes) + 1)[1:nrow(character_changes)] + (diff(seq(from = 0, to = x$branch_duration, length.out = nrow(character_changes) + 1))[1] / 2)))
      
      # If using random then set character change times as random draws from a uniform distribution:
      if (change_times == "random") character_changes <- cbind(character_changes, x$node_age_from_to[1] - stats::runif(n = nrow(character_changes), min = 0, max = x$branch_duration))
      
      # Add column name to change time column:
      colnames(x = character_changes)[ncol(character_changes)] <- "time"
      
      # Subfunction to re-sort character change times so they occur in correct order:
      sort_change_times <- function(character_changes) {
        
        # Sort change time for each character from oldest (first) to youngest (last) and store it:
        character_changes[, "time"] <- unname(unlist(x = lapply(X = as.list(x = unique(x = character_changes[, "character"])), function(x) sort(x = character_changes[which(x = character_changes[, "character"] == x), "time"], decreasing = TRUE))))
        
        # Return sorted character changes:
        character_changes
      }
      
      # Re-sort character change times so they occur in correct order:
      character_changes <- sort_change_times(character_changes)
      
      # Add bin for character change as last column:
      character_changes <- cbind(character_changes, unlist(x = lapply(X = as.list(x = character_changes[, "time"]), function(x) max(which(x = x <= time_bins)))))
      
      # Add column name to change time column:
      colnames(x = character_changes)[ncol(character_changes)] <- "bin"
      
      # Overwrite character changes with new version with changes added:
      x$character_changes <- character_changes
      
      # Return x:
      x
    }
    
    # Add character change times to edge list:
    edge_list <- lapply(X = edge_list, add_change_times, change_times = change_times)
  }
  
  # Subfunction to get edge sections in time bins:
  get_edge_sections_in_bins <- function(x, time_bins = time_bins) {
    
    # Set first appearance datum of edge:
    fad <- x$node_age_from_to[1]
    
    # Set last appearance datum of edge:
    lad <- x$node_age_from_to[2]
    
    # Get any time bin boundaries crossed (can be empty if none are):
    boundaries_crossed <- time_bins[2:(length(x = time_bins) - 1)][intersect(which(x = time_bins[2:(length(x = time_bins) - 1)] > lad), which(x = time_bins[2:(length(x = time_bins) - 1)] < fad))]
    
    # If boundaries are crossed:
    if (length(x = boundaries_crossed) > 0) {
      
      # Break up branch into binned sections as vector of FADs:
      fad <- c(fad, boundaries_crossed)
      
      # Break up branch into binned sections as vector of LADs:
      lad <- c(boundaries_crossed, lad)
    }
    
    # Build matrix of branch sections with FADs and LADs:
    branch_sections <- rbind(fad, lad)
    
    # Add bin number present in to column names:
    colnames(x = branch_sections) <- unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = branch_sections["fad", ]), "<=", time_bins), which), max))
    
    # Add new list section for branch (edge) sections binned by time:
    x$binned_edge_sections <- branch_sections
    
    # Return output:
    x
  }
  
  # Get edge sections in time bins:
  edge_list <- lapply(X = edge_list, get_edge_sections_in_bins, time_bins = time_bins)
  
  # Add binned branch durations to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    branch_durations <- rep(0, length(x = time_bins) - 1)
    branch_durations[as.numeric(colnames(x = x$binned_edge_sections))] <- abs(apply(x$binned_edge_sections, 2, diff))
    x$binned_branch_durations <- branch_durations
    x
  })
  
  # Add proportional binned branch lengths to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$proportional_binned_edge_durations <- x$binned_branch_durations / sum(x$binned_branch_durations)
    x
  })
  
  # Start to build matrix of all changes with list of character changes:
  inferred_character_changes <- lapply(X = edge_list, function(x) x$character_changes)
  
  # Add edge number to each matrix of character changes:
  for (i in 1:length(x = inferred_character_changes)) inferred_character_changes[[i]] <- cbind(rep(i, times = nrow(inferred_character_changes[[i]])), inferred_character_changes[[i]])
  
  # Combine all changes into a single matrix:
  inferred_character_changes <- do.call(what = rbind, args = lapply(X = inferred_character_changes, function(x) {
    colnames(x = x)[1] <- "edge"
    x
  }))
  
  # Remove silly rownames from all changes:
  rownames(x = inferred_character_changes) <- NULL
  
  # Create NULL output variables (to be overwritten if called):
  branch_rates <- clade_rates <- time_rates <- character_rates <- NULL
  
  # If doing some kind of edge test (branch or clade):
  if (!is.null(branch_partitions) || !is.null(clade_partitions)) {
    
    # Get (weighted) number of changes on each edge:
    edge_changes <- unlist(x = lapply(X = edge_list, function(x) sum(x$character_changes[, "steps"] * x$character_changes[, "weight"])))
    
    # Get completeness for each edge:
    edge_completeness <- unlist(x = lapply(X = edge_list, function(x) sum(character_weights[x$comparable_characters]) / sum(character_weights)))
    
    # Get duration of each edge:
    edge_durations <- unlist(x = lapply(X = edge_list, function(x) x$branch_duration))
    
    # Set global rate:
    global_rate <- sum(edge_changes) / sum(edge_completeness * edge_durations)
    
    # If performing branch partition tests:
    if (!is.null(branch_partitions)) {
      
      # Create branch rates for output:
      branch_rates <- lapply(X = list(as.list(x = 1:nrow(time_tree$edge))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]
      
      # Add edge numbers and rates:
      branch_rates <- cbind(edge = 1:nrow(time_tree$edge), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = branch_rates[, "changes"] / branch_rates[, "completeness"])), branch_rates)
      
      # If using Likelihood Ratio Test:
      if (test_type == "lrt") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        partitioned_data <- lapply(X = branch_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))
        
        # Add sampled rate to paritioned data matrices:
        partitioned_data <- lapply(X = partitioned_data, function(x) {
          cbind(rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        })
        
        # Get LRT p-values and combine output as edge test results:
        branch_test_results <- lapply(X = partitioned_data, function(x) {
          list(rates = x[, "rate"], p_value = get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        })
      }
      
      # If using AIC:
      if (test_type == "aic") {
        
        # Build partitioned data for AIC:
        partitioned_data <- lapply(X = branch_partitions, function(x) {
          y <- cbind(partition = rep(NA, length(x = time_tree$edge.length)), rate = rep(NA, length(x = time_tree$edge.length)), changes = edge_changes, completeness = edge_completeness, duration = edge_durations)
          y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(z) rep(sum(y[z, "changes"]) / (sum(y[z, "completeness"] * y[z, "duration"])), length(x = z))))[order(unlist(x = x))]))
          y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
          y
        })
        
        # Get AIC, AICc and rate results:
        branch_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
      }
      
      # Pack branch partitions to test into single strings for output:
      packed_branch_partitions <- pack_partitions(branch_partitions)
      
      # Add packed partitions to results:
      for (i in 1:length(x = branch_test_results)) branch_test_results[[i]]$partition <- packed_branch_partitions[i]
      
      # If not performing branch partition tests:
    } else {
      
      # Make empty branch partition result output:
      branch_test_results <- NULL
    }
    
    # If performing clade partition tests:
    if (!is.null(clade_partitions)) {
      
      # Create clade rates for output:
      clade_rates <- do.call(what = rbind, args = lapply(X = lapply(X = as.list(x = n_tips + c(1:time_tree$Nnode)), function(x) find_descendant_edges(x, tree = time_tree)), function(y) c(changes = sum(edge_changes[y]), completeness = sum(edge_completeness[y] * edge_durations[y]), duration = 1)))
      
      # Add rates and node numbers:
      clade_rates <- cbind(node = n_tips + c(1:time_tree$Nnode), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = clade_rates[, "changes"] / clade_rates[, "completeness"])), clade_rates)
      
      # If using Likelihood Ratio Test:
      if (test_type == "lrt") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        partitioned_data <- lapply(X = clade_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))
        
        # Add sampled rate to paritioned data matrices:
        partitioned_data <- lapply(X = partitioned_data, function(x) {
          cbind(rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        })
        
        # Get LRT p-values and combine output as edge test results:
        clade_test_results <- lapply(X = partitioned_data, function(x) {
          list(rates = x[, "rate"], p_value = get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        })
      }
      
      # If using AIC:
      if (test_type == "aic") {
        
        # Build partitioned data for AIC:
        partitioned_data <- lapply(X = clade_partitions, function(x) {
          y <- cbind(partition = rep(NA, length(x = time_tree$edge.length)), rate = rep(NA, length(x = time_tree$edge.length)), changes = edge_changes, completeness = edge_completeness, duration = edge_durations)
          y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / (sum(y[x, "completeness"] * y[x, "duration"])), length(x = x))))[order(unlist(x = x))]))
          y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
          y
        })
        
        # Get AIC, AICc and rate results:
        clade_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
      }
      
      # Pack clade partitions to test into single strings for output:
      packed_clade_partitions <- pack_partitions(clade_partitions)
      
      # Add packed partitions to results:
      for (i in 1:length(x = clade_test_results)) clade_test_results[[i]]$partition <- packed_clade_partitions[i]
      
      
      # If not performing clade partition tests:
    } else {
      
      # Make empty clade partition result output:
      clade_test_results <- NULL
    }
    
    # If not doing clade OR branch tests:
  } else {
    
    # Make empty branch partition result output:
    branch_test_results <- NULL
    
    # Make empty clade partition result output:
    clade_test_results <- NULL
  }
  
  # If performing character partition tests:
  if (!is.null(character_partitions)) {
    
    # Get vector of (weighted) changes for each character:
    character_changes <- unlist(x = lapply(X = as.list(x = character_numbers), function(x) {
      character_rows <- which(x = inferred_character_changes[, "character"] == x)
      sum(inferred_character_changes[character_rows, "steps"] * inferred_character_changes[character_rows, "weight"])
    }))
    
    # Get vector of weighted durations for each character:
    character_durations <- (character_weights / sum(character_weights)) * sum(time_tree$edge.length)
    
    # Get vector of completness (opportunity to observe changes) for each character:
    character_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) {
      character_presence <- rep(0, times = length(x = character_numbers))
      character_presence[x$comparable_characters] <- 1
      character_presence * x$branch_duration
    })), 2, sum) / sum(time_tree$edge.length)
    
    # Set global rate:
    global_rate <- sum(character_changes) / sum(character_completeness * character_durations)
    
    # Create character rates for output:
    character_rates <- lapply(X = list(as.list(x = 1:max(character_numbers))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(character_changes[y]), sum(character_completeness[y] * character_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]
    
    # Add character numbers and rates:
    character_rates <- cbind(character = 1:max(character_numbers), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = character_rates[, "changes"] / character_rates[, "completeness"])), character_rates)
    
    # If using likelihood ratio test:
    if (test_type == "lrt") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      partitioned_data <- lapply(X = character_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(character_changes[y]), sum(character_completeness[y] * character_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))
      
      # Add sampled rate to paritioned data matrices:
      partitioned_data <- lapply(X = partitioned_data, function(x) {
        x <- cbind(as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        colnames(x = x)[1] <- "rate"
        x
      })
      
      # Get P-Values and combine output as edge test results:
      character_test_results <- lapply(X = partitioned_data, function(x) {
        x <- list(x[, "rate"], get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        names(x) <- c("rates", "p_value")
        x
      })
    }
    
    # If using AIC:
    if (test_type == "aic") {
      
      # Build partitioned data for AIC:
      partitioned_data <- lapply(X = character_partitions, function(x) {
        y <- cbind(partition = rep(NA, max(character_numbers)), rate = rep(NA, max(character_numbers)), changes = character_changes, completeness = character_completeness, duration = character_durations)
        y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / (sum(y[x, "completeness"] * y[x, "duration"])), length(x = x))))[order(unlist(x = x))]))
        y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
        y
      })
      
      # Get AIC, AICc and rate results:
      character_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
    }
    
    # Pack character partitions to test into single strings for output:
    packed_character_partitions <- pack_partitions(character_partitions)
    
    # Add packed partitions to results:
    for (i in 1:length(x = character_test_results)) character_test_results[[i]]$partition <- packed_character_partitions[i]
    
    # If performing branch partition tests:
  } else {
    
    # Make empty character partition result output:
    character_test_results <- NULL
  }
  
  # If performing time bin partition tests:
  if (!is.null(time_partitions)) {
    
    # Get weighted number of changes from each time bin:
    time_bin_changes <- unlist(x = lapply(X = as.list(x = 1:(length(x = time_bins) - 1)), function(x) {
      change_rows <- inferred_character_changes[, "bin"] == x
      sum(inferred_character_changes[change_rows, "steps"] * inferred_character_changes[change_rows, "weight"])
    }))
    
    # If using the Close time bin completeness approach get completeness value for each time bin:
    if (time_binning_approach == "close") time_bin_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$proportional_binned_edge_durations * (sum(character_weights[x$comparable_characters]) / sum(character_weights)))), 2, sum) / apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$proportional_binned_edge_durations)), 2, sum)
    
    # If using the Lloyd time bin completeness approach get completeness value for each time bin::
    if (time_binning_approach == "lloyd") time_bin_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) apply(matrix(character_weights[x$comparable_characters], ncol = 1) %*% x$binned_branch_durations, 2, sum))), 2, sum) / apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) apply(matrix(character_weights, ncol = 1) %*% x$binned_branch_durations, 2, sum))), 2, sum)
    
    # Get durations of edges in each time bin:
    time_bin_durations <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$binned_branch_durations)), 2, sum)
    
    # Set global rate (NB: will differ between Close and Lloyd approaches, but Lloyd approach will match edge or character global rate):
    global_rate <- sum(time_bin_changes) / sum(time_bin_completeness * time_bin_durations)
    
    # Create time rates for output:
    time_rates <- lapply(X = list(as.list(x = 1:(length(x = time_bins) - 1))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(time_bin_changes[y]), sum(time_bin_completeness[y] * time_bin_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]
    
    # Add bin numbers and rates:
    time_rates <- cbind(bin = 1:(length(x = time_bins) - 1), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = time_rates[, "changes"] / time_rates[, "completeness"])), time_rates)
    
    # If using Likelihood Ratio Test to compare partitions:
    if (test_type == "lrt") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later) for LRT:
      partitioned_data <- lapply(X = time_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(time_bin_changes[y]), sum(time_bin_completeness[y] * time_bin_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))
      
      # Add sampled rate to paritioned data matrices for LRT:
      partitioned_data <- lapply(X = partitioned_data, function(x) {
        x <- cbind(as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        colnames(x = x)[1] <- "rate"
        x
      })
      
      # Get P-Values and combine output as edge test results:
      time_test_results <- lapply(X = partitioned_data, function(x) {
        x <- list(x[, "rate"], get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        names(x) <- c("rates", "p_value")
        x
      })
    }
    
    # If using AIC to compare partitions:
    if (test_type == "aic") {
      
      # Build partitioned data for AIC:
      partitioned_data <- lapply(X = time_partitions, function(x) {
        y <- cbind(partition = rep(NA, length(x = time_bin_changes)), rate = rep(NA, length(x = time_bin_changes)), changes = time_bin_changes, completeness = time_bin_completeness * time_bin_durations, duration = rep(1, length(x = time_bin_changes)))
        y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / sum(y[x, "completeness"]), length(x = x))))))
        y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))
        y
      })
      
      # Get AIC, AICc and rate results:
      time_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
    }
    
    # Pack time bin partitions to test into single strings for output:
    packed_time_bin_partitions <- pack_partitions(time_partitions)
    
    # Add packed partitions to results:
    for (i in 1:length(x = time_test_results)) time_test_results[[i]]$partition <- packed_time_bin_partitions[i]
    
    # If not performing time bin partition tests:
  } else {
    
    # Make empty time bin partition result output:
    time_test_results <- NULL
  }
  
  # Set global rate for output:
  global_rate <- sum(unlist(x = lapply(X = edge_list, function(x) sum(x$character_changes[, "steps"] * x$character_changes[, "weight"])))) / sum(unlist(x = lapply(X = edge_list, function(x) sum(character_weights[x$comparable_characters]) / sum(character_weights))) * unlist(x = lapply(X = edge_list, function(x) x$branch_duration)))
  
  # If performing Likelihood Ratio Test:
  if (test_type == "lrt") {
    
    # Subfunction to calculate adjusted alphas for multiple comparison corrections:
    add_cutoffs <- function(test_results, alpha, multiple_comparison_correction = multiple_comparison_correction) {
      
      # Get number of comparisons performed:
      n_comparisons <- length(x = test_results)
      
      # If using the Benjamini-Hochberg false discovery rate approach:
      if (multiple_comparison_correction == "benjaminihochberg") {
        
        # Set cutoff values:
        cutoff_values <- ((1:n_comparisons) / n_comparisons) * alpha
        
        # Get actual p-values found:
        p_values <- unlist(x = lapply(X = test_results, "[[", "p_value"))
        
        # Order cutoffs by p-value rank:
        cutoff_values <- cutoff_values[rank(p_values, ties.method = "random")]
      }
      
      # If using the Bonferroni correction set cutoff values as alpha over N:
      if (multiple_comparison_correction == "bonferroni") cutoff_values <- alpha / n_comparisons
      
      # Add cutoffs to output:
      for (i in 1:length(x = test_results)) test_results[[i]]$CorrectedAlpha <- cutoff_values[i]
      
      # Return modified test results:
      test_results
    }
    
    # If doing branch partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(branch_partitions)) branch_test_results <- add_cutoffs(test_results = branch_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)
    
    # If doing character partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(character_partitions)) character_test_results <- add_cutoffs(test_results = character_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)
    
    # If doing clade partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(clade_partitions)) clade_test_results <- add_cutoffs(test_results = clade_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)
    
    # If doing time bin partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(time_partitions)) time_test_results <- add_cutoffs(test_results = time_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)
  }
  
  # Return compiled output:
  list(time_bins_used = time_bins, inferred_character_changes = inferred_character_changes, mean_rate = global_rate, continuous_characters_discretized = continuous_characters_discretized, branch_test_results = branch_test_results, character_test_results = character_test_results, clade_test_results = clade_test_results, time_test_results = time_test_results, branch_rates = branch_rates, character_rates = character_rates, clade_rates = clade_rates, time_rates = time_rates, time_tree = time_tree)
  }

## date_nodes() ################################################################
date_nodes <- function(time_tree) {
  
  # Need input checks
  
  # Get N tips:
  n_tips <- ape::Ntip(phy = time_tree)
  
  # Get N nodes:
  n_nodes <- ape::Nnode(phy = time_tree)
  
  # Store root node number:
  root_node <- n_tips + 1
  
  # If tree is a complete polytomy:
  if (time_tree$Nnode == 1) {
    
    # Create paths for just tips:
    paths <- as.list(x = 1:n_tips)
    
    # Add root to each path:
    for (i in 1:length(x = paths)) paths[[i]] <- c(paths[[i]], n_tips + 1)
    
    # If tree is not a complete polytomy:
  } else {
    
    # Create initial paths list with end nodes (terminal and internal, excluding the root):
    paths <- split(c(1:n_tips, (n_tips + 2):(n_tips + n_nodes)), f = 1:(n_tips + time_tree$Nnode - 1))
    
    # Strip names:
    names(paths) <- NULL
    
    # For each path:
    for (i in 1:length(x = paths)) {
      
      # Set counter as 1:
      j <- 1
      
      # Identify current node:
      current_node <- paths[[i]][j]
      
      # While current node is not the root (path has not terminated):
      while (current_node != root_node) {
        
        # Update current node and add to path:
        current_node <- paths[[i]][j + 1] <- time_tree$edge[match(current_node, time_tree$edge[, 2]), 1]
        
        # Update counter:
        j <- j + 1
      }
    }
  }
  
  # Create vector to store node ages:
  date_nodes <- vector(mode = "numeric", length = n_tips + time_tree$Nnode)
  
  # For each path:
  for (i in 1:length(x = paths)) {
    
    # Store path lengths from root:
    date_nodes[paths[[i]][1]] <- sum(time_tree$edge.length[match(paths[[i]][1:(length(x = paths[[i]]) - 1)], time_tree$edge[, 2])])
  }
  
  # Subtract path lengths from root time:
  date_nodes <- time_tree$root.time - date_nodes
  
  # Add node numbers:
  names(date_nodes) <- 1:(n_tips + time_tree$Nnode)
  
  # Return node ages:
  return(date_nodes)
}

## find_descendent_edges() #####################################################
find_descendant_edges <- function(n, tree) {
  
  # Find number of tips:
  n_tips <- ape::Ntip(phy = tree)
  
  # Find number of terminals (i.e. stopping point):
  n_terminals <- length(x = strap::FindDescendants(n = n, tree = tree))
  
  # Create vector to store internal nodes:
  nodes <- n
  
  # Create vector to store edge numbers (i.e. row numbers for tree$edge):
  edges <- grep(n, tree$edge[, 1])
  
  # Keep going until all descendant edges are found:
  while (length(x = which(x = tree$edge[edges, 2] <= n_tips)) < n_terminals) {
    
    # Get internal nodes found so far:
    nodes <- tree$edge[edges, 2][which(x = tree$edge[edges, 2] > n_tips)]
    
    # For each node add any new descendant edges:
    for (i in nodes) edges <- sort(x = unique(x = c(edges, which(x = tree$edge[, 1] == i))))
  }
  
  # Return edges vector:
  edges
}

## prune_cladistic_matrix() ####################################################
prune_cladistic_matrix <- function(cladistic_matrix, blocks2prune = c(), characters2prune = c(), taxa2prune = c(), remove_invariant = FALSE) {
  
  # How do blocks and characters to prune interact? (explain to user in manual)
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  # FOLLOWNIG MODIFIED TO SILENCE CHECK
  # if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Subfunction to find length of character types for each character (i.e., unique values excluding polymorphisms but included inapplicables):
  find_length <- function(x) {
    
    # Convert each column of matrix to a list of numeric values:
    x <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = apply(x, 2, as.list), unlist), strsplit, split = "&|/"), unlist), unique), sort), length)
    
    # Return(x):
    return(x)
  }
  
  # Check that something to prune has been specified:
  if (is.null(blocks2prune) && is.null(characters2prune) && is.null(taxa2prune) && remove_invariant == FALSE) stop("No blocks, taxa, or characters to prune specified.")
  
  # Check blocks specified exist and stop and warn :
  if (length(x = setdiff(x = blocks2prune, y = 1:(length(x = cladistic_matrix) - 1))) > 0) stop("Block numbers specified that are not present in data.")
  
  # Check characters specified exist and stop and warn user if not:
  if (length(x = setdiff(x = characters2prune, y = 1:sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol))))) > 0) stop("characters specified that are outside the scope of the matrix. Check and retry.")
  
  # Check taxa specified exist and stop and warn user if not:
  if (length(x = setdiff(x = taxa2prune, y = rownames(x = cladistic_matrix$matrix_1$matrix))) > 0) stop("Taxa specified that are not found in the matrix. Check and retry.")
  
  # Get number of characters:
  n_characters <- sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))
  
  # If there are characters to prune:
  if (!is.null(characters2prune)) {
    
    # Get character blocks for each character in descendng order (as want to work backwards so things match up properly):
    character_blocks <- unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = characters2prune, decreasing = TRUE)), ">", cumsum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))), which), length)) + 1
    
    # Initial build of characters in list form:
    characters_as_list <- lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), function(x) 1:ncol(x))
    
    # Actually form list of character numbers (i.e., renumber characters in second or higher blocks):
    if (length(x = characters_as_list) > 1) for (i in 2:length(x = characters_as_list)) characters_as_list[[i]] <- characters_as_list[[i]] + max(characters_as_list[[(i - 1)]])
    
    # For each unique character block:
    for (i in unique(x = character_blocks)) {
      
      # Find columns to delete in ith matrix:
      columns_to_delete <- match(sort(x = characters2prune, decreasing = TRUE)[character_blocks == i], characters_as_list[[i]])
      
      # Remove characters from matrix:
      cladistic_matrix[[(i + 1)]]$matrix <- cladistic_matrix[[(i + 1)]]$matrix[, -columns_to_delete, drop = FALSE]
      
      # Remove characters from ordering:
      cladistic_matrix[[(i + 1)]]$ordering <- cladistic_matrix[[(i + 1)]]$ordering[-columns_to_delete]
      
      # Remove characters from weights:
      cladistic_matrix[[(i + 1)]]$character_weights <- cladistic_matrix[[(i + 1)]]$character_weights[-columns_to_delete]
      
      # Remove characters from minimum values:
      cladistic_matrix[[(i + 1)]]$minimum_values <- cladistic_matrix[[(i + 1)]]$minimum_values[-columns_to_delete]
      
      # Remove characters from maximum values:
      cladistic_matrix[[(i + 1)]]$maximum_values <- cladistic_matrix[[(i + 1)]]$maximum_values[-columns_to_delete]
    }
  }
  
  # If there are taxa to prune:
  if (!is.null(taxa2prune)) {
    
    # Remove pruned taxa from each block:
    for (i in 2:length(x = cladistic_matrix)) cladistic_matrix[[i]]$matrix <- cladistic_matrix[[i]]$matrix[-match(taxa2prune, rownames(x = cladistic_matrix[[i]]$matrix)), , drop = FALSE]
  }
  
  # If there are blocks to prune:
  if (!is.null(blocks2prune)) {
    
    # Remove blocks to be rpuned:
    cladistic_matrix <- cladistic_matrix[-(blocks2prune + 1)]
    
    # Rename (renumber) remaining matrix blocks:
    names(cladistic_matrix[2:length(x = cladistic_matrix)]) <- paste("matrix_", 1:(length(x = cladistic_matrix) - 1), sep = "")
  }
  
  # If there are invariant characters:
  if (remove_invariant) {
    
    # Find any invariant characters:
    invariants_as_list <- lapply(X = lapply(X = lapply(X = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), find_length), unlist), "<", 1), which)
    
    # If there are invariant characters:
    if (length(x = unlist(x = invariants_as_list)) > 0) {
      
      # For each matrix block:
      for (i in 1:length(x = invariants_as_list)) {
        
        # Only if there are invraints for this block:
        if (length(x = invariants_as_list[[i]]) > 0) {
          
          # Remove characters from matrix:
          cladistic_matrix[[(i + 1)]]$matrix <- cladistic_matrix[[(i + 1)]]$matrix[, -invariants_as_list[[i]], drop = FALSE]
          
          # Remove characters from ordering:
          cladistic_matrix[[(i + 1)]]$ordering <- cladistic_matrix[[(i + 1)]]$ordering[-invariants_as_list[[i]]]
          
          # Remove characters from weights:
          cladistic_matrix[[(i + 1)]]$character_weights <- cladistic_matrix[[(i + 1)]]$character_weights[-invariants_as_list[[i]]]
          
          # Remove characters from minimum values:
          cladistic_matrix[[(i + 1)]]$minimum_values <- cladistic_matrix[[(i + 1)]]$minimum_values[-invariants_as_list[[i]]]
          
          # Remove characters from maximum values:
          cladistic_matrix[[(i + 1)]]$maximum_values <- cladistic_matrix[[(i + 1)]]$maximum_values[-invariants_as_list[[i]]]
        }
      }
    }
  }
  
  # Check for empty blocks and store them as blocks to delete if found:
  new_blocks_to_delete <- which(x = unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)) == 0)
  
  # If there are new blocks to prune:
  if (length(x = new_blocks_to_delete) > 0) {
    
    # Remove blocks to be rpuned:
    cladistic_matrix <- cladistic_matrix[-(new_blocks_to_delete + 1)]
  }
  
  # Rename (renumber) matrix blocks to ensure consistent output:
  names(cladistic_matrix[2:length(x = cladistic_matrix)]) <- paste("matrix_", 1:(length(x = cladistic_matrix) - 1), sep = "")
  
  # Ensure class is set:
  class(cladistic_matrix) <- "cladisticMatrix"
  
  # Return pruned matrix:
  cladistic_matrix
}

## plot_rates_character() ######################################################
plot_rates_character <- function(test_rates_output, model_number, ...) {
  
  # TO DO:
  #
  # - Add more options (plot colours, whetehr to show character numbers etc.)
  # - Add checks that data are present (characters were tested for)
  
  # Reconstruct character partitions list:
  character_partitions <- lapply(X = test_rates_output$character_test_results, function(x) {
    lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
      unlist(x = lapply(X = y, function(z) {
        z <- as.list(x = strsplit(z, split = " ")[[1]])
        unlist(x = lapply(X = z, function(p) {
          if (length(x = grep("-", p)) > 0) {
            p <- strsplit(p, split = "-")[[1]]
            as.numeric(p[1]:as.numeric(p[2]))
          } else {
            as.numeric(p)
          }
        }))
      }))
    })
  })
  
  # Get x values for plotting partitions of model:
  model_x_values <- lapply(X = apply(matrix(c(1, cumsum(unlist(x = lapply(X = character_partitions[[model_number]], function(x) range(1:length(x = x))))[-1])), ncol = 2), 2, list), unlist)
  
  # Get y values for plotting partitions of model:
  model_y_values <- lapply(X = as.list(x = test_rates_output$character_test_results[[model_number]]$rates), rep, 2)
  
  # Make vector of partition colours ready for plotting:
  partition_colours <- unlist(x = unname(mapply(function(x, y) {
    rep(x, length(x = y))
  }, x = hcl.colors(n = length(x = character_partitions[[model_number]]), alpha = 0.5, palette = "viridis"), y = character_partitions[[model_number]])))
  
  # Plot character rates:
  graphics::plot(x = 1:max(unlist(x = character_partitions[[model_number]])), y = test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"], pch = 21, bg = partition_colours, cex = 1.5, xlab = "Character", ylab = "Changes per lineage myr", ylim = c(0, 1.1 * max(test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"])), xaxt = "n", lwd = 0.5, col = "black")
  
  # Add character numbrs to plot:
  graphics::text(x = 1:max(unlist(x = character_partitions[[model_number]])), y = test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"], label = unlist(x = character_partitions[[model_number]]), pos = 1, col = partition_colours, cex = 0.5)
  
  # Add lines representing clustering of requested model to plot:
  for (i in 1:length(x = model_x_values)) graphics::lines(x = c(model_x_values[[i]][1] - 0.5, model_x_values[[i]][2] + 0.5), y = model_y_values[[i]])
  
  # Add model parameters (lambda values) to plot:
  for (i in 1:length(x = model_x_values)) graphics::text(x = mean(model_x_values[[i]]), y = mean(model_y_values[[i]]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3, cex = 1.5)
  
  # Add legend to plot:
  graphics::legend(x = 0, y = 1.1 * max(test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"]), legend = paste0("Partition ", 1:length(x = character_partitions[[model_number]])), pch = rep(21, length(x = character_partitions[[model_number]])), pt.bg = unique(x = partition_colours), col = rep("black", length(x = character_partitions[[model_number]])), pt.lwd = 0.5, pt.cex = 1.5)
}


## plot_rates_tree() ###########################################################
plot_rates_tree <- function(test_rates_output, model_type, model_number, ...) {
  
  # TO DO:
  #
  # - Check model_type is a valid choice
  # - Maybe work out how to plot time trees this way too
  # - Add legend for partitions somehow?
  # - Work out how to drop/edit subtitle from legend
  # - Work out how to not round rates to 1dp (0.005 -> 0 which is useless)
  # - Make example run faster
  
  # Check model type is a valid option:
  if (!model_type %in% c("branch", "clade")) stop("model_type must be one of \"branch\" or \"clade\".")
  
  # Set resolution for plotting (discretisation of continuous rates):
  resolution <- 100
  
  # If requesting branch partitions extract these from rate output:
  if (model_type == "branch") {
    edge_partitions <- lapply(X = test_rates_output$branch_test_results, function(x) {
      lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
        unlist(x = lapply(X = y, function(z) {
          z <- as.list(x = strsplit(z, split = " ")[[1]])
          unlist(x = lapply(X = z, function(p) {
            if (length(x = grep("-", p)) > 0) {
              p <- strsplit(p, split = "-")[[1]]
              as.numeric(p[1]:as.numeric(p[2]))
            } else {
              as.numeric(p)
            }
          }))
        }))
      })
    })
  }
  
  # If requesting clade partitions extract these from rate output:
  if (model_type == "clade") {
    edge_partitions <- lapply(X = test_rates_output$clade_test_results, function(x) {
      lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
        unlist(x = lapply(X = y, function(z) {
          z <- as.list(x = strsplit(z, split = " ")[[1]])
          unlist(x = lapply(X = z, function(p) {
            if (length(x = grep("-", p)) > 0) {
              p <- strsplit(p, split = "-")[[1]]
              as.numeric(p[1]:as.numeric(p[2]))
            } else {
              as.numeric(p)
            }
          }))
        }))
      })
    })
  }
  
  # If requesting branch rates extract these from output:
  if (model_type == "branch") edge_rates <- lapply(X = test_rates_output$branch_test_results, function(x) x$rates)
  
  # If requesting clade rates extract these from output:
  if (model_type == "clade") edge_rates <- lapply(X = test_rates_output$clade_test_results, function(x) x$rates)
  
  # Get discretized vector of edge rates (needed for choosing plot colours):
  discretized_rate_values <- seq(from = 0, to = max(edge_rates[[model_number]]), length.out = resolution)
  
  # Discretize edge rates:
  discretized_edge_rates <- lapply(X = edge_rates, function(x) unlist(x = lapply(X = as.list(x = x), function(y) discretized_rate_values[max(which(x = y >= discretized_rate_values))])))
  
  # Create vector of edge rate values to use in plotting:
  edge_rate_values <- rep(0, nrow(test_rates_output$time_tree$edge))
  
  # Fill vector of edge rate values to use in plotting:
  for (i in 1:length(x = edge_partitions[[model_number]])) edge_rate_values[edge_partitions[[model_number]][[i]]] <- discretized_rate_values[discretized_rate_values == discretized_edge_rates[[model_number]][i]]
  
  # Plot tree with branches colour coded by rate:
  phytools::plotBranchbyTrait(tree = test_rates_output$time_tree, x = edge_rate_values, mode = "edge", xlims = c(0, max(edge_rate_values)), title = "Changes per lineage myr", leg = max(nodeHeights(test_rates_output$time_tree)), ...)
  
  # Dead code attempting to basically do what phytools::plotBranchbyTrait does without calling phytools:
  # names(discretized_rate_values) <- hcl.colors(n = resolution, palette = "viridis")
  # EdgeColours <- rep("white", nrow(test_rates_output$time_tree$edge))
  # for(i in 1:length(x = edge_partitions[[model_number]])) EdgeColours[edge_partitions[[model_number]][[i]]] <- names(discretized_rate_values[discretized_rate_values == discretized_edge_rates[[model_number]][i]])
  # ape::plot.phylo(x = test_rates_output$time_tree, edge.color = EdgeColours, show.tip.label = FALSE, edge.width = 3)
  # cols = names(discretized_rate_values)
  # tree = test_rates_output$time_tree
  # lwd = 4
  # lims = c(0, max(discretized_rate_values))
  # Rounder <- (-1 * (min(c(1, ceiling(log(lims[2], base = 10)))) - 1) + 1)
  # leg <- round(0.5 * max(nodeHeights(tree)), 2)
  # x <- max(nodeHeights(tree)) / 2
  # y <- 10
  # fsize <- 1.0
  # X <- x + cbind(0:(length(x = cols) - 1) / length(x = cols), 1:length(x = cols) / length(x = cols)) * (leg)
  # Y <- cbind(rep(y, length(x = cols)), rep(y, length(x = cols)))
  # lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 2]), lwd = lwd + 2, lend = 2)
  # for(i in 1:length(x = cols)) lines(X[i, ], Y[i, ], col = cols[i], lwd = lwd, lend = 2)
  # text(x = x, y = y, "0", pos = 3, cex = fsize)
  # text(x= x + leg, y = y, round(lims[2], Rounder), pos = 3, cex = fsize)
  # text(x = (2 * x + leg) / 2, y = y, "Changes", pos = 3, cex = fsize)
  # text(x = (2 * x + leg) / 2, y = y, "per lineage myr", pos = 1, cex = fsize)
}

## plot_rates_time() ###########################################################
plot_rates_time <- function(test_rates_output, model_number, ...) {
  
  # TO DO:
  #
  # - Make points size of evidence (i.e., amount of information) with legend.
  # - Add better example that runs.
  
  # Build vector of time bin midpoints for plotting:
  time_bin_midpoints <- (test_rates_output$time_bins_used[2:length(x = test_rates_output$time_bins_used)] + test_rates_output$time_bins_used[1:(length(x = test_rates_output$time_bins_used) - 1)]) / 2
  
  # Get partitions used from results output:
  time_bin_partitions <- lapply(X = test_rates_output$time_test_results, function(x) {
    lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
      if (length(x = grep("-", y)) > 0) {
        z <- strsplit(y, split = "-")[[1]]
        y <- paste0(z[1]:z[2])
      }
      as.numeric(y)
    })
  })
  
  # Get sampled rates for model:
  time_rates <- cbind(lapply(X = time_bin_partitions[model_number], function(x) {
    do.call(what = rbind, args = lapply(X = x, function(y) {
      xs <- c(test_rates_output$time_bins_used[y[1]], test_rates_output$time_bins_used[(y[length(x = y)] + 1)])
    }))
  })[[1]], test_rates_output$time_test_results[[model_number]]$rates, test_rates_output$time_test_results[[model_number]]$rates)
  
  # Create base plot of rates in each time bin with any other requested options paseed as ...:
  geoscale::geoscalePlot(ages = time_bin_midpoints, data = test_rates_output$time_rates[, "rate"], age.lim = c(max(test_rates_output$time_bins_used), min(test_rates_output$time_bins_used)), data.lim = c(0, max(test_rates_output$time_rates[, "rate"]) * 1.1), pch = 20, cex.pt = 2, label = "Character changes per lineage million years", ...)
  
  # Add lines representing clustering of requested model to plot:
  for (i in 1:nrow(time_rates)) lines(x = time_rates[i, 1:2], y = time_rates[i, 3:4])
  
  # Add model parameters (lambda values) to plot:
  for (i in 1:nrow(time_rates)) text(x = mean(as.numeric(time_rates[i, 1:2])), y = as.numeric(time_rates[i, 3]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3)
}


## partition_time_bins() ###########################################################
partition_time_bins <- function(n_time_bins, partition_sizes_to_include = "all") {
  
  # Build time bins vector:
  time_bins <- 1:n_time_bins
  
  # Check there are multiple time bins:
  if (length(x = time_bins) < 2) stop("There must be at least two time bins.")
  
  # Format partition_sizes_to_include as a vector of all possible numbers if input is "all":
  if (any(partition_sizes_to_include == "all")) partition_sizes_to_include <- 1:n_time_bins
  
  # Get number of possible "switches", i.e., positions where a partition can (1) or cannot (0) be:
  n_switches <- length(x = time_bins) - 1
  
  # Calculate expeted number of partitions (if large can stop and warn user):
  n_partitions <- sum(unlist(x = lapply(X = as.list(x = partition_sizes_to_include - 1), function(x) ncol(combn(n_switches, x)))))
  
  # Generate starting splitsiwtches vector:
  split_switches <- as.character(0:1)
  
  # Generate all possible combinations of switches iteratively:
  while (nchar(x = split_switches[1]) < (length(x = time_bins) - 1)) split_switches <- apply(expand.grid(split_switches, as.character(0:1)), 1, paste, collapse = "")
  
  # Work out partition sizes for subsetting below:
  partition_sizes <- unlist(x = lapply(X = strsplit(split_switches, split = ""), function(x) sum(as.numeric(x)))) + 1
  
  # Collpase switches vector to just those of the reequired input size:
  split_switches <- split_switches[!is.na(match(partition_sizes, partition_sizes_to_include))]
  
  # Subfunction to generate partition positions:
  set_partition_positions <- function(switch_sequence) {
    
    # How long should the output vector be?:
    vector_length <- nchar(x = switch_sequence) + 1
    
    # Turn siwtches into split after points (adding end in to complete sequence):
    split_after <- c(which(x = as.numeric(strsplit(switch_sequence, split = "")[[1]]) == 1), vector_length)
    
    # Seed start position:
    start_position <- 1
    
    # Set empty output list:
    output <- list()
    
    # For each split:
    for (i in split_after) {
      
      # Add to output:
      output[[(length(x = output) + 1)]] <- start_position:i
      
      # Update start position:
      start_position <- i + 1
    }
    
    # Return list of vectors:
    output
  }
  
  # Convert split switch strings to lists of vectors of partition elements and return:
  lapply(X = as.list(x = split_switches), function(x) set_partition_positions(x))
}
