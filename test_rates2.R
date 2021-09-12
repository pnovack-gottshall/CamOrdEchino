## Modified version of test_rates() function from 'Claddis' package, Modified to
## accept pre-inferred ancestral states output from
## Claddis::estimate_ancestral_states(). Code by Graeme Lloyd and taken from
## https://github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R ,
## committed 4/19/2021. Thanks, Graeme!

## The subfunction used in calculating AIC, calculate_partition_AIC(), is also
## modified to move the calculation of the log into the dpois() function. This
## allows it to correctly calculate very small density values that otherwise
## yield incalculable AIC values.

## Also includes the code for is.timeBins() and check_timeBins, which are on
## GitHub (committed Oct. 11, 2020 and April 4, 2021, respectively), but not
## apparently included in Claddis 0.6.3 on CRAN.

## All changes are tagged with *** in code

## MODIFIED test_rates() #######################################################
test_rates2 <- function(time_tree,
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
  
  # MAKE RATES OUTPUT MAKE SENSE BY PROPERLY SEPARATING COMPLETENESS AND DURATION
  
  # ADD EXAMPLES OF VISUALISED OUTPUT (MENTION IN MANUAL AN ADD PLOT FUNCTIONS TO SEEALSO)
  
  # TO EXCLUDE OUTGROUP CAN JUST SET UP PARTITIONS THAT EXCLUDE THESE (REQUIRES REMOVING PARTITION FIXER AS STANDARD)
  # COULD ALSO DO INTERNAL AND TERMINAL BRANCHES SEPARATELY THIS WAY
  # ALLOW NAMING OF PARTITIONS AND POTENTIALLY CALLING NAME FOR PLOTS
  
  # AICc BREAKS IF MORE THAN N-2 PARAMETERS (INFINITY OR WRONGLY NEGATIVE OUTPUT CAN OCCUR THIS WAY)
  # GLOBAL RATE NEEDS TO GO INTO LRT OPTION IF NOT USED IN AIC
  # SOMEHOW NEED TO ALLOW MULTIPLE VERSIONS OF MAIN PIEPLINE IF DOING RANDOM ASSIGNMENTS OF CHANGE TIMES
  # NEED TO ADD PARTITIONS USED TO OUTPUT SOMEHOW...
  # STOP REQUIRING TIME BINS TO BE SET IF NOT TESTING TIME PARTITIONS.
  
  # MAYBE DO SOME KIND OF WEIGHTING FOR AIC AS SOME PARTITIONS WILL CONTAIN VERY LITTLE DATA AND BE EXTREME OUTLIERS?
  # NEED TO CHECK FOR SINGLE PARTITION WITH LRT (ALLOWED WITH AIC)
  # NEED TO CHECK FOR FULLY BIFURCATING IF IMPLEMENTING WANG STUDENTS APPROACH? OR IS THAT DONE ANYWAY?
  # IF USING AIC NEED TO CHECK FOR EACH TEST TYPE AT LEAST TWO PARTITIONS ARE SUPPLIED OR ELSE THE RESULTS ARE MEANINGLESS
  # ADD CHECK TO INCLUDE ALL LESS COMPLEX SUBPARTITIONS OF ANY 3 OR MORE SIZE PARTITOINING OF THE DATA
  
  # DESIDERATA (STUFF IT WOULD BE NICE TO ADD IN FUTURE):
  #
  # WRITE SEARCH VERSION FOR FINDING RATE SHIFTS? SHOULD THIS EVEN BE AN OPTION? DOES THIS REQUIRE MODIFYING LRT TO COMPARE E.G. 2-RATE DIRECTLY WITH 3-RATE MODEL? WOULD NEED TO PERMUTE ALL POSSIBLE COMBOS AND NOT SURE HOW LARGE THESE MIGHT GET (VERY FOR EDGES).
  # MAYBE MAKE ANCESTRAL STATE UNCERTAINTY DIFFERENT FOR TIPS THAN NODES? I.E., HOW IT GETS RESOLVED CAN BE DIFFERENT (MORE OPTIONS TO FUNCTION)
  # THESE TWO ARE RELATED: 1. ADD TERMINAL VERSUS INTERNAL OPTION SOMEHOW/SOMEWHERE (BRANCH TYPE ALREADY RECORDED ON EDGE LIST!), 2. ALLOW OPTION TO IGNORE SOME PARTS OF THE TREE FOR PARTITION TESTS? MAKES CALCULATING THE MEAN RATE TRICKIER BUT MIGHT MAKE SENSE E.G. FOR INGROUP ONLY TESTS. EXCLUDE EDGES AFTER DOING ANCESTRAL STATES? OR SET THESE TO ALL NAS TO ENSURE OTHER THINGS WORK FINE? FOR EXAMPLE, USE OUTGROUPS TO SET ANCESTOR THEN EXCLUDE THEM FROM THE PROCESS.
  # EXTRA FUNCTION(S) TO VISUALISE RESULTS MOST LIKELY. DEFO NEEDED! HEAT MAP WITH EDGE BLOCKS?
  # CHECK FOR AUTAPOMORPHIES AND INFORM USER IF FOUND?
  # ADD CONTRIVED EXAMPLES (UNIT TESTS) TO SHOW HOW FUNCTION WORKS, E.G. RATE OF ONE CHANGES PER MILLION YEARS THEN DUPLICATED BLOCK WITH CHARACTER PARTITION TEST.
  # PROBABLY NEED MORE CAREFUL CHECKS FOR ZERO VALUES GENERALLY, E.G., CHARACTER WITH ALL MISSING DATA
  # ALLOW REWEIGHTING OF INAPPLICABLES ZERO AS AN OPTION FOR THEM?
  # HOW TO FORMAT OUTPUT? GET CIS FOR EACH PARTITION FOR VISUALISATION (E.G., BARPLOT OF PARTITION VALUES WITH DASHED LINE FOR MEAN AND ERROR BARS FOR CIS)? STICK WITH LIST OR COLLAPSE TO A MATRIX SOMEHOW?
  # TIME BINS WITH NOTHING IN WILL CAUSE ISSUES AS DIVIDE BY ZEROES WILL OCCUR - ADD CHECK FOR THIS.
  # WHAT IS SIGNIFICANTLY HIGH OR LOW IF THERE ARE THREE OR MORE PARTITIONS? THIS IS NOT EVEN IN OUTPUT YET. PROLLY CANNOT DO FULL STOP NOW PARTITIONS ARE MORE COMPLEX
  
  # *** MODIFIED TO FORCE SETTING CLASS OF cladistic_matrix ***
  class(cladistic_matrix) <- "cladisticMatrix"
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
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
  
  # Check time_bins is in a valid format and stop and warn user if not:
  if (!is.timeBins(x = time_bins)) stop(check_timeBins(time_bins = time_bins))
  
  # Find the Time bin midpoints:
  time_bin_midpoints <- find_time_bin_midpoints(time_bins = time_bins)
  
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
  
    # *** Modified to handle very small dpois values. (In occasional instances,
    # a character that changes at a very high rate will yield a nonsensical
    # density value of 0 due to rounding errors, which results in incalculable
    # AIC values.) ***
  
    # Subfunction to calculate AIC from partition (with columns labelled partition, rate, completeness, duration):
  calculate_partition_aic <- function(partition, aicc = FALSE) {
    
    # *** NEXT LINE MODIFIED FROM ORIGINAL ***
    log_mle <- sum(dpois(round(partition[, "changes"]), partition[, "rate"] * partition[, "completeness"] * partition[, "duration"], log = TRUE))
    
    # Get log maximum likelihood estimate:
    # log_mle <- sum(log(dpois(round(partition[, "changes"]), partition[, "rate"] * partition[, "completeness"] * partition[, "duration"])))

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
  # *** MODIFIED TO ACCEPT PREVIOUSLY ESTIMATED ANCESTRAL MATRIX ***
  ancestral_states <- cladistic_matrix
  
  # Build single matrix of all states in tip label then node number order:
  # *** MODIFIED AS WORKAROUND FOR ERROR IN NEXT LINE ***
  # all_states <- do.call(what = cbind, args = lapply(X = lapply(X = ancestral_states[2:length(x = ancestral_states)], "[[", "matrix"), function(x) x[c(time_tree$tip.label, 1:n_nodes + n_tips), , drop = FALSE]))
  all_states <- cbind(ancestral_states$matrix_1$matrix, ancestral_states$matrix_2$matrix)
  
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
  if (any(ordering == "continuous")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Set default converted continuous characters to TRUE:
    continuous_characters_discretized <- TRUE
    
    # Find out which characters are continuous:
    continuous_characters_found <- which(x = ordering == "continuous")
    
    # Rescale continous characters as zero to one values:
    list_of_continuous_values_rescaled_zero_to_one <- lapply(X = lapply(X = lapply(X = apply(all_states[, continuous_characters_found, drop = FALSE], 2, list), unlist), as.numeric), function(x) {
      x <- x - min(sort(x = x))
      x <- x / max(sort(x = x))
      return(x)
    })
    
    # Now discretize and store these characters (0 to 31 scale):
    all_states[, continuous_characters_found] <- do.call(what = cbind, args = lapply(X = lapply(X = lapply(X = list_of_continuous_values_rescaled_zero_to_one, function(x) as.list(x = x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x = x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    ordering[continuous_characters_found] <- "ordered"
    
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
    if (length(x = character_differences) > 0) character_changes <- rbind(character_changes, cbind(as.numeric(comparable_characters[character_differences]), as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]), as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]), ifelse(comparable_ordering[character_differences] == "unordered", 1, abs(as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]) - as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]))), comparable_weights[character_differences]))
    
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
      character_changes <- cbind(character_changes, unlist(x = lapply(X = as.list(x = character_changes[, "time"]), function(x) max(which(x = x <= time_bins[, "fad"])))))
      
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
    boundaries_crossed <- unname(obj = time_bins[, "lad"][intersect(which(x = time_bins[, "lad"] > lad), which(x = time_bins[, "lad"] < fad))])
    
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
    colnames(x = branch_sections) <- unlist(x = lapply(X = as.list(x = branch_sections["fad", ]), function(x) max(which(x = x <= time_bins[, "fad"]))))
    
    # Add new list section for branch (edge) sections binned by time:
    x$binned_edge_sections <- branch_sections
    
    # Return output:
    x
  }
  
  # Get edge sections in time bins:
  edge_list <- lapply(X = edge_list, get_edge_sections_in_bins, time_bins = time_bins)
  
  # Add binned branch durations to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    branch_durations <- rep(0, nrow(x = time_bins))
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
    time_bin_changes <- unlist(x = lapply(X = as.list(x = 1:nrow(x = time_bins)), function(x) {
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
    time_rates <- lapply(X = list(as.list(x = 1:nrow(x = time_bins))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(time_bin_changes[y]), sum(time_bin_completeness[y] * time_bin_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]
    
    # Add bin numbers and rates:
    time_rates <- cbind(bin = 1:nrow(x = time_bins), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = time_rates[, "changes"] / time_rates[, "completeness"])), time_rates)
    
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
    add_cutoffs <- function(test_results, alpha, multiple_comparison_correction) {
      
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

##### FUNCTIONS NOT INCLUDED ON CRAN VERSION OF CLADDIS 0.6.3 ###################
is.timeBins <- function(x) {
  
  # Get any error messages for time_bins:
  messages <- check_timeBins(time_bins = x)
  
  # Return logical indicating whether object is a valid timeBins object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}

check_timeBins <- function(time_bins) {
  
  # Check time_bins has class timeBins and add error message to output if true:
  if (!inherits(x = time_bins, what = "timeBins")) return("time_bins must be an object of class \"timeBins\".")
  
  # Check time bins are in form of matrix add error message to output if false:
  if (!is.matrix(x = time_bins)) return("time_bins must be in the form of a matrix.")
  
  # Check time_bins has two columns and add error message to output if false:
  if (ncol(x = time_bins) != 2) return("time_bins must have exactly two columns.")
  
  # Check time_bins has mutliple rows (bins) and add error message to output if false:
  if (nrow(x = time_bins) < 2) return("time_bins must have at least two rows.")
  
  # Check time_bins column names are null and add error message to output if true:
  if (is.null(x = colnames(x = time_bins))) return("time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names are null and add error message to output if true:
  if (is.null(x = rownames(x = time_bins))) return("time_bins must have unique row names corresponding to each time bin.")
  
  # Check time_bins column names are correct add error message to output if false:
  if (!all(colnames(x = time_bins) == c("fad", "lad"))) return("time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names have at least one character each and add error message to output if false:
  if (!all(x = nchar(x = rownames(x = time_bins)) > 0)) return("time_bins must have row names formed from at least one character.")
  
  # Check time_bins row names are all unique and add error message to output if false:
  if (any(x = duplicated(x = rownames(x = time_bins)))) return("time_bins must have unique row names.")
  
  # Check time_bins contains only numeric values and add error message to output if false:
  if (!is.numeric(x = time_bins)) return("time_bins must contain numeric values representing millions of years ago (Ma).")
  
  # Check time_bins start before they end and have positive length and add error message to output if false:
  if (!all(x = time_bins[, "fad"] > time_bins[, "lad"])) return("time_bins must consist of \"fad\" values that exceed corresponding \"lad\" values. I.e., each time bin must have an older \"fad\" than \"lad\" and be of positive length.")
  
  # Check time_bins do not overlap each other and add error message to output if true:
  if (any(x = unlist(x = lapply(X = apply(X = time_bins, MARGIN = 1, FUN = list), FUN = function(x) length(x = which(x = apply(X = cbind(x[[1]]["fad"] > time_bins[, "lad"], x[[1]]["lad"] < time_bins[, "fad"]), MARGIN = 1, FUN = all))))) > 1)) return("time_bins must not overlap each other.")
  
  # Check time_bins are in correct order and add error message to output if false:
  if (!all(x = order(x = time_bins[, "fad"]) == seq(from = nrow(x = time_bins), to = 1, by = -1))) return("time_bins must be ordered from oldest (top) to youngest (bottom).")
  
  # Check time_bins have no gaps between them and add error message to output if true:
  if (!all(x = time_bins[1:(nrow(x = time_bins) - 1), "lad"] == time_bins[2:nrow(x = time_bins), "fad"])) return("time_bins must abut each other. I.e., there can be no gaps between successive time bins.")
  
  # Return empty vector:
  vector(mode = "character")
}

find_time_bin_midpoints <- function(time_bins) {
  
  # Check time_bins has class timeBins and stop and warn user if not:
  if (!inherits(x = time_bins, what = "timeBins")) stop("time_bins must be an object of class \"timeBins\".")
  
  # If not a valid timeBins object then stop and provide feedback to user on what is wrong:
  if (!is.timeBins(time_bins)) stop(check_timeBins(time_bins = time_bins)[1])
  
  # Return time bin midpoints:
  apply(X = time_bins, MARGIN = 1, FUN = mean)
}