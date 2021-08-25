## USING CLADDIS TO INFER ANCESTRAL STATES #####################################

## Prior to running, run 1-MakeTimetime_trees.R to build time-calibrated trees
## using 'paleotree' package.

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 7/29/2021, v.
# 0.6.3. Users accustomed with earlier versions will need to either download the
# archived version of the package from GitHub or alter the code accordingly to
# use appropriate argument and function names.

## PREPARATIONS ################################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")
setwd("~/Manuscripts/CamOrdEchinos/Data files/NA reformatted/")

# Load packages
library(beepr)      # v. 1.3
library(paleotree)  # v. 3.3.25
library(parallel)   # v. 4.1.0
library(snowfall)   # v. 1.84-6.1
library(Rcpp)       # v. 1.0.6
library(phytools)    # v. 0.7-83
if(packageVersion("phytools") < "0.7-54")
  stop("wrong version of 'phytools' Get updated version from GitHub.\n")
library(paleotree)  # v. 3.3.25
library(Claddis)    # v. 0.6.3 - Check SI for Lloyd 2018 for walk-through on code functions
if(packageVersion("Claddis") < "0.6.0")
  stop("wrong version of 'Claddis' Get updated version from GitHub.\n")


## IMPORT AND PROCESS TIME TREES ###############################################

## Import time trees saved from prior code
load("~/Manuscripts/CamOrdEchinos/cal3trees")

# Plot median diversity curve with 95%iles
cal3_multiDiv <- paleotree::multiDiv(cal3trees, plot = FALSE)
paleotree::plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))

# Address issue of zero-length branches (ZLBs) (because Claddis' morphological
# disparity functions disallow it).

# Because Claddis disallows using ZLBs, we replace them with arbitrarily small
# lengths, and adjust the root accordingly.
replace.ZLBs <- function(tree, addtime = 0.001) {
  tree.noZLBs <- tree
  if (any(tree.noZLBs$edge.length == 0)) {
    wh.ZLBs <- which(tree.noZLBs$edge.length == 0)
    tree.noZLBs$edge.length[wh.ZLBs] <- addtime
    tree.noZLBs <-
      fixRootTime(tree, tree.noZLBs, fixingMethod = "rescaleUsingTipToRootDist")
  }
  return(tree.noZLBs)
}

# Are polytomies present?
table(!is.binary.multiPhylo(cal3trees)) # No polytomies!

# How many zero-length branches (ZLBs)?
sq <- 1:length(cal3trees)
summary(sapply(sq, function(sq) 
  length(which(cal3trees[[sq]]$edge.length == 0))))
# 43 - 78 ZLBs, mean = 61

# Confirm that replacing ZLBs with 0.001 Myr is negligible
summary(sapply(sq, function(sq) 
  min(cal3trees[[1]]$edge.length[cal3trees[[1]]$edge.length > 0])))
# smallest branch is 0.1 Myr, so 0.001 is negligible

# Effect of removing ZLBs on terminal and node branch lengths
cal3trees.noZLBs <- lapply(cal3trees, replace.ZLBs)
summary(br <- as.vector(sapply(sq, function(sq) 
  compareTermBranches(cal3trees[[sq]], cal3trees.noZLBs[[sq]]))))
summary(nd <- as.vector(sapply(sq, function(sq) 
  compareNodeAges(cal3trees[[sq]], cal3trees.noZLBs[[sq]]))))
100 * round(table(round(nd, 3)) / length(nd), 3)
hist(nd)
# 0 branch adjustments; -0.010 - + 0.017 (median = 0.001) node adjustments

# For ancestral reconstructions below, use the tree without ZLBs
cal3trees <- cal3trees.noZLBs
summary(sapply(sq, function(sq) 
  length(which(cal3trees[[sq]]$edge.length == 0)))) # Now 0 ZLBs


## IMPORT DATA SETS ############################################################

## Import morphological and ecological data sets in NEXUS format
# Note that although in NEXUS format, these files do not contain phylogenetic
# data. All phylogenies are drawn from the time-trees above.
# input <- "EchinoTree_Mode.nex"
# input <- "EchinoTree_Constant.nex"
# input <- "EchinoTree_Raw.nex"
input <- "EchinoTree_Morph.nex"
DataMatrix <- read_nexus_matrix(file_name = input)

# View to confirm imported correctly
str(DataMatrix)
DataMatrix$matrix_1$matrix[1:10, 1:10]


## ESTIMATE ANCESTORS ##########################################################

# Lloyd (2018, Palaeontology) recommends the following ("pre-OASE 1")
# implementation because it minimizes diluting morphological (ecological) signal
# in the data matrix. It only reconstructs ancestral states when tips are known,
# and does not subsequently estimate tips with unknown states. (It's also the
# computationally fastest, among slow algorithms.) The function uses
# phytools:rerootingMethod (which wraps partially around the ancestral character
# estimation method ape::ace) to calculate the marginal ancestral state estimate
# for nodes using Bayesian maximum likelihood (Yang et al., 1995). The method
# uses a continuous-time Markov chain (Mk) model of state changes based on the
# state-transition matrix Q. Using the equal-rate ("ER") model of state changes
# because we have no prior expectation of variable-rate or symmetric changes.

# The code in Claddis is very slow (technically, the code in
# phytools::rerootingMethod is the culprit), so am going to unpack the raw code
# here in order to run the ape::ace() portion 'in parallel'. The relevant line
# of code to reconstruct ancestors is line 334 at
# https://github.com/graemetlloyd/Claddis/blob/master/R/estimate_ancestral_states.R,
# or starting with line 420 here, where 'snowfall' and 'parallel' packages are
# input. Thanks, Graeme!

# Function does following:
#   * Appends tree to first element [[1]].
#   * Appends node reconstructions to the data blocks in $matrix[[2-n]].
# This means it is sensible to run the function on individual blocks of data
# (even portions of a block) in bite-size chunks, then combining the blocks back
# together for analyses of disparity.

# To make the processing slightly faster, we also concatenate all time-scaled
# tree character matrices into a single list after pro-processing, then
# re-seperate them afterwards. This avoids the inevitable load-balancing
# scenario where the parallel processing awaits only one or two slow (=difficult
# to reconstruct) characters. In other words, this allows the process to move on
# to the next tree while solving earlier ones.

## 1. PRO-PROCESSING OF ANCESTRAL DATA #########################################

# 1. Run pre-processing in a loop (to build the data list for each time tree)
num.trees <- length(cal3trees)
pre_all_data_lists <- vector("list", num.trees)

for(tr in 1:num.trees) {
  # Set args:
  time_tree = cal3trees[[tr]]
  cladistic_matrix = DataMatrix
  class(cladistic_matrix) <- "cladisticMatrix" # Forthcoming Claddis v. behavior
  estimate_all_nodes = FALSE
  estimate_tip_values = FALSE
  inapplicables_as_missing = FALSE
  polymorphism_behaviour = "equalp"
  uncertainty_behaviour = "equalp"
  threshold = 0.01
  all_missing_allowed = TRUE

  # Pre-processing troubleshooting, functions, and data preparations: Lines 80 -
  # 369 at https://github.com/graemetlloyd/Claddis/blob/master/R/estimate_ancestral_states.R
  # , committed April 7, 2021
  {
    # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
    if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
    
    # Get number of tips in tree:
    n_tips <- ape::Ntip(phy = time_tree)
    
    # Get number of nodes in tree:
    n_nodes <- ape::Nnode(phy = time_tree)
    
    # Catch problem with trees with no branch lengths:
    if (is.null(time_tree$edge.length)) stop("time_tree must have branch lengths.")
    
    # Catch problem with polytomies:
    if (time_tree$Nnode < (n_tips - 1)) stop("time_tree must be fully bifurcating.")
    
    # Catch problem with zero-length branches:
    if (any(time_tree$edge.length == 0)) stop("time_tree must not have zero-length branches.")
    
    # Check for step matrices and stop and warn if found:
    if (length(x = cladistic_matrix$topper$step_matrices) > 0) stop("Function can not currently deal with step matrices.")
    
    # Check estimate_all_nodes is a logical:
    if (!is.logical(estimate_all_nodes)) stop("estimate_all_nodes must be a logical (TRUE or FALSE).")
    
    # Check estimate_tip_values is a logical:
    if (!is.logical(estimate_tip_values)) stop("estimate_tip_values must be a logical (TRUE or FALSE).")
    
    # Check inapplicables_as_missing is a logical:
    if (!is.logical(inapplicables_as_missing)) stop("inapplicables_as_missing must be a logical (TRUE or FALSE).")
    
    # Check polymorphism_behaviour is a single allowable value:
    if (length(x = polymorphism_behaviour) != 1 || !any(c("equalp", "treatasmissing") == polymorphism_behaviour)) stop("polymorphism_behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
    
    # Check uncertainty_behaviour is a single allowable value:
    if (length(x = uncertainty_behaviour) != 1 || !any(c("equalp", "treatasmissing") == uncertainty_behaviour)) stop("uncertainty_behaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
    
    # Check threshold is a numeric value between the limits of zero and one:
    if (!is.numeric(threshold) || threshold > 0.5 || threshold < 0) stop("threshold must be a numeric value between 0 and 0.5.")
    
    # Collapse matrix to vectors for each character (state and ordering combination):
    collapsed_matrix <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], function(x) apply(rbind(x$matrix, x$ordering), 2, paste, collapse = ""))))
    
    # Isolate ordering elements:
    ordering <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering"))
    
    # Isolate minimum values:
    minimum_values <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values"))
    
    # Isolate maximum values:
    maximum_values <- unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values"))
    
    # Store raw original matrix:
    raw_cladistic_matrix <- cladistic_matrix
    
    # Combine matrix blocks into a single matrix:
    cladistic_matrix <- original_matrix <- do.call(what = cbind, args = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"))
    
    # Find any failed name matches:
    failed_name_matches <- c(setdiff(x = rownames(x = cladistic_matrix), y = time_tree$tip.label), setdiff(x = time_tree$tip.label, y = rownames(x = cladistic_matrix)))
    
    # Check there are no failed name matches and stop and report if found:
    if (length(x = failed_name_matches) > 0) stop(paste("The following names do not match between the tree and matrix: ", paste(sort(x = failed_name_matches), collapse = ", "), ". Check spelling and try again.", sep = ""))
    
    # If treating inapplicables as missing (and there is at least one inapplicable) replace with NA:
    if (inapplicables_as_missing && length(x = which(x = cladistic_matrix == "")) > 0) cladistic_matrix[which(x = cladistic_matrix == "")] <- NA
    
    # If treating polymorphisms as missing:
    if (polymorphism_behaviour == "treatasmissing" && length(x = grep("&", cladistic_matrix)) > 0) cladistic_matrix[grep("&", cladistic_matrix)] <- NA
    
    # If treating uncertainties as missing:
    if (uncertainty_behaviour == "treatasmissing" && length(x = grep("/", cladistic_matrix)) > 0) cladistic_matrix[grep("/", cladistic_matrix)] <- NA
    
    # Get vector of character numbers where all values are NA:
    dataless_characters <- which(x = apply(cladistic_matrix, 2, function(x) all(is.na(x))))
    
    # Look for all missing characters and stop and warn user if found:
    if (!all_missing_allowed && length(x = dataless_characters) > 0) stop(paste0("The following characters are coded as missing across all tips: ", paste0(dataless_characters, collapse = ", "), ". This can arise either because of the input data (in which case it is recommended that the user prune these characters using prune_cladistic_matrix) or because of the chosen options for inapplicables_as_missing, polymorphism_behaviour, and/or uncertainty_behaviour (in which case the user may wish to chose different values for these)."))
    
    # Convert tip states into a list:
    data_list <- apply(cladistic_matrix, 2, function(x) list(tip_states = x))
    
    # For each character:
    for (i in 1:length(x = data_list)) {
      
      # Add minimum value to list:
      data_list[[i]]$minimum_values <- unname(minimum_values[i])
      
      # Add maximum value to list:
      data_list[[i]]$maximum_values <- unname(maximum_values[i])
      
      # Add ordering to list:
      data_list[[i]]$ordering <- unname(ordering[i])
      
      # Add tree to list:
      data_list[[i]]$tree <- time_tree
    }
    
    # If estimating values for all characters (need to set dummy tip states for missing values):
    if (estimate_all_nodes) {
      
      # Subfunction to fill missing values (and inapplicables if desired):
      fill_missing <- function(tip_states) {
        
        # Find which rows correspond to missing states:
        missingRows <- which(x = is.na(tip_states$tip_states))
        
        # If missing states found:
        if (length(x = missingRows) > 0) {
          
          # Build missing state by either forming a polymorphism of all possible tip states, or if continuous the midpoint value:
          fill_states <- ifelse(tip_states$ordering == "continuous", (tip_states$minimum_values + tip_states$maximum_values) / 2, paste(tip_states$minimum_values:tip_states$maximum_values, collapse = "/"))
          
          # Insert missing values:
          tip_states$tip_states[missingRows] <- fill_states
        }
        
        # Return tip states with missing values replaced:
        tip_states
      }
      
      # Apply fill missing function across all characters:
      data_list <- lapply(X = data_list, fill_missing)
    }
    
    # Subfunction to prune tips with missing or inapplicable values:
    prune_tips <- function(x) {
      
      # Find all missing or inapplicable value tip names:
      missing <- names(sort(x = c(which(x = x$tip_states == ""), which(x = is.na(x$tip_states)))))
      
      # Work out how many tips will be left after pruning:
      n_tips_remaining <- length(x = setdiff(x = names(x$tip_states), y = missing))
      
      # If there is at least one missing value:
      if (length(x = missing) > 0) {
        
        # If less than two tips will remain then set tree as NULL:
        if (n_tips_remaining < 2) x$tree <- NULL
        
        # If at least two tips will remain prune missing values from tree:
        if (n_tips_remaining > 1) x$tree <- ape::drop.tip(phy = x$tree, tip = missing)
        
        # Collapse tip states:
        x$tip_states <- x$tip_states[setdiff(x = names(x$tip_states), y = missing)]
      }
      
      # Return pruned output:
      x
    }
    
    # Prune out missing and inapplicable tips:
    data_list <- lapply(X = data_list, prune_tips)
    
    # Subfunction to build tip state matrices:
    convert_tip_states_to_matrix <- function(x) {
      
      # As long asthere is at least one tip state:
      if (length(x = x$tip_states) > 0) {
        
        # If the character is not continuous (i.e., it is some form of discrete character):
        if (x$ordering != "continuous") {
          
          # Temporarily store tip states so matrix format can overwrite the stored version below:
          tip_states <- x$tip_states
          
          # Create matrix of tip state probabilities:
          x$tip_states <- matrix(0, nrow = length(x = x$tip_states), ncol = x$maximum_values - x$minimum_values + 1, dimnames = list(names(x$tip_states), x$minimum_values:x$maximum_values))
          
          # For each character state if a single state is coded store probability as 1:
          for (i in colnames(x = x$tip_states)) x$tip_states[tip_states == i, i] <- 1
          
          # If there are polymorphisms and/or uncertainties:
          if (length(x = grep("&|/", tip_states)) > 0) {
            
            # Get polymorphism locations:
            polymorphisms <- grep("&", tip_states)
            
            # Get uncertainty locations:
            uncertainties <- grep("/", tip_states)
            
            # If there are polymorphisms and using the "equalp" (equal probability of each state) option:
            if (length(x = polymorphisms) > 0 && polymorphism_behaviour == "equalp") {
              
              # For each polymorphisms set each state as equally probable:
              for (i in polymorphisms) x$tip_states[i, strsplit(tip_states[i], split = "&")[[1]]] <- 1 / length(x = strsplit(tip_states[i], split = "&")[[1]])
            }
            
            # If there are uncertainties and using the "equalp" (equal probability of each state) option:
            if (length(x = uncertainties) > 0 && uncertainty_behaviour == "equalp") {
              
              # For each uncertainty set each state as equally probable:
              for (i in uncertainties) x$tip_states[i, strsplit(tip_states[i], split = "/")[[1]]] <- 1 / length(x = strsplit(tip_states[i], split = "/")[[1]])
            }
          }
          
          # If a continuous character:
        } else {
          
          # Simply make tip states the numeric values (should never be a polymorphism) as a vector:
          x$tip_states <- as.numeric(x$tip_states)
        }
        
        # If tip state has no length (all values are missing):
      } else {
        
        # Create row-less tip states matrix:
        x$tip_states <- matrix(nrow = 0, ncol = 1, dimnames = list(c(), "0"))
      }
      
      # Return the revised input in the same list format:
      x
    }
    
    # Reformat tip states ready for ancestral estimation:
    data_list <- lapply(X = data_list, convert_tip_states_to_matrix)
    
    # Subfunction to build character model:
    build_character_model <- function(x) {
      
      # Set default model to equal rates (works for all binary or unordered characters):
      x$model <- "ER"
      
      # If a character is both ordered and has at least three states:
      if ((x$maximum_values - x$minimum_values) > 1 && x$ordering == "ord") {
        
        # Get number of states:
        n_states <- (x$maximum_values - x$minimum_values) + 1
        
        # Build all zero matrix to begin with:
        x$model <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(x$minimum_values:x$maximum_values, x$minimum_values:x$maximum_values))
        
        # for each (just) off-diagonal value store 1 (i.e., N steps to move between adjacent states):
        for (i in 2:n_states) x$model[(i - 1), i] <- x$model[i, (i - 1)] <- 1
      }
      
      # Return full output:
      x
    }
    
    # Add ancestral state model for each character:
    data_list <- lapply(X = data_list, build_character_model)
    
    # Subfunction to get ancestral states:
    estimate_ancestral_state <- function(x, estimate_tip_values, threshold) {
      
      # As long as there is a tree:
      if (!is.null(x$tree)) {
        
        # If character is continuous:
        if (x$ordering == "continuous") {
          
          # Get ancestral states using ace:
          x$ancestral_states <- ace(x = x$tip_states, phy = x$tree)$ace
          
          # If character is discrete:
        } else {
          
          # If invariant character:
          if (ncol(x$tip_states) == 1) {
            
            # Get number of tips:
            n_tips <- ape::Ntip(phy = x$tree)
            
            # Set ancestral states as all the same:
            x$ancestral_states <- matrix(rep(x = 1, times = (n_tips + x$tree$Nnode)), ncol = 1, dimnames = list(c(x$tree$tip.label, (n_tips + 1):(n_tips + x$tree$Nnode)), colnames(x = x$tip_states)))
          }
          
          # If variant character then get ancestral states using rerooting method:
          if (ncol(x$tip_states) > 1) x$ancestral_states <- phytools::rerootingMethod(tree = x$tree, x = x$tip_states, model = x$model)$marginal.anc
          
          # Reformat to most likely state
          x$ancestral_states <- unlist(x = lapply(X = lapply(X = apply(x$ancestral_states, 1, list), unlist), function(x) {
            paste(names(x[x > (max(x) - threshold)]), collapse = "/")
          }))
          
          # If not estimating tip values then prune these:
          if (!estimate_tip_values) x$ancestral_states <- x$ancestral_states[-match(x$tree$tip.label, names(x$ancestral_states))]
        }
        
        # If no tree:
      } else {
        
        # Set ancestral states as NULL:
        x$ancestral_states <- vector(mode = "character")
      }
      
      # Return full output of x:
      x
    }  
  }

  # Store data_list in list of data_lists for next step
  pre_all_data_lists[[tr]] <- data_list
}
beep(3)

# For the morphological data set (ONLY), which required processing as
# individually saved and later re-loaded bites (see below for details), need to
# save and reload several objects from above for post-processing
# save(original_matrix, file = "original_matrix")

# Each list in 'all_data_lists' is a separate character matrix and tree. Here we
# combine them into a single 'data_list' for processing the next step 'in
# parallel.'
nchar <- length(pre_all_data_lists[[1]])
data_list <- unlist(pre_all_data_lists, recursive = FALSE)

# Confirm concatenated correctly
if(identical(data_list[[1]], pre_all_data_lists[[1]][[1]]) |
   identical(data_list[[length(data_list)]],
             pre_all_data_lists[[length(pre_all_data_lists)]][[length(pre_all_data_lists[[length(pre_all_data_lists)]])]]))
  cat( "The 'data_list's were concatenated correctly.\n There are",
    length(data_list), "total characters to reconstruct." ) else
  stop("The 'data_list's were NOT concatenated correctly.")


## BREAK MORPHOLOGICAL DATA INTO BITE-SIZED PIECES FOR PARALLEL PROCESSING #####

# Here there be monsters! For processing the (large character-number)
# morphological data set, we are breaking into five bite-sized pieces to
# ameliorate instances of processing failure. (In other words, if something goes
# awry, like a power failure, we've only lost a day of run time instead of
# having to re-do everything.) Because the objects are large, also removing all
# unnecessary prior objects from memory, which requires re-scripting the
# required sub-functions and args from above code.
if(length(data_list) == 20650L) {
  all_data_lists <- data_list
  # Save for safekeeping (note a large gigabyte object so will take time to
  # save/load)
  # save(all_data_lists, file = "all_data_lists")
  # load("all_data_lists") # Used when restarting for trials 2-5
  beep()
  bite_size <- length(all_data_lists) / 5
  cat("Breaking down into", bite_size, 
      "characters per each of five runs. Make sure to redefine \n below each 'data_list' accordingly when running in parallel.\n")
  bites <- seq(from = 1, to = length(all_data_lists), by = bite_size)
  bite_to <- bites + bite_size - 1
  # To save further memory, it is advisable to manually only define one bite at
  # a time (instead of pre-defining all 5 at once)
  # data_list1 <- all_data_lists[bites[1]:bite_to[1]]
  # data_list2 <- all_data_lists[bites[2]:bite_to[2]]
  # data_list3 <- all_data_lists[bites[3]:bite_to[3]]
  # data_list4 <- all_data_lists[bites[4]:bite_to[4]] # * ERRORS! *
  # data_list4a <- all_data_lists[12391:14455] # But this smaller half-bite works !
  data_list4b <- all_data_lists[14456:16520]
  # data_list4b <- all_data_lists[14456:15456]
  # data_list4c <- all_data_lists[15457:16000]
  # data_list4d <- all_data_lists[16001:16520]
  # data_list5 <- all_data_lists[bites[5]:bite_to[5]]
  # Clean up any no-longer used large objects to save working memory
  rm(list = c("pre_all_data_lists", "data_list", "all_data_lists", "cal3trees", 
              "cal3trees.noZLBs", "DataMatrix", "br", "nd"))
  gc()
  # Re-assign critical args and subfunction used in parallel portion:
  estimate_tip_values = FALSE
  threshold = 0.01
  estimate_ancestral_state <- function(x, estimate_tip_values, threshold) {
    if (!is.null(x$tree)) {
      if (x$ordering == "continuous") {
        x$ancestral_states <- ace(x = x$tip_states, phy = x$tree)$ace
      } else {
        if (ncol(x$tip_states) == 1) {
          n_tips <- ape::Ntip(phy = x$tree)
          x$ancestral_states <- matrix(rep(x = 1, times = (n_tips + x$tree$Nnode)), ncol = 1, dimnames = list(c(x$tree$tip.label, (n_tips + 1):(n_tips + x$tree$Nnode)), colnames(x = x$tip_states)))
        }
        if (ncol(x$tip_states) > 1) x$ancestral_states <- phytools::rerootingMethod(tree = x$tree, x = x$tip_states, model = x$model)$marginal.anc
        x$ancestral_states <- unlist(x = lapply(X = lapply(X = apply(x$ancestral_states, 1, list), unlist), function(x) {
          paste(names(x[x > (max(x) - threshold)]), collapse = "/")
        }))
        if (!estimate_tip_values) x$ancestral_states <- x$ancestral_states[-match(x$tree$tip.label, names(x$ancestral_states))]
      }
    } else {
      x$ancestral_states <- vector(mode = "character")
    }
    x
  }  
}

# Manually update accordingly (only for the morphological data set). If need to
# restart, load the 'all_data_lists' object above and rebuild 'data_list1',
# 'data_list2', etc.
# data_list <- data_list1; rm("data_list1"); gc()
# data_list <- data_list2; rm("data_list2"); gc()
# data_list <- data_list3; rm("data_list3"); gc()
# data_list <- data_list4; rm("data_list4"); gc()
# data_list <- data_list4a; rm("data_list4a"); gc()
# data_list <- data_list4b; rm("data_list4b"); gc()
# data_list <- data_list4c; rm("data_list4c"); gc()
# data_list <- data_list4d; rm("data_list4d"); gc()
# data_list <- data_list5; rm("data_list5"); gc()

## 2. PARALLEL PROCESSING OF ANCESTRAL STATE RECONSTRUCTION ####################

# Parallel implementation to get ancestral states for each character
library(snowfall)
library(parallel)
(t.start0 <- Sys.time())
# Set up computer cluster using all available cores
cpus <- min(parallel::detectCores(), length(data_list))
# Initialize cluster
sfInit(parallel = TRUE, cpus = cpus, slaveOutfile = "initfile")
stopifnot(sfCpus() == cpus)		    # Confirm set up CPUs properly
stopifnot(sfParallel() == TRUE)		# Confirm now running in parallel
sfExportAll()				              # Export all libraries, files, & objects
sfLibrary(phytools)
par.out <- NA
# Version without load-balancing
par.out <- sfLapply(x = data_list, fun = estimate_ancestral_state, 
                    estimate_tip_values = estimate_tip_values, 
                    threshold = threshold)
sfStop()
# Save for safekeeping
save(par.out, file = "par.out")

(Sys.time() - t.start0)
beep(3)
# Timing log:
# 1.88 days for Ecology_Mode and no errors
# 1.60 days for Ecology_Constant and no errors
# 12.8 hrs for Ecology_Raw and no errors (characters 7-8 were all missing and added manually below)
# 23.12 hrs + 22.28 hrs  + 22.36 hrs  + (13.21 + 12.68 = 25.89) hrs + 22.69 hrs =  > 4.8 days for Morph (8 characters were all missing and added manually below)
#  - processed in 5 bites; no error for bites 1-3 & 5, but error in bite 4 which
#      was divided into two half-bites
#  - bite 4 produced an error ("one node produced an error: incorrect number of
#      dimensions" Init file shows 4 instances of "log(comp[1:M + N]) : NaNs
#      produced")
# No errors when run 12391:14455 or 14456:15456 or 15457:16000 or 16001:16520,
#    so problem character must be caused by some combination of stochasticity
#    and/or peculiarities with that bite of characters.
warnings()
str(par.out[[500]])   

# Check for any all-missing states in output above (check any that are
# returned). Following adds a tree and all-missing states to relevant
# characters. Make sure the tree is the one without ZLBs, and that it
# sequentially updates based on the time tree used.

# First check for all-missing states
sq <- seq.int(par.out)
wh.all.missing <-
  which(!unlist(lapply(sq, function(sq) ! is.null(par.out[[sq]]$tree))))
if (length(wh.all.missing) == 0L)
  cat("No characters were skipped when inferring ancestral states.\n")
if (length(wh.all.missing) > 0L)
  cat("Character(s)", wh.all.missing, "were skipped when inferring ancestral states. 
      Trees and missing ancestral states were manually added to these characters.\n")

# Characters 7 & 8 in all trees for 'raw'
# Characters 34, 35, 39, 43, 51, 52, 54, & 410 in all trees for 'morphology'
# If all NAs, do NOT do the manual tree-add below

# And process, if so ...
if (length(wh.all.missing) > 0L) {
  # *** Start by assigning the tree that was first processed. (Will be 1 for the
  # 'raw' treatment and start with seq(from = 1, to = 50, by = 10) for the
  # 5-bite 'morphology' data set.) Also make sure the number of characters for
  # the data set is correctly specified.
  starting.tree <- 36     # *** ! CRITICAL - DO NOT GET THIS WRONG! ***
  nchar <- 413            # *** ! CRITICAL - DO NOT GET THIS WRONG! ***
  tree.index <- floor(wh.all.missing / nchar) + starting.tree
  # Confirm the available list of time-trees lacks ZLBs
  if (!exists("cal3trees"))
    load("~/Manuscripts/CamOrdEchinos/cal3trees")
  sq <- 1:length(cal3trees)
  if (any(sapply(sq, function(sq) length(which(cal3trees[[sq]]$edge.length == 0)))) > 0)
    cal3trees <- lapply(cal3trees, replace.ZLBs)
  # Finally add the correct tree and NAs for all missing characters
  for (ch in 1:length(wh.all.missing)) {
    time_tree <- cal3trees[[tree.index[ch]]]
    par.out[[wh.all.missing[ch]]]$tree <- time_tree
    par.out[[wh.all.missing[ch]]]$ancestral_states <-
      as.character(rep(x = NA, times = Nnode(time_tree)))
    names(par.out[[wh.all.missing[ch]]]$ancestral_states) <-
      as.character((Ntip(time_tree) + 1):(Ntip(time_tree) + Nnode(time_tree)))
  }
}

# Check some all-missing characters (compared to some not)
par.out[[33]]
par.out[[34]]

# ONLY USED FOR MORPHOLOGICAL DATA SETS:
# 1. Redefine and save par.out
# par.out1 <- par.out; save(par.out1, file = "par.out1"); beep()
# par.out2 <- par.out; save(par.out2, file = "par.out2"); beep()
# par.out3 <- par.out; save(par.out3, file = "par.out3"); beep()
# par.out4 <- par.out; save(par.out4, file = "par.out4"); beep()
# par.out4a <- par.out; save(par.out4a, file = "par.out4a"); beep()
# par.out4b <- par.out; save(par.out4b, file = "par.out4b"); beep()
# par.out5 <- par.out; save(par.out5, file = "par.out5"); beep()
# 2. Now reload, combine, and verify worked as intended:
load("par.out1"); load("par.out2"); load("par.out3"); load("par.out4a"); load("par.out4b"); load("par.out5")
beep()
par.out <- c(par.out1, par.out2, par.out3, par.out4a, par.out4b, par.out5)
if(length(par.out) != 413 * 50) stop("'par.out' was not combined correctly!/n")
object.size(par.out) # HUGE! 7.4 gigabytes!
rm(list = c("par.out1", "par.out2", "par.out3", "par.out4a", "par.out4b", "par.out5")); gc() # Clean memory
# 3. Reload and rebuild necessary objects from pre-processing
load("original_matrix"); num.trees <- 50; nchar <- 413
if (!exists("cal3trees")) load("~/Manuscripts/CamOrdEchinos/cal3trees")
sq <- 1:length(cal3trees)
if (any(sapply(sq, function(sq) length(which(cal3trees[[sq]]$edge.length == 0)))) > 0)
  cal3trees <- lapply(cal3trees, replace.ZLBs)

# Reassemble into (tree) list of (character) lists
postpar_data_list <- vector("list", num.trees)
seq.start <- seq(from = 1, to = num.trees * nchar, by = nchar)
seq.end <- seq(from = nchar, to = num.trees * nchar, by = nchar)
for(tr in 1:num.trees) {
  postpar_data_list[[tr]] <- par.out[seq.start[tr]:seq.end[tr]]
}

# Confirm worked as intended
if (length(postpar_data_list) != num.trees)
  stop("The list of lists did NOT rebuild correctly.")
sq <- seq.int(num.trees)
if (any(sapply(sq, function(sq) length(postpar_data_list[[sq]]) != nchar)))
  stop("The list of lists did NOT rebuild correctly.")

# Always efficient to save memory
rm("par.out")
gc()

## 3. POST-PROCESSING OF ANCESTRAL STATE MATRICES ##############################

# Post-processing (lines 373-485 at
# https://github.com/graemetlloyd/Claddis/blob/master/R/estimate_ancestral_states.R)
ancestral_state_matrices <- postpar_data_list

# Reload (if did not run all code above as a continuous workflow)
# input <- "EchinoTree_Mode.nex"
# input <- "EchinoTree_Constant.nex"
# input <- "EchinoTree_Raw.nex"
input <- "EchinoTree_Morph.nex"
raw_cladistic_matrix <- read_nexus_matrix(file_name = input)

for(tr in 1:num.trees) {
  cat("post-processing tree", tr, "of", num.trees, "\n")
  # Redefine necessary objects based on current tree:
  data_list = postpar_data_list[[tr]]
  time_tree = cal3trees[[tr]]
  n_tips <-  ape::Ntip(phy = time_tree)
  n_nodes <- ape::Nnode(phy = time_tree)
  {
    # Get Newick strings of all sampled subtrees (to use to avoid redundancy in tree node mapping):
    newick_strings <- unlist(x = lapply(X = data_list, function(x) ifelse(is.null(x$tree), NA, ape::write.tree(x$tree))))
    
    # Get just unique strings (i.e., the minimum set needded to map everything to the full tree):
    unique_newick_strings <- unique(x = newick_strings[!is.na(newick_strings)])
    
    # Convert unique Newick strings to unique trees:
    unique_trees <- ape::read.tree(text = unique_newick_strings)
    
    # If only a single tree reformat as a list:
    if (inherits(unique_trees, what = "phylo")) unique_trees <- list(unique_trees)
    
    # Subfunction to map nodes from pruned tree to full tree:
    map_to_full_tree <- function(pruned_tree, full_tree) {
      
      # Get number of tips of pruned tree:
      n_tips <- ape::Ntip(phy = pruned_tree)
      
      # Get number of nodes of pruned tree:
      n_nodes <- ape::Nnode(phy = pruned_tree)
      
      # Get all internal node numbers for pruned tree:
      node_numbers <- (n_tips + 1):(n_tips + n_nodes)
      
      # If the pruned tree is different to the full tree:
      if (write.tree(pruned_tree) != write.tree(full_tree)) {
        
        # Get descendants of each node in pruned tree:
        node_descendants <- lapply(X = as.list(x = node_numbers), function(x) pruned_tree$tip.label[strap::FindDescendants(n = x, tree = pruned_tree)])
        
        # Get corresponding ancestral node in full tree:
        ancestral_nodes <- unlist(x = lapply(X = node_descendants, function(x) find_mrca(descendant_names = x, tree = full_tree)))
        
        # If pruned tree is identical to full tree (not pruned at all):
      } else {
        
        # Set ancestors as node numbers:
        ancestral_nodes <- node_numbers
      }
      
      # Output matrix matching node numbers of pruned tree to full tree:
      matrix(c(node_numbers, ancestral_nodes), ncol = 2, dimnames = list(c(), c("pruned_node", "full_node")))
    }
    
    # Get pruned node to full node for each unique tree:
    node_maps <- lapply(X = unique_trees, map_to_full_tree, full_tree = time_tree)
    
    # Build out for all trees (adds in any duplicated trees):
    node_maps <- node_maps[match(newick_strings, unique_newick_strings)]
    
    # Add node maps to data list:
    for (i in 1:length(x = data_list)) data_list[[i]]$node_maps <- node_maps[[i]]
    
    # Get all node names and numbers:
    nodes <- c(rownames(x = original_matrix), (n_tips + 1):(n_tips + n_nodes))
    
    # Renumber nodes of ancestral states:
    data_list <- lapply(X = data_list, function(x) {
      
      # Renumber nodes:
      names(x$ancestral_states)[match(as.character(x$node_maps[, "pruned_node"]), names(x$ancestral_states))] <- as.character(x$node_maps[, "full_node"])
      
      # Return renumbered nodes:
      x
    })
    
    # Collapse down to an ancestral state matrix ready for output:
    ancestral_state_matrix <- do.call(what = cbind, args = lapply(X = data_list, function(x) {
      
      # Get ancestral states for nodes:
      x$ancestral_states <- x$ancestral_states[nodes]
      
      # Add node names:
      names(x$ancestral_states) <- nodes
      
      # Return output:
      x$ancestral_states
    }))
    
    # Isolate estimated tip values:
    tip_matrix <- ancestral_state_matrix[rownames(x = original_matrix), ]
    
    # If there are any missing values:
    if (any(is.na(tip_matrix))) {
      
      # Isolate missing values:
      missing_tip_states <- which(x = is.na(tip_matrix))
      
      # Replace missing values with original (unmodified) input values:
      tip_matrix[missing_tip_states] <- original_matrix[missing_tip_states]
      
      # Add tip values back into full output:
      ancestral_state_matrix[rownames(x = original_matrix), ] <- tip_matrix
    }
    
    # Get column (character) count for each matrix block:
    matrix_columns <- unlist(x = lapply(X = lapply(X = raw_cladistic_matrix[2:length(x = raw_cladistic_matrix)], "[[", "matrix"), ncol))
    
    # For each matrix block:
    for (i in 1:length(x = matrix_columns)) {
      
      # Insert portion of ancestral state estimate into block:
      raw_cladistic_matrix[[(i + 1)]]$matrix <- ancestral_state_matrix[, 1:matrix_columns[i], drop = FALSE]
      
      # Remove that portion from the block:
      ancestral_state_matrix <- ancestral_state_matrix[, -(1:matrix_columns[i]), drop = FALSE]
    }
    
    # Overwrite ancestral state output with updated raw input:
    ancestral_state_matrix <- raw_cladistic_matrix
    
    # Add tree to output:
    ancestral_state_matrix$topper$tree <- time_tree
    
  }
  ancestral_state_matrices[[tr]] <- ancestral_state_matrix
}

# Save processed data
# mode.anc <- ancestral_state_matrices; save(mode.anc, file = "mode.anc"); load("mode.anc")
# constant.anc <- ancestral_state_matrices; save(constant.anc, file = "constant.anc"); load("constant.anc")
# raw.anc <- ancestral_state_matrices; save(raw.anc, file = "raw.anc"); load("raw.anc")
morph.anc <- ancestral_state_matrices; save(morph.anc, file = "morph.anc"); load("morph.anc")

beep(3)

# Confirm worked as intended
if (length(ancestral_state_matrices) != num.trees)
  stop("The list of lists did NOT rebuild and post-process correctly.")
sq <- seq.int(num.trees)
if (any(sapply(sq, function(sq) length(postpar_data_list[[sq]]) != nchar)))
  stop("The list of lists did NOT rebuild and post-process correctly.")

# Confirm everything looks OK (using tree 50)
str(ancestral_state_matrices[[50]])
names(ancestral_state_matrices[[50]]) # Should be 'topper' and 'matrix_1'

# Check states for tips and nodes (using tree 50)
ancestral_state_matrices[[50]]$matrix_1$matrix[c(1:5, 726:731), 1:10]

# Convert ancestral state matrix to csv for viewing outside R (using ONLY tree 50)
# write.csv(mode.anc[[50]]$matrix_1$matrix, file = "modeanc.csv")
# write.csv(constant.anc[[50]]$matrix_1$matrix, file = "constantanc.csv")
# write.csv(raw.anc[[50]]$matrix_1$matrix, file = "rawanc.csv")
# write.csv(morph.anc[[50]]$matrix_1$matrix, file = "morphanc.csv")
