## USING CLADDIS TO INFER ANCESTRAL STATES #####################################

## Prior to running, run 1-MakeTimeTrees.R to build time-calibrated trees using
## 'paleotree' package.

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 7/2/2020, v.
# 0.4.1. Future users will need to either download the archived version of the
# package from GitHub or alter the code accordingly to use the current argument
# and function names.

## PREPARATIONS ################################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")
setwd("~/Manuscripts/CamOrdEchino/")

## Download newest version of 'Claddis' from GitHub
# To implement the Hopkins and St. John (2016) analyses, make sure to update to
# newest verion (0.4.0) only available on GitHub.
library(devtools)
devtools::install_github("graemetlloyd/Claddis") # Downloaded 7/2/2020

## Download 'Rcpp' (dependent for phytools, but requiring binary version)
install.packages("Rcpp", dependencies = TRUE, INSTALL_opts = '--no-lock')

## Download newest version of 'phytools' from GitHub
# Implements a change to rerootingMethod() to allow two-taxon case
devtools::install_github("liamrevell/phytools") # Downloaded 7/6/2020

# Load packages
library(beepr)      # v. 1.3
library(parallel)   # v. 4.0.2
library(snowfall)   # v. 1.84-6.1
library(Rcpp)       # v. 1.0.5
library(phytools)    # v. 0.7-54
if(packageVersion("phytools") < "0.7-54")
  stop("wrong version of 'phytools' Get updated version from GitHub\n")
library(paleotree)  # v. 3.3.25
library(Claddis)    # v. 0.4.1 - Check SI for Lloyd 2018 for walk-through on code functions (in R experiments folder)
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub\n")


## Import time trees saved from prior code
load("equal.tree")
load("cal3.tree")
paleotree::phyloDiv(equal.tree); title(main = "equal tree")
paleotree::phyloDiv(cal3.tree); title(main = "cal3 tree")


# Confirm no zero-length branches present (because Claddis' morphological
# disparity functions disallow it).
any(equal.tree$edge.length == 0)     # FALSE
any(cal3.tree$edge.length == 0)      # TRUE
# CONCLUSION: We will not use the cal3 tree in analyses below. It is used solely
# to demonstrate that the phylogenetically-based diversity curve is not
# significantly different that the simpler-assumption-based "equal" tree used
# here. Because the stratigraphic ranges are essentially the same, the
# corresponding disparity trends ought to be similar. (If wish to rig the
# cal3.tree, can modify paleotree::AddTermBranchLength to add small amount to
# ZLBs. But need to focus on nodes instead of tips, and note there is an error
# in the code that would also need to be fixed.)


## Import morphological and ecological data sets in NEXUS format
# Note that although in NEXUS format, these files do not contain phylogenetic
# data. All phylogenetics are drawn from the time-trees above.
setwd("~/Manuscripts/CamOrdEchino/Data files/NA reformatted/")
input <- "EchinoTree_Mode.nex"
# input <- "EchinoTree_Constant.nex"
# input <- "EchinoTree_Raw.nex"
# input <- "EchinoTree_Morph.nex"
DataMatrix <- ReadMorphNexus(File = input)
str(DataMatrix)

DataMatrix$Matrix_1$Matrix[1:10, 1:10]





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
# phytools::rerootingMethod), so am going to unpack the raw code here in order
# to run in parallel. The relevant line of code to reconstruct ancestors is line
# 334 at github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R.
# Thanks, Graeme!

# Set args:
CladisticMatrix = DataMatrix
Tree = equal.tree
EstimateAllNodes = FALSE
EstimateTipValues = FALSE
InapplicablesAsMissing = FALSE
PolymorphismBehaviour = "equalp"
UncertaintyBehaviour = "equalp"
Threshold = 0.01

# Pre-processing troubleshooting, functions, and data preparations:
# ** Note modifications to lines 331-334 to handle characters with all-missing states
{
  
  # Catch problem with trees with no branch lengths:
  if(is.null(Tree$edge.length)) stop("Tree must have branch lengths.")
  
  # Catch problem with polytomies:
  if(Tree$Nnode < (ape::Ntip(Tree) - 1)) stop("Tree must be fully bifurcating.")
  
  # Catch problem with zero-length branches:
  if(any(Tree$edge.length == 0)) stop("Tree must not have zero-length branches.")
  
  # Check for step matrices and stop and warn if found:
  if(length(CladisticMatrix$Topper$StepMatrices) > 0) stop("Function can not currently deal with step matrices.")
  
  # Check EstimateAllNodes is a logical:
  if(!is.logical(EstimateAllNodes)) stop("EstimateAllNodes must be a logical (TRUE or FALSE).")
  
  # Check EstimateTipValues is a logical:
  if(!is.logical(EstimateTipValues)) stop("EstimateTipValues must be a logical (TRUE or FALSE).")
  
  # Check InapplicablesAsMissing is a logical:
  if(!is.logical(InapplicablesAsMissing)) stop("InapplicablesAsMissing must be a logical (TRUE or FALSE).")
  
  # Check PolymorphismBehaviour is a single allowable value:
  if(length(PolymorphismBehaviour) != 1 || !any(c("equalp", "treatasmissing") == PolymorphismBehaviour)) stop("PolymorphismBehaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
  
  # Check UncertaintyBehaviour is a single allowable value:
  if(length(UncertaintyBehaviour) != 1 || !any(c("equalp", "treatasmissing") == UncertaintyBehaviour)) stop("UncertaintyBehaviour must be a single value of either, \"equalp\" or \"treatasmissing\".")
  
  # Check Threshold is a numeric value between the limits of zero and one:
  if(!is.numeric(Threshold) || Threshold > 0.5 || Threshold < 0) stop("Threshold must be a numeric value between 0 and 1.")
  
  # Collapse matrix to vectors for each character (state and ordering combination):
  collapse.matrix <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], function(x) apply(rbind(x$Matrix, x$Ordering), 2, paste, collapse = ""))))
  
  # Isolate ordering elements:
  ordering <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Ordering"))
  
  # Isolate minimum values:
  min.vals <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MinVals"))
  
  # Isolate maximum values:
  max.vals <- unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MaxVals"))
  
  # Store raw original matrix:
  RawCladisticMatrix <- CladisticMatrix
  
  # Combine matrix blocks into a single matrix:
  CladisticMatrix <- OriginalMatrix <- do.call(cbind, lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"))
  
  # Find any failed name matches:
  FailedNameMatches <- c(setdiff(rownames(CladisticMatrix), Tree$tip.label), setdiff(Tree$tip.label, rownames(CladisticMatrix)))
  
  # Check there are no failed name matches and stop and report if found:
  if(length(FailedNameMatches) > 0) stop(paste("The following names do not match between the tree and matrix: ", paste(sort(FailedNameMatches), collapse = ", "), ". Check spelling and try again.", sep = ""))
  
  # If treating inapplicables as missing (and there is at least one inapplicable) replace with NA:
  if(InapplicablesAsMissing && length(which(CladisticMatrix == "")) > 0) CladisticMatrix[which(CladisticMatrix == "")] <- NA
  
  # If treating polymorphisms as missing:
  if(PolymorphismBehaviour == "treatasmissing" && length(grep("&", CladisticMatrix)) > 0) CladisticMatrix[grep("&", CladisticMatrix)] <- NA
  
  # If treating uncertainties as missing:
  if(UncertaintyBehaviour == "treatasmissing" && length(grep("/", CladisticMatrix)) > 0) CladisticMatrix[grep("/", CladisticMatrix)] <- NA
  
  # Convert tip states into a list:
  DataAsList <- apply(CladisticMatrix, 2, list)
  
  # Add Tip states name to list:
  DataAsList <- lapply(DataAsList, function(x) {names(x) <- "TipStates"; return(x)})
  
  # For each character:
  for(i in 1:length(DataAsList)) {
    
    # Add minimum value to list:
    DataAsList[[i]]$MinVal <- unname(min.vals[i])
    
    # Add maximum value to list:
    DataAsList[[i]]$MaxVal <- unname(max.vals[i])
    
    # Add ordering to list:
    DataAsList[[i]]$Ordering <- unname(ordering[i])
    
    # Add tree to list:
    DataAsList[[i]]$Tree <- Tree
    
  }
  
  # If estimating values for all characters (need to set dummy tip states for missing values):
  if(EstimateAllNodes) {
    
    # Subfunction to fill missing values (and inapplicables if desired):
    FillMissing <- function(TipStates) {
      
      # Find which rows correspond to missing states:
      MissingRows <- which(is.na(TipStates$TipStates))
      
      # If missing states found:
      if(length(MissingRows) > 0) {
        
        # Build missing state by either forming a polymorphism of all possible tip states, or if continuous the midpoint value:
        FillStates <- ifelse(TipStates$Ordering == "cont", (TipStates$MinVal + TipStates$MaxVal) / 2, paste(TipStates$MinVal:TipStates$MaxVal, collapse = "/"))
        
        # Insert missing values:
        TipStates$TipStates[MissingRows] <- FillStates
        
      }
      
      # Return tip states with missing values replaced:
      return(TipStates)
      
    }
    
    # Apply fill missing fucntion across all characters:
    DataAsList <- lapply(DataAsList, FillMissing)
    
  }
  
  # Subfunction to prune tips with missing or inapplicable values:
  PruneTips <- function(x) {
    
    # Find all missing or inapplicable value tip names:
    Missing <- names(sort(c(which(x$TipStates == ""), which(is.na(x$TipStates)))))
    
    # If there is at least one:
    if(length(Missing) > 0) {
      
      # Remove tips from tree:
      x$Tree <- drop.tip(phy = x$Tree, tip = Missing)
      
      # Remove tips from tip states:
      x$TipStates <- x$TipStates[setdiff(names(x$TipStates), Missing)]
      
    }
    
    # Return pruend output:
    return(x)
    
  }
  
  # Prune out missing and inapplicable tips:
  DataAsList <- lapply(DataAsList, PruneTips)
  
  # Subfunction to build tip state matrices:
  TipStateVectorToMatrix <- function(x) {
    
    # If the character is not continuous (i.e., it is some form of discrete character):
    if(x$Ordering != "cont") {
      
      # Temporarily store tip states so matrix format can overwrite the stored version below:
      TipStates <- x$TipStates
      
      # Create marix of tip state probabilities:
      x$TipStates <- matrix(0, nrow = length(x$TipStates), ncol = x$MaxVal - x$MinVal + 1, dimnames = list(names(x$TipStates), x$MinVal:x$MaxVal))
      
      # For each character state if a single state is coded store probability as 1:
      for(i in colnames(x$TipStates)) x$TipStates[TipStates == i, i] <- 1
      
      # If there are polymorphisms and/or uncertainties:
      if(length(grep("&|/", TipStates)) > 0) {
        
        # Get polymorphism locations:
        Polymorphisms <- grep("&", TipStates)
        
        # Get uncertainty locations:
        Uncertainties <- grep("/", TipStates)
        
        # If there are polymorphisms and using the "equalp" (equal probability of each state) option:
        if(length(Polymorphisms) > 0 && PolymorphismBehaviour == "equalp") {
          
          # For each polymorphisms set each state as equally probable:
          for(i in Polymorphisms) x$TipStates[i, strsplit(TipStates[i], split = "&")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "&")[[1]])
          
        }
        
        # If there are uncertainties and using the "equalp" (equal probability of each state) option:
        if(length(Uncertainties) > 0 && UncertaintyBehaviour == "equalp") {
          
          # For each uncertainty set each state as equally probable:
          for(i in Uncertainties) x$TipStates[i, strsplit(TipStates[i], split = "/")[[1]]] <- 1 / length(strsplit(TipStates[i], split = "/")[[1]])
          
        }
        
      }
      
      # If a continuous character:
    } else {
      
      # Simply make tip states the numeric values (should never be a polymorphism) as a vector:
      x$TipStates <- as.numeric(x$TipStates)
      
    }
    
    # Return the revised input in the same list format:
    return(x)
    
  }
  
  # Reformat tip states ready for ancestral estimation:
  DataAsList <- lapply(DataAsList, TipStateVectorToMatrix)
  
  # Subfunction to build tip state matrices:
  ModelBuilder <- function(x) {
    
    # Set default model to equal rates (works for all binary or unordered characters):
    x$Model <- "ER"
    
    # If a character is both ordered and has at least three states:
    if((x$MaxVal - x$MinVal) > 1 && x$Ordering == "ord") {
      
      # Get number of states:
      NStates <- (x$MaxVal - x$MinVal) + 1
      
      # Build all zero matrix to begin with:
      x$Model <- matrix(0, nrow = NStates, ncol = NStates, dimnames = list(x$MinVal:x$MaxVal, x$MinVal:x$MaxVal))
      
      # for each (just) off-diagonal value store 1 (i.e., N steps to move between adjacent states):
      for(i in 2:NStates) x$Model[(i - 1), i] <- x$Model[i, (i - 1)] <- 1
      
    }
    
    # Return full output:
    return(x)
    
  }
  
  # Add ancestral state model for each character:
  DataAsList <- lapply(DataAsList, ModelBuilder)
  
  # Sunfunction to get ancestral states:
  GetAncStates <- function(x, EstimateTipValues, Threshold) {
    
    ## *** MODIFIED FROM FUNCTION ON GITHUB TO FLAG CASES WHERE ALL STATES MISSING ***
    if (is.null(x$Tree)) {
      warning("GetAncStates could not infer ancestral states for all-missing character(s). Check manually.\n")
    } else {
      
      # If character is continuous:
      if(x$Ordering == "cont") {
        
        # Get ancestral states using ace:
        x$AncestralStates <- ace(x = x$TipStates, phy = x$Tree)$ace
        
        # If character is discrete:
      } else {
        
        # If invariant character:
        if(ncol(x$TipStates) == 1) {
          
          # Get number of tips:
          NTreeTips <- ape::Ntip(x$Tree)
          
          # Set ancestral states as all the same:
          x$AncestralStates <- matrix(rep(x = 1, times = (NTreeTips + x$Tree$Nnode)), ncol = 1, dimnames = list(c(x$Tree$tip.label, (NTreeTips + 1):(NTreeTips + x$Tree$Nnode)), colnames(x$TipStates)))
          
        }
        
        # If variant character then get ancestral states using rerooting method:
        if(ncol(x$TipStates) > 1) x$AncestralStates <- rerootingMethod(tree = x$Tree, x = x$TipStates, model = x$Model)$marginal.anc
        
        # Reformat to most likely state
        x$AncestralStates <- unlist(lapply(lapply(apply(x$AncestralStates, 1, list), unlist), function(x) {paste(names(x[x > (max(x) - Threshold)]), collapse = "/")}))
        
        # If not estimating tip values then prune these:
        if (!EstimateTipValues) x$AncestralStates <- x$AncestralStates[-match(x$Tree$tip.label, names(x$AncestralStates))]
        
      }
      
      # Return full output of x:
      return(x)
      
    }
    
  }

}

# Lines 67 - 331 in github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R

# If get following error: "In drop.tip(phy = x$Tree, tip = Missing) ...", it
# means a character is missing states for all taxa. Modified 'GetAncStates' code
# will still process, and code below will post a warning noting which characters
# were affected.

# Parallel implementation to get ancestral states for each character
library(snowfall)
library(parallel)
(t.start0 <- Sys.time())
# Set up computer cluster using all available cores
cpus <- min(parallel::detectCores(), length(DataAsList))
# Initialize cluster
sfInit(parallel = TRUE, cpus = cpus, slaveOutfile = "initfile")
stopifnot(sfCpus() == cpus)		    # Confirm set up CPUs properly
stopifnot(sfParallel() == TRUE)		# Confirm now running in parallel
sfExportAll()				              # Export all libraries, files, & objects
sfLibrary(phytools)
par.out <- NA
# Version without load-balancing
par.out <- sfLapply(x = DataAsList, fun = GetAncStates, 
                    EstimateTipValues = EstimateTipValues, Threshold = Threshold)
sfStop()
beep(3)
(Sys.time() - t.start0)
# Timing log:
# 1.05 hrs for Ecology_Mode and no errors
# 1.04 hrs for Ecology_Constant and no errors
# 0.31 hrs for Ecology_Raw [using modified algorithm above to handle characters with all-missing states]
# 3.28 hrs for Morph and no errors (with phytools v 0.7-43)
warnings()
str(par.out[[1]])
DataAsList <- par.out

# Check for any all-missing states in output above (check any that are
# returned). Check initfile to confirm.
if(length(last.warning) > 0) {
  sq <- seq.int(DataAsList)
  wh.w <- which(!unlist(lapply(sq, function(sq) !is.null(DataAsList[[sq]]$Tree))))
  cat("Check initfile for warnings. Character(s)", wh.w, "were skipped when inferring ancestral states.\n")
}

# ONLY run following if a warning was triggered above. It adds a tree and
# all-missing states to relevant characters.
# wh.missing <- c(7, 8)
# for (m in 1:length(wh.missing)) {
#   DataAsList[[wh.missing[m]]]$Tree <- Tree
#   DataAsList[[wh.missing[m]]]$AncestralStates <-
#     as.character(rep(x = NA, times = Nnode(Tree)))
#   names(DataAsList[[wh.missing[m]]]$AncestralStates) <-
#     as.character((Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree)))
# }


# Post-processing (lines 336-441 at github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R)
{
  # Get Newick strings of all sampled subtrees (to use to avoid reudunancy in tree node mapping):
  NewickStrings <- unlist(lapply(lapply(DataAsList, '[[', "Tree"), ape::write.tree))
  
  # Get just unique strings (i.e., just those trees that need toa ctyal map to full tree):
  UniqueNewickStrings <- unique(NewickStrings)
  
  # Convert unique Newick strings to unique trees:
  UniqueTrees <- read.tree(text = UniqueNewickStrings)
  
  # If only a single tree reformat as a list:
  if(inherits(UniqueTrees, what = "phylo")) UniqueTrees <- list(UniqueTrees)
  
  # Subfunction map nodes from pruned tree to full tree:
  MapPrunedTreeNodesToFullTreeNodes <- function(tree, fulltree) {
    
    # Get number of tips of pruned tree:
    NTips <- ape::Ntip(tree)
    
    # Get number of nodes of pruend tree:
    NNodes <- ape::Nnode(tree)
    
    # Get all internal node numbers for peruend tree:
    NodeNumbers <- (NTips + 1):(NTips + NNodes)
    
    # If the pruned tree is different to the full tree:
    if(write.tree(tree) != write.tree(fulltree)) {
      
      # Get descendants of each node in pruned tree:
      Descendants <- lapply(as.list(NodeNumbers), function(x) tree$tip.label[strap::FindDescendants(x, tree = tree)])
      
      # Get corresponding ancestral node in full tree:
      Ancestors <- unlist(lapply(Descendants, function(x) Claddis::FindAncestor(descs = x, tree = fulltree)))
      
      # If pruned tree is identical to full tree (not pruned at all):
    } else {
      
      # Set ancestors as node numbers:
      Ancestors <- NodeNumbers
      
    }
    
    # Output matrix matching node numbers of pruned tree to full tree:
    return(matrix(c(NodeNumbers, Ancestors), ncol = 2, dimnames = list(c(), c("PrunedNode", "FullNode"))))
    
  }
  
  # Get pruned node to full node for each unique tree:
  NodeMapsList <- lapply(UniqueTrees, MapPrunedTreeNodesToFullTreeNodes, fulltree = Tree)
  
  # Build out for all trees (adds in any duplicated trees):
  NodeMapsList <- NodeMapsList[match(NewickStrings, UniqueNewickStrings)]
  
  # Add node maps to data list:
  for(i in 1:length(DataAsList)) DataAsList[[i]]$NodeMaps <- NodeMapsList[[i]]
  
  # Get number of tips in tree:
  NTips <- ape::Ntip(Tree)
  
  # Get number of nodes in tree:
  NNodes <- ape::Nnode(Tree)
  
  # Get all node names and numbers:
  Nodes <- c(rownames(OriginalMatrix), (NTips + 1):(NTips + NNodes))
  
  # Renumber nodes of ancestral states:
  DataAsList <- lapply(DataAsList, function(x) { names(x$AncestralStates)[match(as.character(x$NodeMaps[, "PrunedNode"]), names(x$AncestralStates))] <- as.character(x$NodeMaps[, "FullNode"]); return(x) })
  
  # Collapse down to an ancestral state matrix ready for output:
  AncestralStateMatrix <- do.call(cbind, lapply(DataAsList, function(x) {x$AncestralStates <- x$AncestralStates[Nodes]; names(x$AncestralStates) <- Nodes; return(x$AncestralStates)}))
  
  # Isolate estimated tip values:
  TipMatrix <- AncestralStateMatrix[rownames(OriginalMatrix), ]
  
  # If there are any missing values:
  if(any(is.na(TipMatrix))) {
    
    # Isolate missing values:
    MissingTipStates <- which(is.na(TipMatrix))
    
    # Replace missing values with original (unmodified) input values:
    TipMatrix[MissingTipStates] <- OriginalMatrix[MissingTipStates]
    
    # Add tip values back into full output:
    AncestralStateMatrix[rownames(OriginalMatrix), ] <- TipMatrix
    
  }
  
  # Get column (character) count for each matrix block:
  MatrixColumns <- unlist(lapply(lapply(RawCladisticMatrix[2:length(RawCladisticMatrix)], '[[', "Matrix"), ncol))
  
  # For each matrix block:
  for(i in 1:length(MatrixColumns)) {
    
    # Insert portion of ancestral state estimate into block:
    RawCladisticMatrix[[(i + 1)]]$Matrix <- AncestralStateMatrix[, 1:MatrixColumns[i], drop = FALSE]
    
    # Remove that portion from the block:
    AncestralStateMatrix <- AncestralStateMatrix[, -(1:MatrixColumns[i]), drop = FALSE]
    
  }
  
  # Overwrite ancestral state output with updated raw input:
  AncestralStateMatrix <- RawCladisticMatrix
  
  # Add tree to output:
  AncestralStateMatrix$Topper$Tree <- Tree
  
}
beep(3)

# Confirm everything looks OK
str(AncestralStateMatrix)

# Check states for tips and nodes
AncestralStateMatrix$Matrix_1$Matrix[c(1:5, 726:731), 1:10]

# Save processed data
# mode.anc <- AncestralStateMatrix; save(mode.anc, file = "mode.anc"); load("mode.anc")
# constant.anc <- AncestralStateMatrix; save(constant.anc, file = "constant.anc"); load("constant.anc")
# raw.anc <- AncestralStateMatrix; save(raw.anc, file = "raw.anc"); load("raw.anc")
# morph.anc <- AncestralStateMatrix; save(morph.anc, file = "morph.anc"); load("morph.anc")

# Convert ancestral state matrix to csv for viewing outside R
write.csv(mode.anc$Matrix_1$Matrix, file = "modeanc.csv")
write.csv(constant.anc$Matrix_1$Matrix, file = "constantanc.csv")
write.csv(raw.anc$Matrix_1$Matrix, file = "rawanc.csv")
write.csv(morph.anc$Matrix_1$Matrix, file = "morphanc.csv")

# Function does following:
#   * Appends tree to first element [[1]].
#   * Appends node reconstructions to the data blocks in $Matrix[[2-n]].
# This means it is sensible to run the function on individual blocks of data
# (even portions of a block) in bite-size chunks, then combining the blocks back
# together for analyses of disparity.
