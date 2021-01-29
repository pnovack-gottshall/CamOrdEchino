## Modified DiscreteCharacterRate() function from 'Claddis' package that accepts
## pre-inferred ancestral states output from Claddis::AncStateEstMatrix(). Code
## by Graeme Lloyd and taken from
## https://github.com/graemetlloyd/Claddis/blob/master/R/AncStateEstMatrix.R.
## Downloaded 7/2/2020. Thanks, Graeme!

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 7/2/2020, v.
# 0.4.1. Future users will need to either download the archived version of the
# package from GitHub or alter the code accordingly to use the current argument
# and function names.

DiscreteCharacterRate2 <- function(tree, CladisticMatrix, TimeBins, 
                                   BranchPartitionsToTest = NULL,
                                   CharacterPartitionsToTest = NULL,
                                   CladePartitionsToTest = NULL,
                                   TimeBinPartitionsToTest = NULL,
                                   ChangeTimes = "random",
                                   LikelihoodTest = "AIC",
                                   Alpha = 0.01,
                                   MultipleComparisonCorrection = "BenjaminiHochberg",
                                   PolymorphismState = "missing",
                                   UncertaintyState = "missing",
                                   InapplicableState = "missing",
                                   TimeBinApproach = "Lloyd",
                                   EnsureAllWeightsAreIntegers = FALSE) {
  
  # Check for step matrices and stop and warn user if found:
  if(is.list(CladisticMatrix$Topper$StepMatrices)) stop("Function cannot currently deal with step matrices.")
  
  # Check tree has branch lengths:
  if(is.null(tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")
  
  # Check tree has root age:
  if(is.null(tree$root.time)) stop("Tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")
  
  # Check ChangeTimes is correctly formatted or stop and warn user:
  if(length(setdiff(ChangeTimes, c("midpoint", "spaced", "random"))) > 0) stop("ChangeTimes must be one of \"midpoint\", \"spaced\", or \"random\".")
  
  # Check MultipleComparisonCorrection is correctly formatted or stop and warn user:
  if(length(setdiff(MultipleComparisonCorrection, c("BenjaminiHochberg", "Bonferroni"))) > 0) stop("MultipleComparisonCorrection must be one of \"BenjaminiHochberg\" or \"Bonferroni\".")
  
  # Check PolymorphismState is correctly formatted or stop and warn user:
  if(length(setdiff(PolymorphismState, c("missing", "random"))) > 0) stop("PolymorphismState must be one of \"missing\" or \"random\".")
  
  # Check UncertaintyState is correctly formatted or stop and warn user:
  if(length(setdiff(UncertaintyState, c("missing", "random"))) > 0) stop("UncertaintyState must be one of \"missing\" or \"random\".")
  
  # Check InapplicableState is correctly formatted or stop and warn user:
  if(length(setdiff(InapplicableState, c("missing"))) > 0) stop("InapplicableState must be \"missing\".")
  
  # Check TimeBinApproach is correctly formatted or stop and warn user:
  if(length(setdiff(TimeBinApproach, c("Close", "Lloyd"))) > 0) stop("TimeBinApproach must be one of \"Close\" or \"Lloyd\".")
  
  # Check partitions are not all NULL values:
  if(is.null(BranchPartitionsToTest) && is.null(CharacterPartitionsToTest) && is.null(CladePartitionsToTest) && is.null(TimeBinPartitionsToTest)) stop("No partitions are requested. Set at least one of BranchPartitionsToTest, CharacterPartitionsToTest, CladePartitionsToTest, or TimeBinPartitionsToTest to a list of appropriate values. Type \"?DiscreteCharacterRate\" for help.")
  
  # Get internal node numbers:
  InternalNodeNumbers <- 1:ape::Nnode(tree) + ape::Ntip(tree)
  
  # Get edge numbers:
  EdgeNumbers <- 1:nrow(tree$edge)
  
  # Get character numbers:
  CharacterNumbers <- 1:sum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)))
  
  # Ensure time bins are in correct order:
  TimeBins <- sort(unique(TimeBins), decreasing = TRUE)
  
  # Find the Time bin midpoints:
  TimeBinMidpoints <- (TimeBins[2:length(TimeBins)] + TimeBins[1:(length(TimeBins) - 1)]) / 2
  
  # Get the numbers for each time bins:
  TimeBinNumbers <- 1:length(TimeBinMidpoints)
  
  # Subfunction to ensure partitions are formatted correctly:
  PartitionFormatter <- function(PartitionsToTest, ValidValues, PartitionName) {
    
    # Check partitions are in the form of a list of lists:
    if(!all(c(all(unlist(lapply(PartitionsToTest, is.list))), is.list(PartitionsToTest)))) stop(paste(PartitionName, " must be in the form of a list of lists.", sep = ""))
    
    # Get a vector of any non-present valid values:
    NonPresentValues <- setdiff(unique(unlist(PartitionsToTest)), ValidValues)
    
    # Check valid values have been used and if not stop and warn user:
    if(length(NonPresentValues) > 0) stop(paste(PartitionName, "Partitions to test must be defined using the valid range of values (", paste(range(ValidValues), collapse = " to "), ") only.", sep = ""))
    
    # Check partitions never overlap and stop and warn user if they do:
    Check <- lapply(PartitionsToTest, function(x) if(any(duplicated(sort(unlist(x))))) stop(paste("Each partition of ", PartitionName, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))
    
    # Subfunction to ad the missing partition (if exists):
    AddMissingPartitions <- function(x, ValidValues) {
      
      # Define any missing values:
      MissingValues <- setdiff(ValidValues, unlist(x))
      
      # If there are missing values add them to list at end:
      if(length(MissingValues) > 0) x[[(length(x) + 1)]] <- MissingValues
      
      # Return x:
      return(x)
      
    }
    
    # Add in missing partitions (if any):
    PartitionsToTest <- lapply(PartitionsToTest, AddMissingPartitions, ValidValues = ValidValues)
    
    # Check partitions are all at least two in size or else no comparison can be made:
    if(any(unlist(lapply(PartitionsToTest, length)) == 1) && LikelihoodTest == "LRT") stop("Partitions must divide the available data into at least two parts if performing likelihood ratio tests.")
    
    # Return formatted partitions to test:
    return(PartitionsToTest)
    
  }
  
  # If performing branch partition test(s) check and reformat branch partitions:
  if(!is.null(BranchPartitionsToTest)) BranchPartitionsToTest <- PartitionFormatter(PartitionsToTest = BranchPartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "BranchPartitionsToTest")
  
  # If performing character partition test(s) check and reformat character partitions:
  if(!is.null(CharacterPartitionsToTest)) CharacterPartitionsToTest <- PartitionFormatter(PartitionsToTest = CharacterPartitionsToTest, ValidValues = CharacterNumbers, PartitionName = "CharacterPartitionsToTest")
  
  # If performing clade partition test(s)
  if(!is.null(CladePartitionsToTest)) {
    
    # Convert clade partitions to edge partitions:
    CladePartitionsToTest <- lapply(CladePartitionsToTest, lapply, GetDescendantEdges, tree = tree)
    
    # Check and reformat clade partitions:
    CladePartitionsToTest <- PartitionFormatter(PartitionsToTest = CladePartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "CladePartitionsToTest")
    
  }
  
  # If performing time bin partition test(s) check and reformat time bin partitions:
  if(!is.null(TimeBinPartitionsToTest)) TimeBinPartitionsToTest <- PartitionFormatter(PartitionsToTest = TimeBinPartitionsToTest, ValidValues = TimeBinNumbers, PartitionName = "TimeBinPartitionsToTest")
  
  # Check LikelihoodTest is correctly formatted or stop and warn user:
  if(length(setdiff(LikelihoodTest, c("AIC", "LRT"))) > 0) stop("LikelihoodTest must be one of \"AIC\" or \"LRT\".")
  
  # Subfunction to calculate maximum likelihood p value:
  GetMaximumLikelihoodPValue <- function(MeanRate, SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
    # Set maximum likelihood numerator:
    MaximumLikelihoodNumerator <- MeanRate
    
    # Set maximum likelihood denominator:
    MaximumLikelihoodDenominator <- SampledRates
    
    # Get log numerator:
    LogNumerator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodNumerator * SampledCompleteness * SampledTime)))
    
    # Get log denominator:
    LogDenominator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodDenominator * SampledCompleteness * SampledTime)))
    
    # Get test statistic:
    TestStatistic <- -2 * (LogNumerator - LogDenominator)
    
    # Calculate position of test statistic in chi-square distribution to get probability:
    PValue <- pchisq(TestStatistic, length(SampledRates) - 1, lower.tail = FALSE)
    
    # Output probability for later alpha comparison:
    return(PValue)
    
  }
  
  # Subfunction to calculate AIC:
  GetAIC <- function(SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
    # Get log maximum likelihood estimate:
    LogMLE <- sum(log(dpois(x = round(SampledChanges), lambda = SampledRates * SampledCompleteness * SampledTime)))
    
    # Calculate AIC:
    AIC <- (2 * length(SampledRates)) - (2 * LogMLE)
    
    # Return AIC:
    return(AIC)
    
  }
  
  # Subfunction to calculate AIC from partition (with columns labelled Partition, Rate, Completeness, Duration):
  GetAICFromPartition <- function(Partition, AICc = FALSE) {
    
    # Get log maximum likelihood estimate:
    LogMLE <- sum(log(dpois(round(Partition[, "Changes"]), Partition[, "Rate"] * Partition[, "Completeness"] * Partition[, "Duration"])))
    
    # Get k (number of parameters) term:
    k <- max(Partition[, "Partition"])
    
    # Calculate AIC:
    AIC <- (2 * k) - (2 * LogMLE)
    
    # If AICc is desired then calculate this and overwrite AIC with it:
    if(AICc) AIC <- AIC + (((2 * (k ^ 2)) + (2 * k)) / (nrow(Partition) - k - 1))
    
    # Return AIC:
    return(AIC)
    
  }
  
  # Get ages for each (tip and internal) node:
  NodeAges <- GetNodeAges(tree)
  
  # Get branch ages (from and to):
  BranchAges <- unname(cbind(NodeAges[as.character(tree$edge[, 1])], NodeAges[as.character(tree$edge[, 2])]))
  
  # Build edge list from node numbers (from-to) for each branch:
  EdgeList <- lapply(apply(tree$edge, 1, list), function(x) {names(x) <- "NodeNumberFromTo"; return(x)})
  
  # Add node ages to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$NodeAgeFromTo <- BranchAges[i, ]
  
  # Add node ages (from-to) to each edge in list:
  EdgeList <- lapply(EdgeList, function(x) {x$BranchDuration <- x$NodeAgeFromTo[1] - x$NodeAgeFromTo[2]; return(x)})
  
  # Get vector of branch types:
  BranchTypes <- gsub("0", "Internal", gsub("1", "Terminal", as.numeric(tree$edge[, 2] <= ape::Ntip(tree))))
  
  # Add branch type to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$BranchType <- BranchTypes[i]
  
  # Find descendant edges for each internal node:
  DescendantEdgesForEachInternalNode <- lapply(as.list(InternalNodeNumbers), GetDescendantEdges, tree = tree)
  
  # Get ancestral character states:
  # *** MODIFIED TO ACCEPT PREVIOUSLY ESTIMATED ANCESTRAL MATRIX
  AncestralStates <- CladisticMatrix
  
  # Build single matrix of all states in tip label then node number order:
  AllStates <- do.call(cbind, lapply(lapply(AncestralStates[2:length(AncestralStates)], '[[', "Matrix"), function(x) x[c(tree$tip.label, 1:ape::Nnode(tree) + ape::Ntip(tree)), , drop = FALSE]))
  
  # Make vector of ordering of characters:
  Ordering <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Ordering")))
  
  # Make vector of weights of characters:
  Weights <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Weights")))
  
  # Make vector of minimum values:
  MinVals <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MinVals")))
  
  # Make vector of maximum values:
  MaxVals <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MaxVals")))
  
  # Find positions in matrix with polymorphisms:
  PolymorphismPositions <- grep("&", AllStates)
  
  # Find positions in matrix with uncertainties:
  UncertaintyPositions <- grep("/", AllStates)
  
  # Find positions in matrix with inapplicables:
  InapplicablePositions <- which(AllStates == "")
  
  # If polymorphisms were found:
  if(length(PolymorphismPositions) > 0) {
    
    # If replacing polymorphsims with missing do so:
    if(PolymorphismState == "missing") AllStates[PolymorphismPositions] <- NA
    
    # If replacing polymorphisms with random values draw and replace:
    if(PolymorphismState == "random") AllStates[PolymorphismPositions] <- unlist(lapply(strsplit(AllStates[PolymorphismPositions], "&"), sample, size = 1))
    
  }
  
  # If uncertainties were found:
  if(length(UncertaintyPositions) > 0) {
    
    # If replacing uncertainties with missing do so:
    if(UncertaintyState == "missing") AllStates[UncertaintyPositions] <- NA
    
    # If replacing uncertainties with random values draw and replace:
    if(UncertaintyState == "random") AllStates[UncertaintyPositions] <- unlist(lapply(strsplit(AllStates[UncertaintyPositions], "/"), sample, size = 1))
    
  }
  
  # If inapplicable states were found:
  if(length(InapplicablePositions) > 0) {
    
    # If replacing inapplicables with missing do so:
    if(InapplicableState == "missing") AllStates[InapplicablePositions] <- NA
    
  }
  
  # Set default converted continuous characters to FALSE:
  ContinuousCharactersConverted <- FALSE
  
  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if(any(Ordering == "cont")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Set default converted continuous characters to TRUE:
    ContinuousCharactersConverted <- TRUE
    
    # Find out which characters are continuous:
    ContinuousCharactersFound <- which(Ordering == "cont")
    
    # Rescale continous characters as zero to one values:
    ListOfContinuousValuesRescaledZeroToOne <- lapply(lapply(lapply(apply(AllStates[, ContinuousCharactersFound, drop = FALSE], 2, list), unlist), as.numeric), function(x) {x <- x - min(sort(x)); x <- x / max(sort(x)); return(x)})
    
    # Now discretize and store these characters (0 to 31 scale):
    AllStates[, ContinuousCharactersFound] <- do.call(cbind, lapply(lapply(lapply(ListOfContinuousValuesRescaledZeroToOne, function(x) as.list(x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    Ordering[ContinuousCharactersFound] <- "ord"
    
    # Convert weights to 1/31:
    Weights[ContinuousCharactersFound] <- 1 / 31
    
    # Set minimum value to zero:
    MinVals[ContinuousCharactersFound] <- 0
    
    # Set maximum value to 31:
    MaxVals[ContinuousCharactersFound] <- 31
    
  }
  
  # If EnsureAllWeightsAreIntegers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if(EnsureAllWeightsAreIntegers) while(is.character(all.equal(sum(Weights %% 1), 0))) Weights <- (1 / (Weights %% 1)[(Weights %% 1) > 0])[1] * Weights
  
  # Add from-to node states for each character to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$CharacterStatesFromTo <- matrix(AllStates[x$NodeNumberFromTo, , drop = FALSE], nrow = 2, dimnames = list(c("From", "To"))); return(x)})
  
  # Subfunction to define character changes:
  BuildChangesMatrix <- function(x) {
    
    # Find only comparable characters (those scored for both from and to states):
    ComparableCharacters <- which(apply(!apply(x$CharacterStatesFromTo, 2, is.na), 2, all))
    
    # Isolate comparable ordering:
    ComparableOrdering <- Ordering[ComparableCharacters]
    
    # Isolate comparable weights:
    ComparableWeights <- Weights[ComparableCharacters]
    
    # Isolate only characters that actually differ (change):
    CharacterDifferences <- which(x$CharacterStatesFromTo["From", ComparableCharacters] != x$CharacterStatesFromTo["To", ComparableCharacters])
    
    # Build character change matrix:
    CharacterChanges <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("Character", "From", "To", "Steps", "Weight")))
    
    # If characters change then make a matrix from them:
    if(length(CharacterDifferences) > 0) CharacterChanges <- rbind(CharacterChanges, cbind(as.numeric(ComparableCharacters[CharacterDifferences]), as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]), as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]), ifelse(ComparableOrdering[CharacterDifferences] == "unord", 1, abs(as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]) - as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]))), ComparableWeights[CharacterDifferences]))
    
    # Store character changes as new sublist for x:
    x$CharacterChanges <- CharacterChanges
    
    # Store comparable characters as new sublist of x:
    x$ComparableCharacters <- ComparableCharacters
    
    # Return x:
    return(x)
    
  }
  
  # Get character changes and comparable characters and add to edge list:
  EdgeList <- lapply(EdgeList, BuildChangesMatrix)
  
  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Subfunction to add change times to character changes:
    AddChangeTimes <- function(x, ChangeTimes) {
      
      # Isolate character changes:
      CharacterChanges <- x$CharacterChanges
      
      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if(any(CharacterChanges[, "Steps"] > 1)) {
        
        # Get multistep character changes:
        MultiStepCharacters <- which(CharacterChanges[, "Steps"] > 1)
        
        # For each multistep character change:
        for(i in rev(MultiStepCharacters)) {
          
          # Isolate other rows:
          OtherRowNumbers <- setdiff(1:nrow(CharacterChanges), i)
          
          # Get unpacked changes (X:Y, e.g., 0:2 would become 0 1 2):
          UnpackedChanges <- CharacterChanges[i, "From"]:CharacterChanges[i, "To"]
          
          # Update character changes with multistep changes unpacked:
          CharacterChanges <- rbind(CharacterChanges[OtherRowNumbers, ], unname(cbind(rep(CharacterChanges[i, "Character"], length.out = length(UnpackedChanges) - 1), UnpackedChanges[1:(length(UnpackedChanges) - 1)], UnpackedChanges[2:length(UnpackedChanges)], rep(1, length.out = length(UnpackedChanges) - 1), rep(CharacterChanges[i, "Weight"], length.out = length(UnpackedChanges) - 1))))
          
        }
        
        # Resort character changes by character number:
        CharacterChanges <- CharacterChanges[order(CharacterChanges[, "Character"]), ]
        
      }
      
      # If using midpoint option set character change times as midpoint of branch:
      if(ChangeTimes == "midpoint") CharacterChanges <- cbind(CharacterChanges, rep(x$NodeAgeFromTo[1] - (x$BranchDuration / 2), length.out = nrow(CharacterChanges)))
      
      # If using spaced then set character change times as equally spaced along branch:
      if(ChangeTimes == "spaced") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - (seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1)[1:nrow(CharacterChanges)] + (diff(seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1))[1] / 2)))
      
      # If using random then set character change times as random draws from a uniform distribution:
      if(ChangeTimes == "random") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - stats::runif(n = nrow(CharacterChanges), min = 0, max = x$BranchDuration))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Time"
      
      # Subfunction to re-sort character change times so they occur in correct order:
      SortChangeTimes <- function(CharacterChanges) {
        
        # Sort change time for each character from oldest (first) to youngest (last) and store it:
        CharacterChanges[, "Time"] <- unname(unlist(lapply(as.list(unique(CharacterChanges[, "Character"])), function(x) sort(CharacterChanges[which(CharacterChanges[, "Character"] == x), "Time"], decreasing = TRUE))))
        
        # Return sorted character changes:
        return(CharacterChanges)
        
      }
      
      # Re-sort character change times so they occur in correct order:
      CharacterChanges <- SortChangeTimes(CharacterChanges)
      
      # Add bin for character change as last column:
      CharacterChanges <- cbind(CharacterChanges, unlist(lapply(as.list(CharacterChanges[, "Time"]), function(x) max(which(x <= TimeBins)))))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Bin"
      
      # Overwrite character changes with new version with changes added:
      x$CharacterChanges <- CharacterChanges
      
      # Return x:
      return(x)
      
    }
    
    # Add character change times to edge list:
    EdgeList <- lapply(EdgeList, AddChangeTimes, ChangeTimes = ChangeTimes)
    
  }
  
  # Subfunction to get edge sections in time bins:
  EdgeSectionsInBins <- function(x, TimeBins = TimeBins) {
    
    # Set first appearance datum of edge:
    FAD <- x$NodeAgeFromTo[1]
    
    # Set last appearance datum of edge:
    LAD <- x$NodeAgeFromTo[2]
    
    # Get any time bin boundaries crossed (an be empty if none are):
    BoundariesCrossed <- TimeBins[2:(length(TimeBins) - 1)][intersect(which(TimeBins[2:(length(TimeBins) - 1)] > LAD), which(TimeBins[2:(length(TimeBins) - 1)] < FAD))]
    
    # If boundaries are crossed:
    if(length(BoundariesCrossed) > 0) {
      
      # Break up branch into binned sections as vector of FADs:
      FAD <- c(FAD, BoundariesCrossed)
      
      # Break up branch into binned sections as vector of LADs:
      LAD <- c(BoundariesCrossed, LAD)
      
    }
    
    # Build matrix of branch sections with FADs and LADs:
    BranchSections <- rbind(FAD, LAD)
    
    # Add bin number present in to column names:
    colnames(BranchSections) <- unlist(lapply(lapply(lapply(as.list(BranchSections["FAD", ]), '<=', TimeBins), which), max))
    
    # Add new list section for branch (edge) sections binned by time:
    x$BinnedEdgeSections <- BranchSections
    
    # Return output:
    return(x)
    
  }
  
  # Get edge sections in time bins:
  EdgeList <- lapply(EdgeList, EdgeSectionsInBins, TimeBins = TimeBins)
  
  # Add binned branch durations to edge list:
  EdgeList <- lapply(EdgeList, function(x) {BranchDurations <- rep(0, length(TimeBins) - 1); BranchDurations[as.numeric(colnames(x$BinnedEdgeSections))] <- abs(apply(x$BinnedEdgeSections, 2, diff)); x$BinnedBranchDurations <- BranchDurations; return(x)})
  
  # Add proportional binned branch lengths to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$ProportionalBinnedEdgeDurations <- x$BinnedBranchDurations / sum(x$BinnedBranchDurations); return(x)})
  
  # Start to build matrix of all changes with list of character changes:
  AllChanges <- lapply(EdgeList, function(x) x$CharacterChanges)
  
  # Add edge number to each matrix of character changes:
  for(i in 1:length(AllChanges)) AllChanges[[i]] <- cbind(rep(i, times = nrow(AllChanges[[i]])), AllChanges[[i]])
  
  # Combine all changes into a single matrix:
  AllChanges <- do.call(rbind, lapply(AllChanges, function(x) {colnames(x)[1] <- "Edge"; x}))
  
  # Remove silly rownames from all changes:
  rownames(AllChanges) <- NULL
  
  # If doing some kind of edge test (branch or clade):
  if(!is.null(BranchPartitionsToTest) || !is.null(CladePartitionsToTest)) {
    
    # Get (weighted) number of changes on each edge:
    EdgeChanges <- unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))
    
    # Get completeness for each edge:
    EdgeCompleteness <- unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]) / sum(Weights)))
    
    # Get duration of each edge:
    EdgeDurations <- unlist(lapply(EdgeList, function(x) x$BranchDuration))
    
    # Set global rate:
    GlobalRate <- sum(EdgeChanges) / sum(EdgeCompleteness * EdgeDurations)
    
    # If performing branch partition tests:
    if(!is.null(BranchPartitionsToTest)) {
      
      # If using Likelihood Ratio Test:
      if(LikelihoodTest == "LRT") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        PartitionedData <- lapply(BranchPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
        
        # Add sampled rate to paritioned data matrices:
        PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
        
        # Get LRT p-values and combine output as edge test results:
        BranchPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
        
      }
      
      # If using AIC:
      if(LikelihoodTest == "AIC") {
        
        # Build partitioned data for AIC:
        PartitionedData <- lapply(BranchPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(tree$edge.length)), Rate = rep(NA, length(tree$edge.length)), Changes = EdgeChanges, Completeness = EdgeCompleteness, Duration = EdgeDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
        
        # Get AIC, AICc and rate results:
        BranchPartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = GetAICFromPartition(x), AICc = GetAICFromPartition(x, AICc = TRUE)))
        
      }
      
      # If not performing branch partition tests:
    } else {
      
      # Make empty branch partition result output:
      BranchPartitionTestResults <- NULL
      
    }
    
    # If performing clade partition tests:
    if(!is.null(CladePartitionsToTest)) {
      
      # If using Likelihood Ratio Test:
      if(LikelihoodTest == "LRT") {
        
        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        PartitionedData <- lapply(CladePartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
        
        # Add sampled rate to paritioned data matrices:
        PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
        
        # Get LRT p-values and combine output as edge test results:
        CladePartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
        
      }
      
      # If using AIC:
      if(LikelihoodTest == "AIC") {
        
        # Build partitioned data for AIC:
        PartitionedData <- lapply(CladePartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(tree$edge.length)), Rate = rep(NA, length(tree$edge.length)), Changes = EdgeChanges, Completeness = EdgeCompleteness, Duration = EdgeDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
        
        # Get AIC, AICc and rate results:
        CladePartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = GetAICFromPartition(x), AICc = GetAICFromPartition(x, AICc = TRUE)))
        
      }
      
      # If not performing clade partition tests:
    } else {
      
      # Make empty clade partition result output:
      CladePartitionTestResults <- NULL
      
    }
    
    # If not doing clade OR branch tests:
  } else {
    
    # Make empty branch partition result output:
    BranchPartitionTestResults <- NULL
    
    # Make empty clade partition result output:
    CladePartitionTestResults <- NULL
    
  }
  
  # If performing branch partition tests:
  if(!is.null(CharacterPartitionsToTest)) {
    
    # Get vector of (weighted) changes for each character:
    CharacterChanges <- unlist(lapply(as.list(CharacterNumbers), function(x) {CharacterRows <- which(AllChanges[, "Character"] == x); sum(AllChanges[CharacterRows, "Steps"] * AllChanges[CharacterRows, "Weight"])}))
    
    # Get vector of weighted durations for each character:
    CharacterDurations <- (Weights / sum(Weights)) * sum(tree$edge.length)
    
    # Get vector of completness (opportunity to observe changes) for each character:
    CharacterCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) {CharacterPresence <- rep(0, times = length(CharacterNumbers)); CharacterPresence[x$ComparableCharacters] <- 1; CharacterPresence * x$BranchDuration})), 2, sum) / sum(tree$edge.length)
    
    # Set global rate:
    GlobalRate <- sum(CharacterChanges) / sum(CharacterCompleteness * CharacterDurations)
    
    # If using likelihood ratio test:
    if(LikelihoodTest== "LRT") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      PartitionedData <- lapply(CharacterPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(CharacterChanges[y]), sum(CharacterCompleteness[y] * CharacterDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      CharacterPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
      
    }
    
    # If using AIC:
    if(LikelihoodTest== "AIC") {
      
      # Build partitioned data for AIC:
      PartitionedData <- lapply(CharacterPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, max(CharacterNumbers)), Rate = rep(NA, max(CharacterNumbers)), Changes = CharacterChanges, Completeness = CharacterCompleteness, Duration = CharacterDurations); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / (sum(y[x, "Completeness"]) * sum(y[x, "Duration"])), length(x))))[order(unlist(x))])); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length)))[order(unlist(x))]; y})
      
      # Get AIC, AICc and rate results:
      CharacterPartitionTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = GetAICFromPartition(x), AICc = GetAICFromPartition(x, AICc = TRUE)))
      
    }
    
    # If performing branch partition tests:
  } else {
    
    # Make empty character partition result output:
    CharacterPartitionTestResults <- NULL
    
  }
  
  # If performing time bin partition tests:
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Get weighted number of changes from each time bin:
    TimeBinChanges <- unlist(lapply(as.list(1:(length(TimeBins) - 1)), function(x) {ChangeRows <- AllChanges[, "Bin"] == x; sum(AllChanges[ChangeRows, "Steps"] * AllChanges[ChangeRows, "Weight"])}))
    
    # If using the Close time bin completeness approach get completeness value for each time bin:
    if(TimeBinApproach == "Close") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations * (sum(Weights[x$ComparableCharacters]) / sum(Weights)))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations)), 2, sum)
    
    # If using the Lloyd time bin completeness approach get completeness value for each time bin::
    if(TimeBinApproach == "Lloyd") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(Weights[x$ComparableCharacters], ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(Weights, ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum)
    
    # Get durations of edges in each time bin:
    TimeBinDurations <- apply(do.call(rbind, lapply(EdgeList, function(x) x$BinnedBranchDurations)), 2, sum)
    
    # Set global rate (NB: will differ between Close and Lloyd approaches, but Lloyd approach will match edge or character global rate):
    GlobalRate <- sum(TimeBinChanges) / sum(TimeBinCompleteness * TimeBinDurations)
    
    # If using Likelihood Ratio Test to compare partitions:
    if(LikelihoodTest == "LRT") {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later) for LRT:
      PartitionedData <- lapply(TimeBinPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(TimeBinChanges[y]), sum(TimeBinCompleteness[y] * TimeBinDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices for LRT:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      TimeBinTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
      
    }
    
    # If using AIC to compare partitions:
    if(LikelihoodTest == "AIC") {
      
      # Build partitioned data for AIC:
      PartitionedData <- lapply(TimeBinPartitionsToTest, function(x) {y <- cbind(Partition = rep(NA, length(TimeBinChanges)), Rate = rep(NA, length(TimeBinChanges)), Changes = TimeBinChanges, Completeness = TimeBinCompleteness * TimeBinDurations, Duration = rep(1, length(TimeBinChanges))); y[, "Rate"] <- as.numeric(gsub(NaN, 0, unlist(lapply(x, function(x) rep(sum(y[x, "Changes"]) / sum(y[x, "Completeness"]), length(x)))))); y[, "Partition"] <- rep(1:length(x), unlist(lapply(x, length))); y})
      
      # Get AIC, AICc and rate results:
      TimeBinTestResults <- lapply(PartitionedData, function(x) list(Rates = unname(unlist(lapply(as.list(unique(x[, "Partition"])), function(y) x[x[, "Partition"] == y, "Rate"][1]))), AIC = GetAICFromPartition(x), AICc = GetAICFromPartition(x, AICc = TRUE)))
      
    }
    
    # If not performing time bin partition tests:
  } else {
    
    # Make empty time bin partition result output:
    TimeBinTestResults <- NULL
    
  }
  
  # Set global rate for output:
  GlobalRate <- sum(unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))) / sum(unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]) / sum(Weights))) * unlist(lapply(EdgeList, function(x) x$BranchDuration)))
  
  # If performing Likelihood Ratio Test:
  if(LikelihoodTest == "LRT") {
    
    # Subfunction to calculate adjusted alphas for multiple comparison corrections:
    AddMultipleComparisonCorrectionCutoffs <- function(TestResults, Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection) {
      
      # Get number of comparisons performed:
      NComparisons <- length(TestResults)
      
      # If using the Benjamini-Hochberg false discovery rate approach:
      if(MultipleComparisonCorrection == "BenjaminiHochberg") {
        
        # Set cutoff values:
        CutoffValues <- ((1:NComparisons) / NComparisons) * Alpha
        
        # Get actual p-values found:
        PValues <- unlist(lapply(TestResults, '[[', "PValue"))
        
        # Order cutoffs by p-value rank:
        CutoffValues <- CutoffValues[rank(PValues, ties.method = "random")]
        
      }
      
      # If using the Bonferroni correction set cutoff values as alpha over N:
      if(MultipleComparisonCorrection == "Bonferroni") CutoffValues <- Alpha / NComparisons
      
      # Add cutoffs to output:
      for(i in 1:length(TestResults)) TestResults[[i]]$CorrectedAlpha <- CutoffValues[i]
      
      # Return modified test results:
      return(TestResults)
      
    }
    
    # If doing branch partition tests then add multiple comparison alpha cutoffs:
    if(!is.null(BranchPartitionsToTest)) BranchPartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = BranchPartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing character partition tests then add multiple comparison alpha cutoffs:
    if(!is.null(CharacterPartitionsToTest)) CharacterPartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = CharacterPartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing clade partition tests then add multiple comparison alpha cutoffs:
    if(!is.null(CladePartitionsToTest)) CladePartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = CladePartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
    # If doing time bin partition tests then add multiple comparison alpha cutoffs:
    if(!is.null(TimeBinPartitionsToTest)) TimeBinTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = TimeBinTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
    
  }
  
  # Compile output:
  Output <- list(InferredCharacterChanges = AllChanges, IntrinsicCharacterRate = GlobalRate, ContinuousCharactersConvertedToDiscrete = ContinuousCharactersConverted, BranchPartitionResults = BranchPartitionTestResults, CharacterPartitionResults = CharacterPartitionTestResults, CladePartitionResults = CladePartitionTestResults, TimeBinResults = TimeBinTestResults)
  
  # Return output:
  return(Output)
  
}