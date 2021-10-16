## CONVERGENCE ANALYSES (STAYTON 2015) #########################################

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' package. The package 'convevol' assumes Brownian motion ancestral
# reconstructions, which are not appropriate for the discrete characters we are
# using. Instead, we will use Claddis::AncStateEstMatrix(), which wraps around
# phytools::rerootingMethod(), to reconstruct ancestral states using the
# continuous-time Markov chain (Mk) model of discrete state changes. 'convevol'
# also measures distances between tips and nodes using raw Euclidean distance.
# Here (and like in Lloyd, 2018), we use the Wills generalized Euclidean
# distance matrix (with alpha = 0.5 for hierarchical character dependencies,
# c.f., Hopkins and St. John 2018) from 3-DisparityDistances.R

# Stayton's (2015: p. 2144) C1 statistic measures "the proportion of the maximum
# distance between two lineages that has been 'closed' by subsequent evolution"
# with a value of 1 meaning 100% of characters are convergent and 0 meaning 0%
# of characters are convergent. The C2 statistic measures the absolute magnitude
# of convergence, with greater values interpreted as greater convergence. (Note
# that the C2 statistic should not be compared across different data sets, as
# the values are dependent on the data used.)

# *** The code for calculating C1 and C2 is modified from that in Lloyd (2018)
# supplementary information. Thanks, Graeme! ***

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
setwd("~/Manuscripts/CamOrdEchinos/Data files/NA Reformatted")

# Load packages
library(beepr)      # v. 1.3
library(doParallel) # v. 1.0.16
library(Claddis)    # v. 0.6.3
library(ade4)       # v. 1.7-17
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub\n")




## IMPORT FILES ################################################################

# Import output from 3-DisparityDistances.R
load("morph.distances.GED.5")
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")

# Import time trees saved from 2-InferAncestralStates.R (because the original
# cal3 trees from 1-MakeTimeTrees.R contains zero-length branches but they were
# removed from those appended to ancestral states.) Using the trees from raw.anc
# because the trees are identical across x.anc, and this one is the smallest
# object.
load("raw.anc")



## FUNCTIONS ###################################################################

# Function to calculate the D_tip and D_max values used in Stayton (2015) to
# calculate C1 and C2 metrics. The function is modified from that in Lloyd
# (2018) supplementary information, including removing the conversion of
# distance matrix into a Euclidean distance matrix. Thanks, Graeme!
#   dist.matrix = output from Claddis::MorphDistMatrix (or our modified version 
#                       MorphDistMatrix2).
#   tree = tree (in phylo format)
#   taxon.pairs = two-column matrix of taxon tip pairs to calculate statistics 
#                       for.
StaytonConvergence <- function(dist.matrix, tree, taxon.pairs) {

  if (!inherits(tree, "phylo"))
    stop("your tree must be class 'phylo.'")
  
  # Confirm all taxa in taxon.pairs are in 'tree'
  all.taxa <- unique(as.vector(taxon.pairs))
  if (!all(all.taxa %in% tree$tip.label))
    stop("not all taxa are in 'tree'. Check spelling of desired taxa.\n")
  
  # Reformat distance matrix *** MODIFIED: ORIGINAL CONVERTED TO EUCLIDEAN DISTANCE MATRIX
  dist.matrix <- as.matrix(as.dist(dist.matrix))
  
  LineageDistances <- matrix(nrow = 0, ncol = 2)

  for (i in 1:nrow(taxon.pairs)) {
    CurrentMRCA <- Claddis::find_mrca(descendant_names = taxon.pairs[i, ], tree = tree)
    LineageNodes <- CurrentNodes <- c(which(tree$tip.label == taxon.pairs[i, 1]), 
                                      which(tree$tip.label == taxon.pairs[i, 2]))
    AncestorNodes <- c(tree$edge[which(tree$edge[, 2] == CurrentNodes[1]), 1], 
                       tree$edge[which(tree$edge[, 2] == CurrentNodes[2]), 1])
    while (length(which(AncestorNodes == CurrentMRCA)) < 2) {
      CurrentNodes[which(AncestorNodes != CurrentMRCA)] <-
        AncestorNodes[which(AncestorNodes != CurrentMRCA)]
      AncestorNodes <- c(tree$edge[which(tree$edge[, 2] == CurrentNodes[1]), 1], 
                         tree$edge[which(tree$edge[, 2] == CurrentNodes[2]), 1])
      LineageNodes <- c(LineageNodes, CurrentNodes)
    }
    LineageNodes <-
      c(taxon.pairs[i, ], setdiff(unique(LineageNodes), 1:Ntip(tree)))
    CurrentRows <- match(LineageNodes, rownames(dist.matrix))
    D_max <- max(dist.matrix[CurrentRows, CurrentRows])
    D_tip <- dist.matrix[CurrentRows[1], CurrentRows[2]]
    LineageDistances <- rbind(LineageDistances, c(D_tip, D_max))
  }
  
  C1 <- 1 - (LineageDistances[, 1] / LineageDistances[, 2])
  C2 <- LineageDistances[, 2] - LineageDistances[, 1]
  out <- as.data.frame(cbind(taxon.pairs, LineageDistances, C1, C2))
  colnames(out) <- c("Taxon1", "Taxon2", "Dtip", "Dmax", "C1", "C2")
  out[, 3:6] <- sapply(out[, 3:6], as.numeric)
  
  return(out)
}

# Function to extract branching history (MRCA and all subsequent nodes in
# phylospace) for pair of tips. (Clunkier version than used in Tristan Stayton's
# convevol::convrat() function, but consistent with code below.)
#   tree = tree (in phylo format)
#   t1 = character name of tip genus 1
#   t2 = character name of tip genus 2
ancestral.lineages <- function(tree, t1, t2) {
  # Find ancestral edges for pair of taxa until united in MRCA:
  wh.tip <- which(tree$tip.label == t1 | tree$tip.label == t2)
  wh <- wh.tip
  repeat {
    wh.nodes <- tree$edge[which(tree$edge[, 2] %in% wh), 1]
    wh <- c(wh, wh.nodes)
    if (any(duplicated(wh.nodes)))
      break
  }
  # Trim off any pre-MRCA branchs (caused when one taxon has shorter path)
  # Using because faster, in this case, than: ape::getMRCA(tree, c(t1, t2))
  wh.dupl <- wh.nodes[which(duplicated(wh.nodes))]
  branching.history <- tree$edge[which(tree$edge[, 2] %in% wh), ]
  branching.history <- branching.history[which(branching.history[, 1] >= wh.dupl), ]
  wh <- c(wh.tip, branching.history[, 1])
  out <- branching.history
  return(out)
}

# Function to use output from ancestral.lineages to calculate phylogenetic
# distance traveled by two tips since MRCA. (Clunkier version than used in
# Tristan Stayton's convevol::convrat() function, but allows use of
# pre-calculated distance matrix.)
#   branching.history = matrix, output from ancestral.lineages() above. First 
#                       cell must be the MRCA.
#   dist.matrix = output from Claddis::MorphDistMatrix (or our modified version 
#                       MorphDistMatrix2).
branch.distance <- function(branching.history, dist.matrix) {
  dist.matrix <- as.matrix(as.dist(dist.matrix))
  MRCA <- branching.history[1]
  break.point <- which(branching.history[, 1] == MRCA)
  hist1 <- matrix(branching.history[1:diff(break.point), ], ncol = 2)
  hist2 <-
    matrix(branching.history[break.point[2]:nrow(branching.history), ], ncol = 2)
  dist1 <- dist2 <- 0
  for (r in 1:nrow(hist1)) {
    dist1 <- dist1 + dist.matrix[hist1[r, 1], hist1[r, 2]]
  }
  for (r in 1:nrow(hist2)) {
    dist2 <- dist2 + dist.matrix[hist2[r, 1], hist2[r, 2]]
  }
  pair.dist <- sum(dist1, dist2)
  return(pair.dist)
}





## CALCULATE CONVERGENCE STATISTICS C1, C2, and C3 OF STAYTON (2015) ###########

# Get names for all pairs of tips for each tree. Although the order of
# $tip.label is different in each time tree, no need to build separately for
# each tree because the tip taxa (i.e., the rownames of $ranges.used) are
# identical across trees.
tree <- raw.anc[[1]]$topper$tree
taxa <- rownames(tree$ranges.used)
taxon.pairs <- matrix(nrow = 0, ncol = 2)
for (i in 1:(Ntip(tree) - 1)) {
  taxon.pairs <- rbind(taxon.pairs, cbind(rep(taxa[i], (Ntip(tree) - 1)), 
                                          taxa[taxa != taxa[i]]))
}
beep()
# save(taxon.pairs, file = "taxon.pairs")
load("taxon.pairs")



## Calculate C1 & C2 convergence statistics for these pairs of taxa.

# Running in parallel with load-balancing because different tree structures may
# have differences in calculation time for statistics. Note that, because the
# same tip character states are used across trees, resulting $Dtip will be
# identical across trees. The convergence statistics will vary based on Dmax and
# branching lengths, which depend on the ancestral states and varying branching
# edge.lengths and toplogies. Not using load-balancing because run times are
# approximately equal in duration.

ntrees <- length(morph.distances.GED.5)
morph.conv <- mode.conv <- constant.conv <- raw.conv <- vector("list", ntrees)
CPUs <- detectCores(logical = TRUE)
cl <- makeCluster(CPUs)
registerDoParallel(cl)

# Morphology
(start <- Sys.time())
cat("starting 'morphology' at", format(start, "%H:%M:%S"), "\n")
morph.conv <- foreach(i = 1:ntrees, .inorder = TRUE, .packages = "FD") %dopar% {
  StaytonConvergence(dist.matrix = morph.distances.GED.5[[i]]$distance_matrix,
                     tree = raw.anc[[i]]$topper$tree, taxon.pairs = taxon.pairs)
  }

# Mode
cat("starting 'mode' at", format(Sys.time(), "%H:%M:%S"), "\n")
mode.conv <- foreach(i = 1:ntrees, .inorder = TRUE, .packages = "FD") %dopar% {
  StaytonConvergence(dist.matrix = mode.distances.GED.5[[i]]$distance_matrix,
                     tree = raw.anc[[i]]$topper$tree, taxon.pairs = taxon.pairs)
}

# Constant
cat("starting 'constant' at", format(Sys.time(), "%H:%M:%S"), "\n")
constant.conv <- foreach(i = 1:ntrees, .inorder = TRUE, .packages = "FD") %dopar% {
  StaytonConvergence(dist.matrix = constant.distances.GED.5[[i]]$distance_matrix,
                     tree = raw.anc[[i]]$topper$tree, taxon.pairs = taxon.pairs)
}

# Raw
cat("starting 'raw' at", format(Sys.time(), "%H:%M:%S"), "\n")
raw.conv <- foreach(i = 1:ntrees, .inorder = TRUE, .packages = "FD") %dopar% {
  StaytonConvergence(dist.matrix = raw.distances.GED.5[[i]]$distance_matrix,
                     tree = raw.anc[[i]]$topper$tree, taxon.pairs = taxon.pairs)
}

stopCluster(cl)
Sys.time() - start     # 3.9 hrs (~ 1 hour per data set), on 8-core laptop
# save(morph.conv, file = "morph.conv")
# save(mode.conv, file = "mode.conv")
# save(constant.conv, file = "constant.conv")
# save(raw.conv, file = "raw.conv")
beep(3)



## Calculate branch distances (used for C3) for these pairs of taxa (running
## each tree as a loop, and processing each tree's 'taxon.pairs' in parallel
## because that's the fastest way.)
(start <- Sys.time())
ntrees <- length(morph.distances.GED.5)
sq <- 1:nrow(taxon.pairs)
CPUs <- detectCores(logical = TRUE)
for(t in 1:ntrees){
  # To monitor progress
  cat("processing all data sets for tree", t, "at", 
      format(Sys.time(), "%H:%M:%S"), "\n")
  
  # Store results for each individual tree's taxon.pair statistics
  morph.BD <- mode.BD <- constant.BD <- raw.BD <- rep(NA, nrow(taxon.pairs))

  # Morphology
  cl <- makeCluster(CPUs)
  registerDoParallel(cl)
  morph.BD <- foreach(i = sq, .combine = c) %dopar% {
    BH <- ancestral.lineages(tree = raw.anc[[t]]$topper$tree,
                             t1 = taxon.pairs[i, 1], t2 = taxon.pairs[i, 2])
    morph.BD.pair <- 
      branch.distance(BH, morph.distances.GED.5[[t]]$distance_matrix)
  }
  stopCluster(cl)
  # Calculate C3 ( = C2 / branch distance) and append as new column
  # Append branch-length distances as new column
  morph.conv[[t]]$C3 <- morph.conv[[t]]$C2 / morph.BD
  morph.conv[[t]]$branch.dist <- morph.BD
  
  # Mode
  cl <- makeCluster(CPUs)
  registerDoParallel(cl)
  mode.BD <- foreach(i = sq, .combine = c) %dopar% {
    BH <- ancestral.lineages(tree = raw.anc[[t]]$topper$tree, 
                             t1 = taxon.pairs[i, 1], t2 = taxon.pairs[i, 2])
    mode.BD.pair <- 
      branch.distance(BH, mode.distances.GED.5[[t]]$distance_matrix)
  }
  stopCluster(cl)
  mode.conv[[t]]$C3 <- mode.conv[[t]]$C2 / mode.BD
  mode.conv[[t]]$branch.dist <- mode.BD
  
  # Constant
  cl <- makeCluster(CPUs)
  registerDoParallel(cl)
  constant.BD <- foreach(i = sq, .combine = c) %dopar% {
    BH <- ancestral.lineages(tree = raw.anc[[t]]$topper$tree, 
                             t1 = taxon.pairs[i, 1], t2 = taxon.pairs[i, 2])
    constant.BD.pair <- 
      branch.distance(BH, constant.distances.GED.5[[t]]$distance_matrix)
  }
  stopCluster(cl)
  constant.conv[[t]]$C3 <- constant.conv[[t]]$C2 / constant.BD
  constant.conv[[t]]$branch.dist <- constant.BD
  
  # Raw
  cl <- makeCluster(CPUs)
  registerDoParallel(cl)
  raw.BD <- foreach(i = sq, .combine = c) %dopar% {
    BH <- ancestral.lineages(tree = raw.anc[[t]]$topper$tree, 
                             t1 = taxon.pairs[i, 1], t2 = taxon.pairs[i, 2])
    raw.BD.pair <- 
      branch.distance(BH, raw.distances.GED.5[[t]]$distance_matrix)
  }
  stopCluster(cl)
  raw.conv[[t]]$C3 <- raw.conv[[t]]$C2 / raw.BD
  raw.conv[[t]]$branch.dist <- raw.BD
}
Sys.time() - start     # 120.8 hrs on 8-core laptop

# Save output
# save(mode.conv, file = "mode.conv")
# save(constant.conv, file = "constant.conv")
# save(raw.conv, file = "raw.conv")
# save(morph.conv, file = "morph.conv")

beep(3)




## WHAT TAXONOMIC RANK ARE THE PAIRINGS? #######################################

# In other words, what is the taxonomic distance among convergent pairs? Because
# based on taxonomy and not phylogeny, the taxonomic distance will be identical
# across trees and data sets. Will use the morphological data set to identify
# matches, then append to the prior convergence output for all trees and data
# sets.

# Import ecological data set, which has all taxonomic assignments:
data <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", 
                 header = TRUE, stringsAsFactors = FALSE)
data[1:10, 1:5]
# Modify names for two subgenera treated as genera
data$Genus[which(data$Genus == "Anatifopsis")] <-
  c("Anatiferocystis", "Guichenocarpos")

# Assign rank for morphological data set (which is same for all data sets)
rank <- rep(NA, nrow(morph.conv[[1]]))
for (r in 1:length(rank)) {
  pair <- morph.conv[[1]][r, 1:2]
  taxa <- data[which(data$Genus %in% pair), 1:9]
  shared <- which(apply(taxa, 2, duplicated)[2, ])
  matches <- length(shared)
  shared.ranks <- names(shared)
  if (matches == 1L)
    rank[r] <- shared.ranks
  else {
    for (c in matches:1) {
      if (taxa[1, shared.ranks[c]] == "" |
          taxa[1, shared.ranks[c]] == "UNCERTAIN")
        next
      else {
        rank[r] <- shared.ranks[c]
        break
      }
    }
  }
}
beep()
sort(table(rank))


# Append to earlier output and save
for (t in 1:length(morph.conv)) {
  morph.conv[[t]]$PairRank <- rank
  mode.conv[[t]]$PairRank <- rank
  constant.conv[[t]]$PairRank <- rank
  raw.conv[[t]]$PairRank <- rank
}

# save(mode.conv, file = "mode.conv")
# save(constant.conv, file = "constant.conv")
# save(raw.conv, file = "raw.conv")
# save(morph.conv, file = "morph.conv")




## CALCULATE THE MOST FREQUENT CONVERGENCE STATISTICS ACROSS ALL 50 TREES ######
load("mode.conv")
load("constant.conv")
load("raw.conv")
load("morph.conv")
beep()

## FUNCTION TO CALCULATE THE MODE AMONG OBSERVED STATISTICS
# NAs (and missing "" data) can be omitted (with default to include them, to
# match behavior of mean() and median(). The standard 'mode' solution in case of
# ties (an intermediate value) is nonsensical when dealing with discrete states.
# Here, ties are won based on which state is listed first in the input data,
# which is a pseudorandom way to 'flip a coin' to break a tie.
Mode <- function(x, na.rm = FALSE) {
  if (na.rm)
    ux <- unique(x[!is.na(x) & x != ""])
  else
    ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Calculate median convergence statistics (using median because not normally
# distributed). Because x.conv is a data frame, less easy than using 
# 'apply(simplify2array(x.conv), 1:2, Mode)' directly.

# Remember that the first two and the last one columns are constant across trees
# (and data sets), so only need to calculate colums 3-8.
Mode.mode.conv <- Mode.constant.conv <- Mode.raw.conv <- Mode.morph.conv <- 
  mode.conv[[1]]
sq <- 1:length(mode.conv)

# Convert numeric columns to matrix and calculate Mode
tmp.conv <- lapply(sq, function(sq) as.matrix(mode.conv[[sq]][, 3:8]))
Mode.mode.conv[, 3:8] <- apply(simplify2array(tmp.conv), 1:2, Mode, na.rm = TRUE)

tmp.conv <- lapply(sq, function(sq) as.matrix(constant.conv[[sq]][, 3:8]))
Mode.constant.conv[, 3:8] <- apply(simplify2array(tmp.conv), 1:2, Mode, na.rm = TRUE)

tmp.conv <- lapply(sq, function(sq) as.matrix(raw.conv[[sq]][, 3:8]))
Mode.raw.conv[, 3:8] <- apply(simplify2array(tmp.conv), 1:2, Mode, na.rm = TRUE)

tmp.conv <- lapply(sq, function(sq) as.matrix(morph.conv[[sq]][, 3:8]))
Mode.morph.conv[, 3:8] <- apply(simplify2array(tmp.conv), 1:2, Mode, na.rm = TRUE)
beep(3)

# Save
# save(Mode.morph.conv, file = "Mode.morph.conv")
# save(Mode.mode.conv, file = "Mode.mode.conv")
# save(Mode.constant.conv, file = "Mode.constant.conv")
# save(Mode.raw.conv, file = "Mode.raw.conv")
rm("tmp.conv")
gc()


## A curious outcome explained:

# Stayton (2015) claims C3 should have a maximum value of 1 because it is the
# proportion of C2 convergence accounted for in branch evolution. But 82 taxon
# pairs return C3 values greater than 1. A few are easily explained as rounding
# errors (e.g., C3 = 1.11e-15). Most of the remaining "impossible" values here
# occur as a result of a few taxa, especially pairings with Wellerocystis,
# Canadocystis, Implicaticystis, Apektocrinus, Holocystites, Oklahomacystis,
# Petalocystites, Titanocrinus, and Cryptocrinites. Infinite and "> 1" C3 values
# occur when the tips, MRCA, and intervening lineages all have zero distances.
# In some of these cases, they may not be exactly identical (and may have
# different coordinates in PCoA space) because of pairwise distance handling for
# missing character states. Overall, there 96 such pairings in the mode
# treatment, 480 in constant, and 146 in the raw, accounting for 0.07 - 0.36% of
# the cases, and these can be ignored. There are 0 cases in the morphological
# data set. The higher proportion in the constant is additional confirmation
# these are caused by missing character states.
rounders <- which(Mode.mode.conv$C3 > 1 & Mode.mode.conv$C3 < 1.001)
(Mode.mode.conv$C3[rounders] - 1)

oddballs <- which(Mode.mode.conv$C3 > 1.001 & is.finite(Mode.mode.conv$C3))
length(oddballs)
sort(table(c(Mode.mode.conv$Taxon1[oddballs], Mode.mode.conv$Taxon2[oddballs])))
Mode.mode.conv[oddballs, ]

odder.oddballs <- which(is.infinite(Mode.mode.conv$C3))
sort(table(c(Mode.mode.conv$Taxon1[odder.oddballs], Mode.mode.conv$Taxon2[odder.oddballs])))
Mode.mode.conv[odder.oddballs, ]

length(unique(c(oddballs, odder.oddballs)))
100 * round(length(unique(c(oddballs, odder.oddballs))) / nrow(Mode.mode.conv), 4)




## ANALYSE PATTERNS ############################################################

# Plot C1 histograms:
par(mfrow = c(2,2))
breaks <- seq(from = 0, to = 1, by = 0.05)
h.morph1 <- hist(Mode.morph.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "morphology")
h.mode1 <- hist(Mode.mode.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (mode treatment)")
h.constant1 <- hist(Mode.constant.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (constant treatment)")
h.raw1 <- hist(Mode.raw.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (raw treatment)")

# Plot C2 histograms:
breaks2 <- seq(from = 0, to = 6.25, by = 0.25)
h.morph2 <- hist(Mode.morph.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "morphology")
h.mode2 <- hist(Mode.mode.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (mode treatment)")
h.constant2 <- hist(Mode.constant.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (constant treatment)")
h.raw2 <- hist(Mode.raw.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (raw treatment)")

# Plot C3 histograms (truncated for errors caused by missing data when
# calculating pairwise distance measures):
breaks3 <- seq(from = 0, to = 1, by = 0.05)
# Remove errors caused by rounding and missing data problems (see above)
h.morph3 <- hist(Mode.morph.conv$C3[Mode.morph.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "morphology")
h.mode3 <- hist(Mode.mode.conv$C3[Mode.mode.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (mode treatment)")
h.constant3 <- hist(Mode.constant.conv$C3[Mode.constant.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (constant treatment)")
h.raw3 <- hist(Mode.raw.conv$C3[Mode.raw.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (raw treatment)")

# Plot branch-length-distance histograms:
breaks4 <- seq(from = 0, to = 40, by = 2)
h.morph4 <- hist(Mode.morph.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "morphology")
h.mode4 <- hist(Mode.mode.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (mode treatment)")
h.constant4 <- hist(Mode.constant.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (constant treatment)")
h.raw4 <- hist(na.omit(Mode.raw.conv$branch.dist), breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (raw treatment)")

# Are the convergence statistics correlated between ecological and morphological
# data sets?
cor(Mode.morph.conv$C1, Mode.mode.conv$C1, use = "complete.obs") # r = 0.0413
cor(Mode.morph.conv$C2, Mode.mode.conv$C2)                       # r = 0.0467
wh.inf <- which(is.infinite(Mode.mode.conv$C3))
cor(Mode.morph.conv$C3[-wh.inf], Mode.mode.conv$C3[-wh.inf], use = "complete.obs") # r = 0.031
cor(Mode.morph.conv$branch.dist, Mode.mode.conv$branch.dist)     # r = 0.578
par(mfrow = c(1,1))
plot(Mode.morph.conv$branch.dist, Mode.mode.conv$branch.dist)    # but lots of noise



# pdf(file = "convergence_statistics.pdf")
par(mfrow = c(2,2))
breaks <- seq(from = 0, to = 1, by = 0.05)
h.morph1 <- hist(Mode.morph.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "morphology")
h.mode1 <- hist(Mode.mode.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (mode treatment)")
breaks3 <- seq(from = 0, to = 1, by = 0.05)
h.morph3 <- hist(Mode.morph.conv$C3[Mode.morph.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "morphology")
h.mode3 <- hist(Mode.mode.conv$C3[Mode.mode.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (mode treatment)")
# dev.off()

# Median pairwise convergence:
round(apply(Mode.morph.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(Mode.mode.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(Mode.constant.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(Mode.raw.conv[, 5:8], 2, median, na.rm = TRUE), 2)
# Morphology: median C1 = 0.09, C3 = 0.03, branch.distance = 19.79
# Eco.(mode): median C1 = 0.11, C3 = 0.03, branch.distance = 11.18
# Eco.(con.): median C1 = 0.21, C3 = 0.06, branch.distance =  9.46
# Eco.(raw):  median C1 = 0.14, C3 = 0.04, branch.distance =  7.58

# Are the distributions different?
ks.test(Mode.morph.conv$C1, Mode.mode.conv$C1)     # D = 0.082, p < 2.2e-16
ks.test(Mode.morph.conv$C1, Mode.constant.conv$C1) # D = 0.237, p < 2.2e-16
ks.test(Mode.morph.conv$C1, Mode.raw.conv$C1)      # D = 0.153, p < 2.2e-16

ks.test(Mode.morph.conv$C3, Mode.mode.conv$C3)     # D = 0.098, p < 2.2e-16
ks.test(Mode.morph.conv$C3, Mode.constant.conv$C3) # D = 0.269, p < 2.2e-16
ks.test(Mode.morph.conv$C3, Mode.raw.conv$C3)      # D = 0.203, p < 2.2e-16

ks.test(Mode.morph.conv$branch.dist, Mode.mode.conv$branch.dist)     # D = 0.513, p < 2.2e-16
ks.test(Mode.morph.conv$branch.dist, Mode.constant.conv$branch.dist) # D = 0.614, p < 2.2e-16
ks.test(Mode.morph.conv$branch.dist, Mode.raw.conv$branch.dist)      # D = 0.784, p < 2.2e-16


# Histogram of branch distances for convergent pairings
# pdf(file = "branch_distances.pdf")
par(mfrow = c(1,1))
breaks <-
  pretty(c(Mode.morph.conv$branch.dist, Mode.mode.conv$branch.dist), 10)
hist(Mode.mode.conv$branch.dist, main = "Phylogenetic distances", 
     xlab = "branch distance among taxon pairs", ylab = "Density", 
     breaks = breaks, col = "transparent", border = "transparent", prob = TRUE, 
     cex.main = 1)
hist(Mode.mode.conv$branch.dist, add = T, border = "white", col = "darkgray", breaks = breaks, prob = TRUE)
hist(Mode.morph.conv$branch.dist, add = T, border = "black", col = "transparent", breaks = breaks, prob = TRUE)
legend("topright", inset = 0, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2, bty = "n")
# dev.off()



# What proportion of taxon pairs are more than 90% convergent?
round(100 * length(which(Mode.morph.conv$C1 >= 0.9)) / nrow(Mode.morph.conv), 2)
round(100 * length(which(Mode.mode.conv$C1 >= 0.9)) / nrow(Mode.mode.conv), 2)
round(100 * length(which(Mode.constant.conv$C1 >= 0.9)) / nrow(Mode.constant.conv), 2)
round(100 * length(which(Mode.raw.conv$C1 >= 0.9)) / nrow(Mode.raw.conv), 2)
# morphology: 0.52% (696 pairs)
#    ecology: 1.74% (mode), 6.60% (constant), 1.92% (raw)

# What proportion of taxon pairs are 100% convergent?
round(100 * length(which(Mode.morph.conv$C1 == 1)) / nrow(Mode.morph.conv), 3)
round(100 * length(which(Mode.mode.conv$C1 == 1)) / nrow(Mode.mode.conv), 3)
round(100 * length(which(Mode.constant.conv$C1 == 1)) / nrow(Mode.constant.conv), 3)
round(100 * length(which(Mode.raw.conv$C1 == 1)) / nrow(Mode.raw.conv), 3)
# morphology: 0.005% (6 pairs)
#    ecology: 1.74% (mode), 6.60% (constant), 1.92% (raw)

# Among highly convergent taxa, are they convergent both anatomically and
# ecologically?
cor(Mode.morph.conv$C1[which(Mode.morph.conv$C1 >= 0.9)],
    Mode.mode.conv$C1[which(Mode.morph.conv$C1 >= 0.9)], use = "complete.obs")
# r = 0.0766

lm.C1 <-
  lm(Mode.mode.conv$C1[which(Mode.morph.conv$C1 >= 0.9)] ~ Mode.morph.conv$C1[which(Mode.morph.conv$C1 >= 0.9)])
plot(Mode.morph.conv$C1[which(Mode.morph.conv$C1 >= 0.9)], 
    Mode.mode.conv$C1[which(Mode.morph.conv$C1 >= 0.9)])
summary(lm.C1) # p = 0.0435, r2 = 0.004
abline(lm.C1, col = "red")
# Although statistically correlated, the low r and r2 mean essentially
# independent.

# Which taxa are 100% convergent?
wh <- which(Mode.mode.conv$C1 == 1)
(convs <- unique(unlist(c(Mode.mode.conv[wh , 1:2]))))    # 295 genera in mode
# for examples, only include those coded at sp/gen-level
data <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", 
                 header = TRUE, stringsAsFactors = FALSE)
gen.sps <- data$Genus[which(data$EcologyScale == "Genus" | data$EcologyScale ==
                              "Subgenus" | data$EcologyScale == "Species")]
examples <- Mode.mode.conv[wh, 1:2]
examples[,1] <- examples[,1] %in% gen.sps
examples[,2] <- examples[,2] %in% gen.sps
good.ones <- which(apply(examples, 1, sum) == 2L)
(examples <- Mode.mode.conv[names(good.ones), ])

wh <- which(Mode.constant.conv$C1 == 1)
unique(unlist(c(Mode.constant.conv[wh , 1:2]))) # 349 genera in constant

wh <- which(Mode.raw.conv$C1 == 1)
unique(unlist(c(Mode.raw.conv[wh , 1:2])))      # 209 genera in raw treatment

# Confirm some examples (showing role of some missing states):
gens <- c("Cincinnaticrinus", "Doliocrinus", "Ramseyocrinus", "Serendipocrinus")
data[match(gens, data$Genus), c(10, 21:60)]

# 42 major clusters of functionally identical (or nearly so, depending on
# missing states) but phylogenetically unrelated convergences (mode treatment).
# Note that these results have not been updated since original version, using a
# single 'equal'-method time-scale tree. However, the matchings are essentially
# unchanged using the 50 cal3 trees used below.

# Phylum: Aristocystites = Bizarroglobus = Bockia = Cigara = Edrioaster = Kailidiscus = Kinzercystis = Vyscystis = Walcottidiscus
# Phylum: Cardiocystites = Columbicrinus = Merocrinus
# Phylum: Cheirocystis = Chirocrinus = Streptaster      *** Large phylogenetic branch distance ***
# Phylum: Paradiabolocrinus = Cotyacrinus = Rhopalocystis
# Phylum: Eknomocrinus = Cnemecrinus = Lyracystis
# Phylum: Felbabkacystis = Titanocrinus
# Phylum: Myeinocystites = Vizcainocarpus
# Phylum: Ponticulocarpus = Syringocrinus
# Subphylum: Amecystis = Belemnocystites                *** Large phylogenetic branch distance ***
# Subphylum: Asturicystis = Cardiocystella = Ctenocystis = Etoctenocystis = Jugoszovia
# Subphylum: Anulocrinus = Calceocrinus = Charactocrinus = Cyathocystis = Sprinkleocystis
# Subphylum: Cnemidactis = Stenaster
# Subphylum: Persiadiskos = Picassocrinus
# Subphylum: Apektocrinus = Cambraster = Carneyella = Edriodiscus = Felbahacystis = Petalocystites
# Subphylum: Anedriophus = Apektocrinus = Balangicystis = Echinosphaerites = Gogia = Guizhoueocrinus = Helicocystis = Helicoplacus = Hybocrinus = Hybocystis = Mandalacystis = Oklahomacystis = Sphaeronites = Tatonkacystis = Trachelocrinus = Turbanicystis = Wudingeocrinus
# Class: Abludoglyptocrinus = Anisocrinus = Haimacystis = Haptocrinus = Porocrinus = Tryssocrinus
# Class: Alisocrinus = Caleidocrinus = Dystactocrinus
# Class: Archaeocrinus = Dendrocrinus = Grenprisia
# Class: Barroubiocystis = Dibrachicystis = Vizcainoia
# Class: Cigara = Comarocystites = Lepidocystis = Llanocystis = Nolichuckia
# Class: Dendrocystites = Iowacystis
# Class: Dibrachicystis = Sanducystis = Velieuxicystis = Vizcainoia
# Class: Eopatelliocrinus = Morenacrinus
# Class: Glenocrinus = Reteocrinus
# Class: Macrostylocrinus = Praecupulocrinus
# Class: Ottawacrinus = Trichinocrinus
# Subclass: Cleiocrinus = Clidochirus = Habrotecrinus
# Subclass: Coralcrinus = Cornucrinus = Ectenocrinus = Glyptocrinus = Stromatocystites = Totiglobus
# Subclass: Eustenocrinus or Heviacrinus = Praecursoricrinus = Tunguskocrinus
# Order: Aspidocarpus = Balanocystites = Castericystis = Ctenoimbricata = Davidocinctus = Graciacystis = Gyrocystis = Lagynocystis = Nelegerocystis = Pahvanticystis = Rozanovicystis = Scalenocystis = Sucocystis = Trochocystites
# Order: Astakocrinus = Xenocrinus
# Order: Ceratocystis = Cothurnocystis
# Order: Delgadocrinus = Periglyptocrinus
# Order: Nevadaecystis = Reticulocarpos
# Order: Archaepyrgus = Chatsworthia = Coleicarpus = Hadrodiscus = Macurdablastus = Picassocrinus
# Suborder: Ateleocystites = Guichenocarpos = Mitrocystella
# Suborder: Eodimerocrinites = Gaurocrinus = Goyacrinus
# Superfamily: Dalicrinus = Rhaphanocrinus = Visocrinus
# Family: Cotylacrinna = Paradiabolocrinus
# Family: Akadocrinus = Alanisicystis = Cambroblastus = Echinoencrinites = Eustypocystis = Globoeocrinus = Lichenoides = Marjumicystis = Pareocrinus = Picassocrinus = Sinoeocrinus = Ubaghsicystis
# Family: Ambonacrinus = Diabolocrinus = Fombuenacrinus = Proexenocrinus = Pycnocrinus = Ursucrinus

# Genera that are more than 98% morphologically identical
wh <- which(Mode.morph.conv$C1 >= 0.98)
unique(unlist(c(Mode.morph.conv[wh , 1:2])))      #  31 genera in morphological data set
Mode.morph.conv[wh[order(Mode.morph.conv$C1[wh], decreasing = TRUE)], ]
# Class: Eopetalocrinus & Grammocrinus & Pentamerocrinus & Putilovocrinus *** = 1 ***
#        *** NOTE THESE ARE FALSE POSITIVES, CAUSED BY PAIRINGS WITH
#        EOPETALOCRINUS, WHICH HAS 27 MISSING STATES, WHICH THE WILLS GED
#        DISTANCE METRIC FILLS IN WITH THE PAIRED KNOWN STATE! ***
# Class: ... also cluster near Alphacrinus & Cataraquicrinus & Cefnocrinus & Cincinnaticrinus & Coralcrinus & Daedalocrinus & Delgadocrinus & Haptocrinus & Heviacrinus & Inyocrinus & Iocrinus & Maennilicrinus & Morenacrinus & Ohiocrinus & Peltacrinus & Penicillicrinus & Peniculocrinus & Pogonipocrinus & Ristnacrinus & Tunguskocrinus
# Class: Alisocrinus & Anisocrinus & Canistrocrinus & Compsocrinus & Dalicrinus & Neoarchaeocrinus & Visocrinus
# Class: Euspirocrinus & Heviacrinus & Hoplocrinus (& Picassocrinus at Subclass)
# Class: Palaeaster & Schuchertia
# Subclass: Cluster of Euspirocrinus & Heviacrinus & Hoplocrinus & Hybocrinus & Isotomocrinus & Picassocrinus & Praecursoricrinus & Tornatilicrinus
# Suborder: Ambonacrinus & Eodimerocrinites 






# What taxonomic ranks are highly convergent pairings?
rank.order <- c("Subfamily", "Family", "Superfamily", "Suborder", "Order",
                "Subclass", "Class", "Subphylum", "Phylum")
abbr.ranks <- c("Subf", "Fam", "Supf", "Subor", "Or", "Subcl", "Cl", "Subph", 
                "Phy")
min.conv <- 0.90
wh.morph <- which(Mode.morph.conv$C1 >= min.conv)
tab.morph <- table(factor(Mode.morph.conv$PairRank[wh.morph], levels = rank.order))
wh.mode <- which(Mode.mode.conv$C1 >= min.conv)
tab.mode <- table(factor(Mode.mode.conv$PairRank[wh.mode], levels = rank.order))
wh.constant <- which(Mode.constant.conv$C1 >= min.conv)
tab.constant <- table(factor(Mode.constant.conv$PairRank[wh.constant], levels = rank.order))
wh.raw <- which(Mode.raw.conv$C1 >= min.conv)
tab.raw <- table(factor(Mode.raw.conv$PairRank[wh.raw], levels = rank.order))

# Summary statistics
tab.morph
tab.mode
round(100 * tab.morph / sum(tab.morph), 1)
round(100 * tab.mode / sum(tab.mode), 1)
round(100 * cumsum(rev(tab.morph / sum(tab.morph))), 1)
round(100 * cumsum(rev(tab.mode / sum(tab.mode))), 1)

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1))
barplot(tab.morph, names.arg = abbr.ranks, cex.names = .4, main = "morphology")
barplot(tab.mode, names.arg = abbr.ranks, cex.names = .4, main = "mode")
barplot(tab.constant, names.arg = abbr.ranks, cex.names = .4, main = "constant")
barplot(tab.raw, names.arg = abbr.ranks, cex.names = .4, main = "raw")
par(op)

# pdf(file = "convergence_ranks.pdf")
par(mfrow = c(2, 1), mar = c(4, 4, 2, .1))
barplot(tab.morph, names.arg = abbr.ranks, cex.names = .8, 
        xlab = "rank of genus pairs", ylab = "Number of genus pairs",
        main = "morphological convergences")
barplot(tab.mode, names.arg = abbr.ranks, cex.names = .8, 
        xlab = "rank of genus pairs", ylab = "Number of genus pairs",
        main = "ecological convergences")
par(op)
# dev.off()


# Histogram of branch distances for convergent pairings
# pdf(file = "convergence_phylo_distances.pdf")
par(mfrow = c(1,1))
breaks <-
  pretty(c(Mode.morph.conv[wh.morph, "branch.dist"], Mode.mode.conv[wh.mode, "branch.dist"]), 10)
leg.eco <- paste0("ecology (n = ", length(wh.mode), " pairs)")
leg.morph <- paste0("morphology (n = ", length(wh.morph), " pairs)")
hist(c(Mode.morph.conv[wh.morph, "branch.dist"], Mode.mode.conv[wh.mode, "branch.dist"]),
     main = "Phylogenetic distance among highly convergent taxon pairs", 
     xlab = "branch distance", ylab = "Density", breaks = breaks, 
     col = "transparent", border = "transparent", prob = TRUE, cex.main = 1)
hist(Mode.mode.conv[wh.mode, "branch.dist"], add = T, border = "white", col = "darkgray", breaks = breaks, prob = TRUE)
hist(Mode.morph.conv[wh.morph, "branch.dist"], add = T, border = "black", col = "transparent", breaks = breaks, prob = TRUE)
legend("topright", inset = 0, c(leg.eco, leg.morph), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 0.9, pt.cex = 2, bty = "n")
# dev.off()

# Are these distributions different?
summary(Mode.morph.conv[wh.morph, "branch.dist"])       # Median = 11.96
summary(Mode.mode.conv[wh.mode, "branch.dist"])         # Median =  5.89
summary(Mode.constant.conv[wh.constant, "branch.dist"]) # Median =  5.49
summary(Mode.raw.conv[wh.raw, "branch.dist"])           # Median =  5.92

ks.test(Mode.morph.conv[wh.morph, "branch.dist"], 
        Mode.mode.conv[wh.mode, "branch.dist"])
# Very much so: D = 0.578, p-value < 2.2e-16 ***

# Remove known problematic false-positive pairings with Eopetalocrinus
wh.morph.no.Eop <- which((Mode.morph.conv$Taxon1 != "Eopetalocrinus" & 
                            Mode.morph.conv$Taxon2 != "Eopetalocrinus")
                         & Mode.morph.conv$C1 >= min.conv)
summary(Mode.morph.conv[wh.morph.no.Eop, "branch.dist"])
# Median = 12.14 (this genus has negligible effect overall on results)

# Use Mantel test to confirm there is not a bias caused by differences in the
# two distance matrices
load("morph.distances.GED.5")
load("mode.distances.GED.5")

ntrees <- length(mode.distances.GED.5)
mantel.results <- vector("list", ntrees)
(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
mantel.results <- foreach(t = 1:ntrees, .inorder = TRUE, .packages = "ade4") %dopar% { 
  out <- rep(NA, 4)
  set.seed(314)  # Set RNG seed to allow replication
  d1 <- morph.distances.GED.5[[t]]$distance_matrix
  d2 <- mode.distances.GED.5[[t]]$distance_matrix
  diag(d1) <- diag(d2) <- NA
  d1 <- as.dist(d1)
  d2 <- as.dist(d2)
  # In case of incalculable errors (which occur for trees 31-35)
  results <- try(ade4::mantel.rtest(d1, d2, nrepet = 999))
  if (!inherits(results, "try-error"))
    out <- c(results$obs, results$pvalue, max(as.vector(d1)), max(as.vector(d2)))
  return(out)
}
stopCluster(cl)
(Sys.time() - start) # 6.3 minutes on 8-core laptop
beep(3)
# save(mantel.results, file = "mantel.results")

# Summarize using mean values
mean.mantel <- matrix(apply(simplify2array(mantel.results), 1, mean, na.rm = TRUE), nrow = 1)
colnames(mean.mantel) <- c("Mantel.r", "p-val", "morph.max", "mode.max")
mean.mantel

# r = 0.5529, p-value = 0.001, therefore significantly positively correlated ***
# and they span similar ranges (0 - 8.1 for morph and 0 - 7.9 for mode)



## EXAMPLES ####################################################################

## IMPROVED ALL.EQUAL THAT ONLY PRINTS COLUMNS WHEN ANY NUMBER OF PAIRS ARE
## DISSIMILAR.
## a    = minimum of two integers (corresponding to rows in a matrix)
## data = data matrix (where rows = tip & ancestors and cols = characters)
multi.all.equal <- function(a = NULL, data = NULL) {
  nr <- length(a)
  if (nr < 2)
    stop("'a' should have more than 1 item to compare.\n")
  nc <- seq.int(ncol(data))
  row.names.data <- row.names(data[a, ])
  # Convert gaps to NAs (but allow polymorphisms)
  data1 <- matrix(data[a, ], nrow = 1)
  data1[which(data1 == "")] <- NA
  # Note drops unused rows in next line!
  data3 <- matrix(data1, nrow = nr, 
                  dimnames = list(row.names.data, nc))
  wh.diff <- 
    sapply(nc, function(nc) length(na.omit(unique(data3[, nc]))) > 1L)
  return(matrix(data3[, which(wh.diff)], nrow = nr,
                dimnames = list(row.names.data, nc[which(wh.diff)])))
}

# Load necessary objects:
load("mode.pcoa"); load("Mode.mode.conv"); load("morph.anc")
load("morph.pcoa"); load("Mode.morph.conv"); load("mode.anc")

# The order here needs to be in alphabetical order to work
tp <- matrix(c("Anedriophus", "Gogia"), nrow = 1)
# tp <- matrix(c("Amecystis", "Belemnocystites"), nrow = 1)
# tp <- matrix(c("Carneyella", "Isorophus"), nrow = 1)
# tp <- matrix(c("Cambroblastus", "Estonocystis"), nrow = 1)
# tp <- matrix(c("Edrioaster", "Kinzercystis"), nrow = 1)
# tp <- matrix(c("Cheirocystis", "Streptaster"), nrow = 1)
# tp <- matrix(c("Cnemidactis", "Stenaster"), nrow = 1)
# tp <- matrix(c("Cheirocystis", "Cheirocrinus"), nrow = 1)
# tp <- matrix(c("Caleidocrinus", "Glaucocrinus"), nrow = 1)
# tp <- matrix(c("Cnemecrinus", "Picassocrinus"), nrow = 1)
# tp <- matrix(c("Isotomocrinus", "Picassocrinus"), nrow = 1)
# tp <- matrix(c("Alisocrinus", "Anisocrinus"), nrow = 1)

# Other interesting pairings:
# Edrioasteroid Anedriophus & eocrinoid Gogia
#     100% ecologically & 4% morphologically convergent
# Rhombiferan Amecystis & solute Belemnocystites
#     100% ecologically & 47% morphologically convergent with long branch distances
# Edrioasteroids Carneyella & Isorophus
#     0% ecologically & 92% morphologically convergent
# Edrioasteroid Cambroblastus & 'diploporitan' Estonocystis
#     0% ecologically & 21% morphologically convergent
# Edrioasteroid Edrioaster & eocrinoid Kinzercystis:
#     100% ecologically & 12% morphologically convergent
# Edrioasteroid Streptaster & rhombiferan Cheirocystis
#     100% ecologically & 0% morphologically convergent with long branch distances
# Disparid crinoids in different suborders Caleidocrinus & Glaucocrinus
#     29% ecologically and 98% morphologically convergent
# Crinoids in different subclasses Cnemecrinus & Picassocrinus
#     22% ecologically and 92% morphologically convergent
# Crinoids in different subclasses Isotomocrinus & Picassocrinus
#     32% ecologically and 98% morphologically convergent

# Confirm visually on phylomorphospace & phyloecospace:
# pdf(file = "ConvPair2.pdf")
par(mar = c(5, 4, 2, 2))

# Use tree #50 (which is the most similar to the consensus tree). Make sure that
# the same tree is used for the same corresponding distance matrix, as the
# identity of ancestral nodes (e.g., anc367, anc 453) matches for the same list
# item, but may not match across list items. (The tips match across all, but the
# order of tree$tip.label may not match.)
t <- 50
tree <- mode.pcoa[[t]]$tree
# Stayton statistics
Mode.morph.conv[which(Mode.morph.conv$Taxon1 == tp[1] & Mode.morph.conv$Taxon2 == tp[2]), ]
Mode.mode.conv[which(Mode.mode.conv$Taxon1 == tp[1] & Mode.mode.conv$Taxon2 == tp[2]), ]

# Plotting settings (needing to match tree tip names with pcoa rownames b/c
# different order for tips)
wh.tip <- which(rownames(morph.pcoa[[t]]$vectors.cor) == tp[1] | 
                  rownames(morph.pcoa[[t]]$vectors.cor) == tp[2])
root <- Ntip(tree) + 1
tip.seq <- 1:Ntip(tree)
node.seq <- root:(Ntip(tree) + Nnode(tree))
con <- list(col.edge = setNames(rep("lightgray", nrow(tree$edge)), 
                                as.character(tree$edge[, 2])))
# Find ancestral edges for pair of taxa until united in MRCA (function defined
# above). The MRCA should be the first [1,1] and then be repeated for second
# lineage.
(branching.history <- ancestral.lineages(tree, t1 = tp[1], t2 = tp[2]))
# Identify edges since MRCA
MRCA <- branching.history[1]
break.point <- which(branching.history[, 1] == MRCA)
hist.first <- matrix(branching.history[1:diff(break.point), ], ncol = 2)
hist.second <- matrix(branching.history[break.point[2]:nrow(branching.history), ], ncol = 2)
# Switch so that labels and branch histories are the same color
if (hist.first[nrow(hist.first), 2] == wh.tip[1]) {
  hist1 <- hist.first
  hist2 <- hist.second
} else {
  hist1 <- hist.second
  hist2 <- hist.first
}
cols <- c("black", "red")[order(wh.tip)]
# Plot
par(mfrow = c(2, 1), mar = c(4, 4, 1, .5))
phytools::phylomorphospace(tree = tree, 
                           X = morph.pcoa[[t]]$vectors.cor[tip.seq, 1:2], 
                           A = morph.pcoa[[t]]$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = "PCoA 1", 
                           ylab = "PCoA 2", pch = NA)
mtext("Phylomorphospace", 3)
text(-3, 4, tp[1], pos = 4, col = cols[1])
text(5, 4, tp[2], pos = 2, col = cols[2])
points(x = morph.pcoa[[t]]$vectors.cor[branching.history[1], 1], 
       y = morph.pcoa[[t]]$vectors.cor[branching.history[1], 2], col = "indianred1", pch = 15, 
       cex = 1.5) # MRCA
points(x = morph.pcoa[[t]]$vectors.cor[wh.tip, 1], 
       y = morph.pcoa[[t]]$vectors.cor[wh.tip, 2], 
       col = cols, pch = c(16, 17), cex = 1.25)
for(r in 1:nrow(hist1)) {
  anc.pt <- hist1[r, 1]
  desc.pt <- hist1[r, 2]
  # Switch needed for tips b/c order of tree$tip.label doesn't match order in pcoa
  if (desc.pt <= Ntip(tree)) desc.pt <- wh.tip[1]
  lines(x = morph.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = morph.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[1], lwd = 1)
}
for(r in 1:nrow(hist2)) {
  anc.pt <- hist2[r, 1]
  desc.pt <- hist2[r, 2]
  if (desc.pt <= Ntip(tree)) desc.pt <- wh.tip[2]
  lines(x = morph.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = morph.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[2], lwd = 1)
}
# Same for phyloecospace
phytools::phylomorphospace(tree = tree, 
                           X = mode.pcoa[[t]]$vectors.cor[tip.seq, 1:2], 
                           A = mode.pcoa[[t]]$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = "PCoA 1", 
                           ylab = "PCoA 2", pch = NA)
mtext("Phyloecospace", 3)
points(x = mode.pcoa[[t]]$vectors.cor[branching.history[1], 1], 
       y = mode.pcoa[[t]]$vectors.cor[branching.history[1], 2], col = "indianred1", pch = 15, 
       cex = 1.5) # MRCA
points(x = mode.pcoa[[t]]$vectors.cor[wh.tip, 1], 
       y = mode.pcoa[[t]]$vectors.cor[wh.tip, 2], 
       col = cols, pch = c(16, 17), cex = 1.25)
for(r in 1:nrow(hist1)) {
  anc.pt <- hist1[r, 1]
  desc.pt <- hist1[r, 2]
  if (desc.pt <= Ntip(tree)) desc.pt <- wh.tip[1]
  lines(x = mode.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = mode.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[1], lwd = 1)
}
for(r in 1:nrow(hist2)) {
  anc.pt <- hist2[r, 1]
  desc.pt <- hist2[r, 2]
  if (desc.pt <= Ntip(tree)) desc.pt <- wh.tip[2]
  lines(x = mode.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = mode.pcoa[[t]]$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[2], lwd = 1)
}
par(op)
# dev.off()


# Show tip and MRCA states for comparison of how changed
(morph.MRCA <- 
    multi.all.equal(a = c(wh.tip, MRCA), data = morph.anc[[t]]$matrix_1$matrix))
(eco.MRCA <- 
    multi.all.equal(a = c(wh.tip, MRCA), data = mode.anc[[t]]$matrix_1$matrix))

# Only show characters where converged
morph.MRCA[, which(morph.MRCA[1, ] == morph.MRCA[2, ])]
eco.MRCA[, which(eco.MRCA[1, ] == eco.MRCA[2, ])]

