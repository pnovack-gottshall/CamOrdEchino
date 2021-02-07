## CONVERGENCE ANALYSES (STAYTON 2015) #########################################

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' package. 'convevol' assumes Brownian motion ancestral
# reconstructions, which are not appropriate for the discrete characters we are
# using. Instead, we will use Claddis::AncStateEstMatrix, which wraps around
# phytools::rerootingMethod, to reconstruct ancestral states using the
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
setwd("~/Manuscripts/CamOrdEchinos/Data files/NA Reformatted")

# Load packages
library(beepr)      # v. 1.3
library(doParallel) # v. 1.0.15
library(Claddis)    # v. 0.4.1
library(ade4)       # v. 1.7-15
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub\n")




## IMPORT FILES ################################################################

# Import output from 3-DisparityDistances.R
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")

# Import time trees saved from 1-MakeTimeTrees.R
load("~/Manuscripts/CamOrdEchinos/equal.tree")
tree <- equal.tree



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
    CurrentMRCA <- FindAncestor(descs = taxon.pairs[i, ], tree = tree)
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

# Get names for all pairs of tips (surprisingly slow)
taxon.pairs <- matrix(nrow = 0, ncol = 2)
for (i in 1:(Ntip(tree) - 1)) {
  for (j in (i + 1):Ntip(tree)) {
    taxon.pairs <- rbind(taxon.pairs, c(tree$tip.label[i], tree$tip.label[j]))
  }
}
beep()
# save(taxon.pairs, file = "taxon.pairs")
# load("taxon.pairs")



# Calculate C1 & C2 convergence statistics for these pairs of taxa
morph.conv <-
  StaytonConvergence(dist.matrix = morph.distances.GED.5$DistanceMatrix,
                     tree = tree, taxon.pairs = taxon.pairs)
mode.conv <-
  StaytonConvergence(dist.matrix = mode.distances.GED.5$DistanceMatrix, 
                     tree = tree, taxon.pairs = taxon.pairs)
constant.conv <- 
  StaytonConvergence(dist.matrix = constant.distances.GED.5$DistanceMatrix,
                     tree = tree, taxon.pairs = taxon.pairs)
raw.conv <-
  StaytonConvergence(dist.matrix = raw.distances.GED.5$DistanceMatrix, 
                     tree = tree, taxon.pairs = taxon.pairs)
beep(3)


# Calculate branch distances (used for C3) for these pairs of taxa (in parallel)
sq <- 1:nrow(taxon.pairs)
morph.BD <- mode.BD <- constant.BD <- raw.BD <- rep(NA, nrow(taxon.pairs))
CPUs <- detectCores(logical = TRUE)
cl <- makeCluster(CPUs)
registerDoParallel(cl)

# Morphology
morph.BD <- foreach(i = sq, .combine = c) %dopar% {
              BH <- ancestral.lineages(tree = tree, t1 = taxon.pairs[i, 1],
                                       t2 = taxon.pairs[i, 2])
              morph.BD.pair <- branch.distance(BH, morph.distances.GED.5$DistanceMatrix)
}

# Mode
mode.BD <- foreach(i = sq, .combine = c) %dopar% {
              BH <- ancestral.lineages(tree = tree, t1 = taxon.pairs[i, 1],
                                       t2 = taxon.pairs[i, 2])
              mode.BD.pair <- branch.distance(BH, mode.distances.GED.5$DistanceMatrix)
}

# Constant
constant.BD <- foreach(i = sq, .combine = c) %dopar% {
                BH <- ancestral.lineages(tree = tree, t1 = taxon.pairs[i, 1],
                                         t2 = taxon.pairs[i, 2])
                constant.BD.pair <- branch.distance(BH, constant.distances.GED.5$DistanceMatrix)
}

# Raw
raw.BD <- foreach(i = sq, .combine = c) %dopar% {
            BH <- ancestral.lineages(tree = tree, t1 = taxon.pairs[i, 1],
                                     t2 = taxon.pairs[i, 2])
            raw.BD.pair <- branch.distance(BH, raw.distances.GED.5$DistanceMatrix)
}
stopCluster(cl)
beep(3)

# Calculate C3 ( = C2 / branch distance) and append as new column
morph.conv$C3 <- morph.conv$C2 / morph.BD
mode.conv$C3 <- mode.conv$C2 / mode.BD
constant.conv$C3 <- constant.conv$C2 / constant.BD
raw.conv$C3 <- raw.conv$C2 / raw.BD

# Append branch-length distances as new column
morph.conv$branch.dist <- morph.BD
mode.conv$branch.dist <- mode.BD
constant.conv$branch.dist <- constant.BD
raw.conv$branch.dist <- raw.BD

# Save output
# save(mode.conv, file = "mode.conv")
# save(constant.conv, file = "constant.conv")
# save(raw.conv, file = "raw.conv")
# save(morph.conv, file = "morph.conv")




## WHAT TAXONOMIC RANK ARE THE PAIRINGS? #######################################
# In other words, what is the taxonomic distance among convergent pairs?

# Import ecological data set, which has all taxonomic assignments:
data <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", 
                 header = TRUE, stringsAsFactors = FALSE)
data[1:10, 1:5]
# Modify names for two subgenera treated as genera
data$Genus[which(data$Genus == "Anatifopsis")] <-
  c("Anatiferocystis", "Guichenocarpos")

# Assign rank for morphological data set (which is same for all data sets)
rank <- rep(NA, nrow(morph.conv))
for (r in 1:length(rank)) {
  pair <- morph.conv[r, 1:2]
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
morph.conv$PairRank <- rank
mode.conv$PairRank <- rank
constant.conv$PairRank <- rank
raw.conv$PairRank <- rank

# save(mode.conv, file = "mode.conv")
# save(constant.conv, file = "constant.conv")
# save(raw.conv, file = "raw.conv")
# save(morph.conv, file = "morph.conv")
load("mode.conv")
load("constant.conv")
load("raw.conv")
load("morph.conv")


## A curious outcome explained:

# Stayton (2015) claims C3 should have a maximum value of 1 because it is the
# proportion of C2 convergence accounted for in branch evolution. But 82 taxon
# pairs return C3 values greater than 1. A few are easily explained as rounding
# errors (e.g., C3 = 1 + 3.5e-15). Most of the remaining "impossible" values
# here occur as a result of a few taxa, especially pairings with Apektocrinus,
# Haplosphaeronis, and Cryptocrinites. Infinite and "> 1" C3 values occur
# when the tips, MRCA, and intervening lineages all have zero distances. In some
# of these cases, they may not be exactly identical (and may have different
# coordinates in PCoA space) because of pairwise distance handling for missing
# character states. Overall, there 71 such pairings in the mode treatment, 9 in
# the morphology, 825 in constant, and 70 in the raw, accounting for < 0.1% of
# the cases, and these can be ignored. The higher proportion in the constant is
# additional confirmation these are caused by missing character states.
rounders <- which(mode.conv$C3 > 1 & mode.conv$C3 < 1.001)
(mode.conv$C3[rounders] - 1)

oddballs <- which(mode.conv$C3 > 1.001 & is.finite(mode.conv$C3))
length(oddballs)
sort(table(c(mode.conv$Taxon1[oddballs], mode.conv$Taxon2[oddballs])))
mode.conv[oddballs, ]

odder.oddballs <- which(is.infinite(mode.conv$C3))
sort(table(c(mode.conv$Taxon1[odder.oddballs], mode.conv$Taxon2[odder.oddballs])))
mode.conv[odder.oddballs, ]





## ANALYSE PATTERNS ############################################################

# Plot C1 histograms:
par(mfrow = c(2,2))
breaks <- seq(from = 0, to = 1, by = 0.05)
h.morph1 <- hist(morph.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "morphology")
h.mode1 <- hist(mode.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (mode treatment)")
h.constant1 <- hist(constant.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (constant treatment)")
h.raw1 <- hist(raw.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (raw treatment)")

# Plot C2 histograms:
breaks2 <- seq(from = 0, to = 5, by = 0.25)
h.morph2 <- hist(morph.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "morphology")
h.mode2 <- hist(mode.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (mode treatment)")
h.constant2 <- hist(constant.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (constant treatment)")
h.raw2 <- hist(raw.conv$C2, breaks = breaks2, prob = TRUE, xlab = "C2", main = "ecology (raw treatment)")

# Plot C3 histograms (truncated for errors caused by missing data when
# calculating pairwise distance measures):
breaks3 <- seq(from = 0, to = 1, by = 0.05)
# Remove errors caused by rounding and missing data problems (see below)
h.morph3 <- hist(morph.conv$C3[morph.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "morphology")
h.mode3 <- hist(mode.conv$C3[mode.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (mode treatment)")
h.constant3 <- hist(constant.conv$C3[constant.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (constant treatment)")
h.raw3 <- hist(raw.conv$C3[raw.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (raw treatment)")

# Plot branch-length-distance histograms:
breaks4 <- seq(from = 0, to = 30, by = 1.5)
h.morph4 <- hist(morph.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "morphology")
h.mode4 <- hist(mode.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (mode treatment)")
h.constant4 <- hist(constant.conv$branch.dist, breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (constant treatment)")
h.raw4 <- hist(na.omit(raw.conv$branch.dist), breaks = breaks4, prob = TRUE, xlab = "branch-length distance", main = "ecology (raw treatment)")

# Are the convergence statistics correlated between ecological and morphological
# data sets?
cor(morph.conv$C1, mode.conv$C1, use = "complete.obs") # r = 0.0588
cor(morph.conv$C2, mode.conv$C2)                       # r = 0.0582
wh.inf <- which(is.infinite(mode.conv$C3))
cor(morph.conv$C3[-wh.inf], mode.conv$C3[-wh.inf], use = "complete.obs") # r = 0.029
cor(morph.conv$branch.dist, mode.conv$branch.dist)     # r = 0.526
plot(morph.conv$branch.dist, mode.conv$branch.dist)    # but lots of noise



# pdf(file = "convergence_statistics.pdf")
par(mfrow = c(2,2))
breaks <- seq(from = 0, to = 1, by = 0.05)
h.morph1 <- hist(morph.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "morphology")
h.mode1 <- hist(mode.conv$C1, breaks = breaks, prob = TRUE, xlab = "C1", main = "ecology (mode treatment)")
breaks3 <- seq(from = 0, to = 1, by = 0.05)
h.morph3 <- hist(morph.conv$C3[morph.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "morphology")
h.mode3 <- hist(mode.conv$C3[mode.conv$C3 <= 1.00001], breaks = breaks3, prob = TRUE, xlab = "C3", main = "ecology (mode treatment)")
# dev.off()

# Median pairwise convergence:
round(apply(morph.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(mode.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(constant.conv[, 5:8], 2, median, na.rm = TRUE), 2)
round(apply(raw.conv[, 5:8], 2, median, na.rm = TRUE), 2)
# Morphology: median C1 = 0.01, C3 = 0.03, branch.distance = 13.7
# Eco.(mode): median C1 = 0.10, C3 = 0.03, branch.distance =  7.06

# Are the distributions different?
ks.test(morph.conv$C1, mode.conv$C1)     # D = 0.141, p < 2.2e-16
ks.test(morph.conv$C1, constant.conv$C1) # D = 0.240, p < 2.2e-16
ks.test(morph.conv$C1, raw.conv$C1)      # D = 0.167, p < 2.2e-16

ks.test(morph.conv$C3, mode.conv$C3)     # D = 0.141, p < 2.2e-16
ks.test(morph.conv$C3, constant.conv$C3) # D = 0.306, p < 2.2e-16
ks.test(morph.conv$C3, raw.conv$C3)      # D = 0.202, p < 2.2e-16

ks.test(morph.conv$branch.dist, mode.conv$branch.dist)     # D = 0.439, p < 2.2e-16
ks.test(morph.conv$branch.dist, constant.conv$branch.dist) # D = 0.475, p < 2.2e-16
ks.test(morph.conv$branch.dist, raw.conv$branch.dist)      # D = 0.679, p < 2.2e-16


# Histogram of branch distances for convergent pairings
# pdf(file = "branch_distances.pdf")
par(mfrow = c(1,1))
breaks <-
  pretty(c(morph.conv$branch.dist, mode.conv$branch.dist), 10)
hist(c(morph.conv$branch.dist, mode.conv$branch.dist),
     main = "Phylogenetic distances", xlab = "branch distance among taxon pairs", ylab = "Density", breaks = breaks, 
     col = "transparent", border = "transparent", prob = TRUE, cex.main = 1)
hist(mode.conv$branch.dist, add = T, border = "white", col = "darkgray", breaks = breaks, prob = TRUE)
hist(morph.conv$branch.dist, add = T, border = "black", col = "transparent", breaks = breaks, prob = TRUE)
legend("topright", inset = 0, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2, bty = "n")
# dev.off()



# What proportion of taxon pairs are more than 90% convergent?
round(100 * length(which(morph.conv$C1 >= 0.9)) / nrow(morph.conv), 2)
round(100 * length(which(mode.conv$C1 >= 0.9)) / nrow(mode.conv), 2)
round(100 * length(which(constant.conv$C1 >= 0.9)) / nrow(constant.conv), 2)
round(100 * length(which(raw.conv$C1 >= 0.9)) / nrow(raw.conv), 2)
# morphology: 0.79% 
#    ecology: 1.33% (mode), 5.96% (constant), 0.60% (raw)

# What proportion of taxon pairs are 100% convergent?
round(100 * length(which(morph.conv$C1 == 1)) / nrow(morph.conv), 3)
round(100 * length(which(mode.conv$C1 == 1)) / nrow(mode.conv), 3)
round(100 * length(which(constant.conv$C1 == 1)) / nrow(constant.conv), 3)
round(100 * length(which(raw.conv$C1 == 1)) / nrow(raw.conv), 3)
# morphology: 0.004% (3 pairs)
#    ecology: 1.33% (mode), 5.96% (constant), 0.60% (raw)

# Among highly convergent taxa, are they convergent both anatomically and
# ecologically?
cor(morph.conv$C1[which(morph.conv$C1 >= 0.9)],
    mode.conv$C1[which(morph.conv$C1 >= 0.9)], use = "complete.obs") # r = 0.1309

lm.C1 <-
  lm(mode.conv$C1[which(morph.conv$C1 >= 0.9)] ~ morph.conv$C1[which(morph.conv$C1 >= 0.9)])
plot(morph.conv$C1[which(morph.conv$C1 >= 0.9)], 
    mode.conv$C1[which(morph.conv$C1 >= 0.9)])
summary(lm.C1) # p = 0.00274, r2 = 0.017
abline(lm.C1, col = "red")
# Although statistically correlated, the low r and r2 mean essentially
# independent.

# Which taxa are 100% convergent?
wh <- which(mode.conv$C1 == 1)
(convs <- unique(unlist(c(mode.conv[wh , 1:2]))))    # 278 genera in mode
# for examples, only include those coded at sp/gen-level
gen.sps <- data$Genus[which(data$EcologyScale == "Genus" | data$EcologyScale ==
                              "Subgenus" | data$EcologyScale == "Species")]
examples <- mode.conv[wh, 1:2]
examples[,1] <- examples[,1] %in% gen.sps
examples[,2] <- examples[,2] %in% gen.sps
good.ones <- which(apply(examples, 1, sum) == 2L)
(examples <- mode.conv[names(good.ones), ])

wh <- which(constant.conv$C1 == 1)
unique(unlist(c(constant.conv[wh , 1:2]))) # 342 genera in constant

wh <- which(raw.conv$C1 == 1)
unique(unlist(c(raw.conv[wh , 1:2])))      # 167 genera in raw treatment

gens <- c("Cincinnaticrinus", "Doliocrinus", "Ramseyocrinus", "Serendipocrinus")
data[match(gens, data$Genus), c(10, 21:60)]

# 42 major clusters of functionally identical (or nearly so, depending on
# missing states) but phylogenetically unrelated convergences (mode treatment):
# Phylum: Aristocystites = Bizarroglobus = Bockia = Cigara = Edrioaster = Kailidiscus = Kinzercystis = Vyscystis = Walcottidiscus
# Phylum: Cardiocystites = Columbicrinus = Merocrinus
# Phylum: Cheirocystis = Chirocrinus = Streptaster      *** Large phylogenetic branch distance ***
# Phylum: Paradiabolocrinus = Rhopalocystis
# Phylum: Eknomocrinus = Lyracystis
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
# Suborder: Anatiferocystis = Ateleocystites = Guichenocarpos = Mitrocystella
# Suborder: Eodimerocrinites = Gaurocrinus = Goyacrinus
# Superfamily: Dalicrinus = Rhaphanocrinus = Visocrinus
# Family: Cotylacrinna = Paradiabolocrinus
# Family: Akadocrinus = Alanisicystis = Cambroblastus = Echinoencrinites = Eustypocystis = Globoeocrinus = Lichenoides = Marjumicystis = Pareocrinus = Picassocrinus = Sinoeocrinus = Ubaghsicystis
# Family: Ambonacrinus = Diabolocrinus = Fombuenacrinus = Proexenocrinus = Pycnocrinus = Ursucrinus

# Genera that are more than 98% morphologically identical
wh <- which(morph.conv$C1 >= 0.98)
unique(unlist(c(morph.conv[wh , 1:2])))      #  45 genera in morphological data set
morph.conv[wh[order(morph.conv$C1[wh], decreasing = TRUE)], ]
# Class: Eopetalocrinus & Grammocrinus & Pentamerocrinus & Putilovocrinus *** = 1 ***
#        *** NOTE THESE ARE FALSE POSITIVES, CAUSED BY PAIRINGS WITH
#        EOPETALOCRINUS, WHICH HAS 27 MISSING STATES, WHICH THE WILLS GED
#        DISTANCE METRIC FILLS IN WITH THE PAIRED KNOWN STATE! ***
# Class: ... also cluster near Alphacrinus & Cataraquicrinus & Cefnocrinus & Cincinnaticrinus & Coralcrinus & Daedalocrinus & Delgadocrinus & Haptocrinus & Heviacrinus & Inyocrinus & Iocrinus & Maennilicrinus & Morenacrinus & Ohiocrinus & Peltacrinus & Penicillicrinus & Peniculocrinus & Pogonipocrinus & Ristnacrinus & Tunguskocrinus
# Class: Alisocrinus & Anisocrinus & Canistrocrinus & Compsocrinus & Dalicrinus & Neoarchaeocrinus & Visocrinus
# Class: Cnemecrinus & Isotomocrinus & Picassocrinus & Tornatilicrinus
# Class: Palaeaster & Schuchertia
# Subclass:  Caleidocrinus & Glaucocrinus
# Subclass: Cluster of Euspirocrinus & Heviacrinus & Hoplocrinus & Hybocrinus & Isotomocrinus & Picassocrinus & Praecursoricrinus & Tornatilicrinus
# Suborder: Ambonacrinus & Eodimerocrinites 






# What taxonomic ranks are highly convergent pairings?
rank.order <- c("Subfamily", "Family", "Superfamily", "Suborder", "Order",
                "Subclass", "Class", "Subphylum", "Phylum")
abbr.ranks <- c("Subf", "Fam", "Supf", "Subor", "Or", "Subcl", "Cl", "Subph", 
                "Phy")
min.conv <- 0.90
wh.morph <- which(morph.conv$C1 >= min.conv)
tab.morph <- table(factor(morph.conv$PairRank[wh.morph], levels = rank.order))
wh.mode <- which(mode.conv$C1 >= min.conv)
tab.mode <- table(factor(mode.conv$PairRank[wh.mode], levels = rank.order))
wh.constant <- which(constant.conv$C1 >= min.conv)
tab.constant <- table(factor(constant.conv$PairRank[wh.constant], levels = rank.order))
wh.raw <- which(raw.conv$C1 >= min.conv)
tab.raw <- table(factor(raw.conv$PairRank[wh.raw], levels = rank.order))

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
  pretty(c(morph.conv[wh.morph, "branch.dist"], mode.conv[wh.mode, "branch.dist"]), 10)
leg.eco <- paste0("ecology (n = ", length(wh.mode), " pairs)")
leg.morph <- paste0("morphology (n = ", length(wh.morph), " pairs)")
hist(c(morph.conv[wh.morph, "branch.dist"], mode.conv[wh.mode, "branch.dist"]),
     main = "Phylogenetic distance among highly convergent taxon pairs", 
     xlab = "branch distance", ylab = "Density", breaks = breaks, 
     col = "transparent", border = "transparent", prob = TRUE, cex.main = 1)
hist(mode.conv[wh.mode, "branch.dist"], add = T, border = "white", col = "darkgray", breaks = breaks, prob = TRUE)
hist(morph.conv[wh.morph, "branch.dist"], add = T, border = "black", col = "transparent", breaks = breaks, prob = TRUE)
legend("topright", inset = 0, c(leg.eco, leg.morph), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 0.9, pt.cex = 2, bty = "n")
# dev.off()

# Are these distributions different?
summary(morph.conv[wh.morph, "branch.dist"]) # Median = 7.3
summary(mode.conv[wh.mode, "branch.dist"])   # Median = 3.0
summary(constant.conv[wh.constant, "branch.dist"]) # Median = 3.5
summary(raw.conv[wh.raw, "branch.dist"])   # Median = 2.0

ks.test(morph.conv[wh.morph, "branch.dist"], mode.conv[wh.mode, "branch.dist"])
# Very much so: D = 0.425, p-value < 2.2e-16 ***

# Remove known problematic false-positive pairings with Eopetalocrinus
wh.morph.no.Eop <- which((morph.conv$Taxon1 != "Eopetalocrinus" & 
                            morph.conv$Taxon2 != "Eopetalocrinus")
                         & morph.conv$C1 >= min.con)
summary(morph.conv[wh.morph.no.Eop, "branch.dist"])
# Median = 7.4 (this genus has negligible effect)

# Use Mantel test to confirm there is not a bias caused by differences in the
# two distance matrices
set.seed(314)  # Set RNG seed to allow replication
d1 <- morph.distances.GED.5$DistanceMatrix
d2 <- mode.distances.GED.5$DistanceMatrix
diag(d1) <- diag(d2) <- NA
d1 <- as.dist(d1)
d2 <- as.dist(d2)
ade4::mantel.rtest(d1, d2, nrepet = 999)
# r = 0.5595, p-value = 0.001, therefore significantly positively correlated
summary(as.vector(d1))
summary(as.vector(d2))
# and they span similar ranges (0 - 7.9 for both)




# First pass to remove those that are all NAs
d1.na <- apply(d1, 1, function(x) all(is.na(x)))
d2.na <- apply(d2, 1, function(x) all(is.na(x)))
d1 <- d1[!d1.na, !d1.na]
d2 <- d2[!d2.na, !d2.na]
diag(d1) <- diag(d2) <- 0
# Second pass to remove those (fewer) that still have a few NAs
d1.na <- apply(d1, 1, function(x) any(is.na(x)))
d2.na <- apply(d2, 1, function(x) any(is.na(x)))
d1 <- d1[!d1.na,!d1.na]
d2 <- d2[!d2.na,!d2.na]
}
d1 <- as.dist(d1)
d2 <- as.dist(d2)
if (remove.zeros) {
  # Replace zeros with small values (half the smallest positive value)
  d1 <- replace(d1, d1 == 0, min(d1[d1 > 0]) / 2)
  d2 <- replace(d2, d2 == 0, min(d2[d2 > 0]) / 2)
  if (length(d1) != length(d2))
    stop("inconsisent lengths when NAs removed.\n")
}
dist.dists[r, c] <- ade4::mantel.rtest(d1, d2, nrepet = 99)$obs





## EXAMPLES ####################################################################
# The order here needs to be in alphabetical order to work
tp <- matrix(c("Anedriophus", "Gogia"), nrow = 1)
# tp <- matrix(c("Amecystis", "Belemnocystites"), nrow = 1)
# tp <- matrix(c("Carneyella", "Isorophus"), nrow = 1)
# tp <- matrix(c("Edrioaster", "Kinzercystis"), nrow = 1)
# tp <- matrix(c("Cheirocystis", "Streptaster"), nrow = 1)
# tp <- matrix(c("Caleidocrinus", "Glaucocrinus"), nrow = 1)
# tp <- matrix(c("Cnemecrinus", "Picassocrinus"), nrow = 1)
tp <- matrix(c("Isotomocrinus", "Picassocrinus"), nrow = 1)

# Other interesting pairings:
# Edrioasteroid Anedriophus & eocrinoid Gogia
#     100% ecologically & 4% morphologically convergent
# Rhombiferan Amecystis & solute Belemnocystites
#     100% ecologically & 45% morphologically convergent with long branch distances
# Edrioasteroids Carneyella & Isorophus
#     0% ecologically & 94% morphologically convergent
# Edrioasteroid Edrioaster & eocrinoid Kinzercystis:
#     100% ecologically & 22% morphologically convergent
# Edrioasteroid Streptaster & rhombiferan Cheirocystis
#     100% ecologically & 0% morphologically convergent with long branch distances
# Disparid crinoids in different suborders Caleidocrinus & Glaucocrinus
#     55% ecologically and 98% morphologically convergent
# Crinoids in different subclasses Cnemecrinus & Picassocrinus
#     38% ecologically and 99% morphologically convergent
# Crinoids in different subclasses Isotomocrinus & Picassocrinus
#     46% ecologically and 98% morphologically convergent

# Confirm visually on phylomorphospace & phyloecospace:
# pdf(file = "ConvPair2.pdf")
par(mar = c(5, 4, 2, 2))
load("~/Manuscripts/CamOrdEchinos/equal.tree"); tree <- equal.tree
load("mode.pcoa"); load("mode.conv")
load("morph.pcoa"); load("morph.conv")
mode.conv[which(mode.conv$Taxon1 == tp[1] & mode.conv$Taxon2 == tp[2]), ]
morph.conv[which(morph.conv$Taxon1 == tp[1] & morph.conv$Taxon2 == tp[2]), ]
wh.tip <- which(tree$tip.label == tp[1] | tree$tip.label == tp[2])
root <- Ntip(morph.pcoa$Tree) + 1
tip.seq <- 1:Ntip(morph.pcoa$Tree)
node.seq <- root:(Ntip(morph.pcoa$Tree) + Nnode(morph.pcoa$Tree))
con <- list(col.edge = setNames(rep("lightgray", nrow(morph.pcoa$Tree$edge)), 
                                as.character(morph.pcoa$Tree$edge[, 2])))
# Find ancestral edges for pair of taxa until united in MRCA:
branching.history <- ancestral.lineages(tree, t1 = tp[1], t2 = tp[2])
# To confirm works as intended:
branching.history
MRCA <- branching.history[1]
break.point <- which(branching.history[, 1] == MRCA)
hist.first <- matrix(branching.history[1:diff(break.point), ], ncol = 2)
hist.second <- matrix(branching.history[break.point[2]:nrow(branching.history), ], ncol = 2)
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
phytools::phylomorphospace(tree = morph.pcoa$Tree, 
                           X = morph.pcoa$vectors.cor[tip.seq, 1:2], 
                           A = morph.pcoa$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = "PCoA 1", 
                           ylab = "PCoA 2", pch = NA)
mtext("Phylomorphospace", 3)
text(-32, 13, tp[1], pos = 4, col = cols[1])
text(20, 13, tp[2], pos = 2, col = cols[2])
points(x = morph.pcoa$vectors.cor[branching.history[1], 1], 
       y = morph.pcoa$vectors.cor[branching.history[1], 2], col = "indianred1", pch = 15, 
       cex = 1.5) # MRCA
points(x = morph.pcoa$vectors.cor[wh.tip, 1], y = morph.pcoa$vectors.cor[wh.tip, 2], 
       col = cols, pch = c(16, 17), cex = 1.25)
for(r in 1:nrow(hist1)) {
  anc.pt <- hist1[r, 1]
  desc.pt <- hist1[r, 2]
  lines(x = morph.pcoa$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = morph.pcoa$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[1], lwd = 1)
}
for(r in 1:nrow(hist2)) {
  anc.pt <- hist2[r, 1]
  desc.pt <- hist2[r, 2]
  lines(x = morph.pcoa$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = morph.pcoa$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[2], lwd = 1)
}
# Same for phyloecospace
phytools::phylomorphospace(tree = mode.pcoa$Tree, 
                           X = mode.pcoa$vectors.cor[tip.seq, 1:2], 
                           A = mode.pcoa$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = "PCoA 1", 
                           ylab = "PCoA 2", pch = NA)
mtext("Phyloecospace", 3)
points(x = mode.pcoa$vectors.cor[branching.history[1], 1], 
       y = mode.pcoa$vectors.cor[branching.history[1], 2], col = "indianred1", pch = 15, 
       cex = 1.5) # MRCA
points(x = mode.pcoa$vectors.cor[wh.tip, 1], y = mode.pcoa$vectors.cor[wh.tip, 2], 
       col = cols, pch = c(16, 17), cex = 1.25)
for(r in 1:nrow(hist1)) {
  anc.pt <- hist1[r, 1]
  desc.pt <- hist1[r, 2]
  lines(x = mode.pcoa$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = mode.pcoa$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[1], lwd = 1)
}
for(r in 1:nrow(hist2)) {
  anc.pt <- hist2[r, 1]
  desc.pt <- hist2[r, 2]
  lines(x = mode.pcoa$vectors.cor[c(anc.pt, desc.pt), 1], 
        y = mode.pcoa$vectors.cor[c(anc.pt, desc.pt), 2], col = cols[2], lwd = 1)
}
par(op)
# dev.off()


