## PALEOTREE ANALYSES: PREPARE THE TIME-CALIBRATED TREES #######################

## Prior to running, run code in UpdateAges&DivCurve.R in Database/Maintenance &
## update R scripts folder to download and update stratigraphic ages for the
## genera from the PBDB. Once downloaded, manually change first Anatifopsis to
## its subgenus Anatiferocystis and the second one to the subgenus
## Guichenocarpus so tip names match those used in the NEXUS tree files.

## PREPARATIONS ################################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")

library(paleotree)  # v. 3.3.25
library(phytools)   # v. 0.7-47
library(phangorn)   # v. 2.5.5
library(beepr)      # v. 1.3

## Import cladogram in Nexus format
tree <- read.nexus("echinodermtree.nex")
if (class(tree) != "phylo") stop("Nexus file is not class 'phylo")
# pdf(file = "EchinoPhylogeny_supertree.pdf", height = 20)
# plot(tree, cex = 0.4); axisPhylo()
# dev.off()


## Process stratigraphic ranges

# Import strat ranges (previously downloaded fom PBDB)
ranges <- read.csv("~/Manuscripts/CamOrdEchinos/Data files/GenusStratRanges.csv", header = TRUE)
rownames(ranges) <- ranges$Genus
ranges <- ranges[, -1]
if (!identical(tree$tip.label, rownames(ranges)))
  stop("tip labels in cladogram need to match the names in the strat ranges")

# bin_timePaleoPhy() also requires the bins used (here, from PBDB)
strat_names <-
  read.csv("https://www.paleobiodb.org/data1.2/intervals/list.csv?all_records&vocab=pbdb")
head(strat_names)
# Eons are level 1, eras=level 2, periods=3, subperiods=4, epochs=5
l5s <- strat_names[which(strat_names$scale_level == 5),]
l5s[, 1:5]
if(any(ranges$max_ma > max(l5s$max_ma)))
  stop("genus ranges extend below Phanerozoic. Add older Ediacaran bins.\n")

# paleoTree's binning 'timeList' requires first matrix to contain strat interval
# ranges and second matrix to list the bins each genus occurs. The first matrix
# should have the oldest bin at top, with youngest at bottom.
l5s <- l5s[rev(seq.int(nrow(l5s))), ]
rownames(l5s) <- seq.int(nrow(l5s))
taxon.times <- ranges
colnames(taxon.times) <- c("max_bin", "min_bin")
for(t in 1:nrow(l5s)) {
  wh.max <- which(ranges$max_ma > l5s$min_ma[t] & 
                    ranges$max_ma <= l5s$max_ma[t])
  wh.min <- which(ranges$min_ma >= l5s$min_ma[t] & 
                    ranges$min_ma < l5s$max_ma[t])
  taxon.times$max_bin[wh.max] <- t
  taxon.times$min_bin[wh.min] <- t
}
timeList <- list(int.times = l5s[, 9:10], taxon.times = taxon.times)
# save(timeList, file = "timeList")
# load("timeList")



## CALCULATE SAMPLING, EXTINCTION, AND ORIGINATION RATES FOR CAL3 TIME TREE

# Use fossil record to calculate FreqRat, the frequency ratio method of Foote
# and Raup (1996) for estimating sampling probability.
pres1 <- freqRat(timeList, calcExtinction = TRUE)
pres1
# FreqRat = 0.805, an ~80% preservation probability per sampling interval, which
# is quite high. (Examples in Foote and Raup ranged from low of 0.25 for North
# American Cenozoic mammals to a high of 0.87 for Jurassic bivalves.)
# Per-interval extinction rate equals 0.25.

## Updated maximum-likelihood of the observed frequency of taxon durations (from
## Foote 1997).
likFun <- make_durationFreqDisc(timeList)
pres2 <- optim(parInit(likFun), likFun, lower = parLower(likFun),
               upper = parUpper(likFun), method = "L-BFGS-B", 
               control = list(maxit = 1000000))
pres2$par
# Per-interval taxonomic sampling probability (= R, probability of preservation
# at least once in a stratigraphic bin) equals 1, which is perfect (and probably
# not correct). Instantaneous per-capita extinction rate equals 0.29.
qsProb2Comp(R = pres2$par[2], q = pres2$par[1])
# Convert the R and extinction rate to fossil completeness: 100% of taxa are
# sampled in this clade per sampled interval, which seems unrealistic. This is
# likely caused by presence of zero-length branches.

# These methods assume the duration of intervals is constant. Is that
# approximately true? (Limiting to intervals with sampled genera)
Int.durations <-
  -apply(timeList[[1]][seq.int(max(timeList[[2]]$min_bin)), ], 1, diff)
hist(Int.durations, n = 20)
meanInt <- mean(Int.durations)
meanInt
# 6.22 Myr, approximately the same if restrict to Cambrian-Ordovician (5.74 Myr)

# Convert sampling probability to mean sampling rate:
sRate <- sProb2sRate(pres2$par[2], int.length = meanInt)
sRate
# Infinite sampling rate per lineage-million years??? Bapst and Hopkins (2017)
# note that this can occur as by-product of zero-length branches. Because this
# is clearly unrealistic, setting sRate as the value of 0.10 calculated by Foote
# and Raup (1996) for early Paleozoic crinoids (as converted to same units in
# Bapst and Hopkins 2017: p. 54, table 2). Note that the cal3 time trees
# produced below using sRate = Inf are not substantially different in their
# output.
sRate <- 0.10

# To get the extinction rate (also the branching rate), divide the extinction
# rate by the interval length.
divRate <- pres2$par[1] / meanInt
divRate # 0.047 genus extinctions per million years



## Introduction to paleoTS time trees ##########################################
set.seed(312)
timeTree <- timePaleoPhy(tree = tree, timeData = ranges)
bin_timeTree <- bin_timePaleoPhy(tree = tree, timeList = timeList)
bin_equaltimeTree <- bin_timePaleoPhy(tree = tree, timeList = timeList, 
                                      type = "equal", vartime = 1)

# Compare trees (with diversity curve)
phyloDiv(tree)
# plot(tree, cex=.1, main = "raw tree scaled to character changes"); axisPhylo()
phyloDiv(timeTree)
# plot(timeTree, cex=.1, main = "basic time tree"); axisPhylo()
phyloDiv(bin_timeTree)
# plot(bin_timeTree, cex=.1, main = "binned basic time tree"); axisPhylo()
phyloDiv(bin_equaltimeTree)
# plot(bin_timeTree, cex=.1, main = "binned equal time tree"); axisPhylo()

# Note that the cladogram and time-trees differ primarily in having different
# branch lengths (and the time-trees add a specified root time.)
identical(tree$tip.label, timeTree$tip.label)             # TRUE
identical(timeTree$tip.label, bin_timeTree$tip.label)     # TRUE
identical(tree$edge, timeTree$edge)                       # TRUE
identical(tree$edge.length, timeTree$edge.length)         # FALSE

# Observe how the ranges were modified by the binning resampling from uniform
# distribution:
head(ranges)
head(bin_timeTree$ranges.used)
head(bin_equaltimeTree$ranges.used)



## BUILD EQUAL TIME TREE #######################################################

# Using this tree because of prior common usage (e.g., Lloyd 2018, Button, et
# al., 2017; Puttick, et al., 2017; Allen, et al., 2019) and because has few
# zero-length branches (compared to "basic" tree). The "equal" method, first
# used in Brusette, et al. (2008), uses fossil occurrence FADs to date the
# first node in a phylogeny, and with zero-length branches avoided by converting
# the ZLB to a small length and reducing the branch leading to that node by a
# similar amount.

# Create 100 time trees and then generate a median diversity curve ('randres'
# randomly resolves polytomies). Using vartime = 1 million years, which is a
# small but not negligible value. (Trying different values did not alter
# results.)
set.seed(312)
equal.trees <- bin_timePaleoPhy(tree = tree, timeList = timeList, ntrees = 100, 
                               type = "equal", vartime = 1, randres = TRUE, 
                               add.term = TRUE)
beep(3)
# save(equal.trees, file = "equal.trees")
# load("equal.trees")

# Compare first 9
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for (i in 1:9) {
  plot(ladderize(equal.trees[[i]]), show.tip.label = FALSE, no.margin = TRUE)
}
par(op)
bin_multiDiv <- multiDiv(equal.trees, plot = FALSE)
plotMultiDiv(bin_multiDiv, timelims = c(550, 440))
# abline(v = l5s$max_ma) # Interval boundaries

# Create "consensus" tree. (Not sure this is a thing, but relatively easy to do
# and seems reasonable to me). The preferable way would be to use
# phytools::averageTree() but the available methods are not computationally
# tractable with our ensemble of large time trees. ls.consensus() implements the
# patristic distance method of Lapointe & Cucumel (1997). This tree is used as a
# sensitivity test to confirm that the cal3 trees and resulting chosen trees
# used in subsequent disparity analyses are relatively similar.
(t.start <- Sys.time())
consensus.equal.tree <- phytools::ls.consensus(trees = equal.trees)
(Sys.time() - t.start) # 5.3 minutes
beep(3)
# save(consensus.equal.tree, file = "consensus.equal.tree")
# load("consensus.equal.tree")

# This consensus tree lacks time-calibrated branch lengths (and other output
# from paleoTS). Use minTreeDist() to find the tree that is most similar to use
# as the consensus. (See ?treedist for descriptions of each method.) All 4
# available methods choose the same tree, so using the default quadratic
# (=weighted) path difference method of Steel and Penny (1993) that uses branch
# lengths (which are most important for distinguishing time trees, given the
# underlying phylogenetic topology is more-or-less constant). All methods are
# relatively slow, taking several hours. The resulting chosen tree is identical
# for all four methods.
(t.start <- Sys.time())
closest.equal <- phytools::minTreeDist(tree = consensus.equal.tree, 
                                       trees = equal.trees,
                                       method = "quadratic.path.difference")
(Sys.time() - t.start) # ~ 2 hours
beep(3)
# save(closest.equal, file = "closest.equal")
# load("closest.equal")
# Confirm that it is a close match (they're identical):
round(phangorn::treedist(tree1 = consensus.equal.tree, tree2 = closest.equal), 6)

# Which tree in the equal.trees sample is this most like? Here we identify
# close matches, then pick the match with the closest lineage-richness trend
# line for analyses.
sq <- 1:length(equal.trees)
treedists <- matrix(ncol = 4, nrow = max(sq))
for (s in sq) {
  treedists[s,] <- treedist(equal.trees[[s]], closest.equal)
}
bests <- apply(treedists, 2, which.min)
bests
# RF/sym = 1, BSD/KF = 71, path = 46, QPD = 90
# Use the one with most similar diversity curve:
p.1 <- phyloDiv(equal.trees[[1]], int.times = bin_multiDiv$int.times)
p.71 <- phyloDiv(equal.trees[[71]], int.times = bin_multiDiv$int.times)
p.46 <- phyloDiv(equal.trees[[46]], int.times = bin_multiDiv$int.times)
p.90 <- phyloDiv(equal.trees[[90]], int.times = bin_multiDiv$int.times)

# Plot all on same graph, limiting to Cambro-Ordovician
plotMultiDiv(bin_multiDiv, timelims = c(550, 440))
lines(p.1[, c(1, 3)], col = "yellow", lwd = 2)
lines(p.71[, c(1, 3)], col = "red", lwd = 2)
lines(p.46[, c(1, 3)], col = "blue", lwd = 2)
lines(p.90[, c(1, 3)], col = "orange", lwd = 2)

# Limit correlation analysis to Cambro-Ordovician
wh.overlap <- which(bin_multiDiv$int.times[, 1] <= 541 & 
                      bin_multiDiv$int.times[, 1] >= 443.4)
# First difference correlation coefficient (simplified because equal bins)
diff.equal <- diff(bin_multiDiv$median.curve[wh.overlap, 1])
diff.p.1 <- diff(p.1[wh.overlap, 3])
diff.p.71 <- diff(p.71[wh.overlap, 3])
diff.p.46 <- diff(p.46[wh.overlap, 3])
diff.p.90 <- diff(p.90[wh.overlap, 3])
round(cor(diff.equal, diff.p.1), 3)   # 0.949 *** most correlated
round(cor(diff.equal, diff.p.71), 3)  # 0.926
round(cor(diff.equal, diff.p.46), 3)  # 0.942
round(cor(diff.equal, diff.p.90), 3)  # 0.924


# Confirm no zero-length branches present (because Claddis' morphological
# disparity functions disallow it).
sq <- 1:length(equal.trees)
# How many have ZLBs?
ZLBs <- sapply(sq, function(sq) any(equal.trees[[sq]]$edge.length == 0))
table(ZLBs) # None do!

# Confirm that the time tree producing the most correlated phylogenetic
# diversity is a close match to the consensus tree
round(phangorn::treedist(tree1 = consensus.equal.tree, tree2 = equal.trees[[1]]), 6)

# While at it, what's the distribution of estimated root ages?
roots <- sapply(sq, function(sq) equal.trees[[sq]]$root.time)
summary(roots)                     # Mean age is early Cambrian
length(roots < max(ranges$max_ma)) # with all 100 younger than oldest FAD, 
# b/c paleotree "equal" method assumes this
hist(roots)

# Save equal tree for Claddis analyses
equal.tree <- equal.trees[[1]]
phyloDiv(equal.tree)
# save(equal.tree, file = "equal.tree")
# load("equal.tree")
# pdf(file = "EchinoPhylogeny_equal.pdf", height = 20)
# plot(equal.tree, cex = 0.4); axisPhylo()
# dev.off()







## BUILD CAL3 TIME TREE ########################################################
set.seed(312)
(t.start <- Sys.time())
cal3tree <- bin_cal3TimePaleoPhy(tree, timeList, brRate = divRate, 
                                 extRate = divRate, sampRate = sRate, 
                                 ntrees = 1, plot = TRUE)
Sys.time() - t.start
beep(3)
phyloDiv(cal3tree)

# Build 100 in parallel

# (Was not done in summer 2020, but if need to run again, refer to
# 6-Phylogenetic inertia.R for better way to implement load-balancing where the
# L-Ecuyer RNG seed can allow replication.)
set.seed(312)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
(t.start <- Sys.time())
nreps <- 100
par.cal3trees <- foreach(i = 1:nreps, .packages = "paleotree") %dopar% {
  bin_cal3TimePaleoPhy(tree, timeList, brRate = divRate, extRate = divRate, 
                       sampRate = sRate, ntrees = 1, plot = FALSE)
}
Sys.time() - t.start
stopCluster(cl)
beep(3)

# Save trees
# save(par.cal3trees, file = "par.cal3trees")
# load("par.cal3trees")

# Plot first tree
phyloDiv(par.cal3trees[[1]])

# Compare first 9
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for (i in 1:9) {
  plot(ladderize(par.cal3trees[[i]]), show.tip.label = FALSE, no.margin = TRUE)
}
par(op)

# Plot median diversity curve with 95%iles
cal3_multiDiv <- multiDiv(par.cal3trees, plot = FALSE)
plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))
abline(v = l5s$max_ma) # Interval boundaries

# Create "consensus" cal3 tree. (See comments above.) 
(t.start <- Sys.time())
consensus.cal3.tree <- phytools::ls.consensus(trees = par.cal3trees)
(Sys.time() - t.start)
beep(3)
# save(consensus.cal3.tree, file = "consensus.cal3.tree")
# load("consensus.cal3.tree")

# This consensus tree lacks time-calibrated branch lengths (and other output
# from paleoTS). Use minTreeDist() to find the tree that is most similar to use
# as the consensus. (See ?treedist for descriptions of each method.) All 4
# available methods choose the same tree, so using the default quadratic
# (=weighted) path difference method of Steel and Penny (1993) that uses branch
# lengths. The method is relatively slow (symmetric = 1 hr, path = 1.3 hrs,
# quadratic = 1.6 hrs, and branch = 4.3 hrs). The resulting chosen trees are
# identical.
(t.start <- Sys.time())
closest.cal3 <- phytools::minTreeDist(tree = consensus.cal3.tree, 
                                      trees = par.cal3trees, 
                                      method = "quadratic.path.difference")
(Sys.time() - t.start)
beep(3)
# save(closest.cal3, file = "closest.cal3")
# load("closest.cal3")
# Confirm that it is a close match (they're identical):
round(phangorn::treedist(tree1 = consensus.cal3.tree, tree2 = closest.cal3), 6)

# Which tree in the par.cal3trees sample is this most like? Here we identify
# close matches, then pick the match with the closest lineage-richness trend
# line for analyses.
sq <- 1:length(par.cal3trees)
treedists <- matrix(ncol = 4, nrow = max(sq))
for (s in sq) {
  treedists[s,] <- treedist(par.cal3trees[[s]], closest.cal3)
}
bests <- apply(treedists, 2, which.min)
bests
# RF/sym = 68, BSD/KF = 2, path = 31, QPD = 80
# Use the one with most similar diversity curve:
p.68 <- phyloDiv(par.cal3trees[[68]], int.times = cal3_multiDiv$int.times)
p.2 <- phyloDiv(par.cal3trees[[2]], int.times = cal3_multiDiv$int.times)
p.31 <- phyloDiv(par.cal3trees[[31]], int.times = cal3_multiDiv$int.times)
p.80 <- phyloDiv(par.cal3trees[[80]], int.times = cal3_multiDiv$int.times)

# Plot all on same graph, limiting to Cambro-Ordovician
plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))
lines(p.68[, c(1, 3)], col = "yellow", lwd = 2)
lines(p.2[, c(1, 3)], col = "red", lwd = 2)
lines(p.31[, c(1, 3)], col = "blue", lwd = 2)
lines(p.80[, c(1, 3)], col = "orange", lwd = 2)

# Limit correlation analysis to Cambro-Ordovician
wh.overlap <- which(cal3_multiDiv$int.times[, 1] <= 541 & 
                      cal3_multiDiv$int.times[, 1] >= 443.4)
# First difference correlation coefficient (simplified because equal bins)
diff.cal3 <- diff(cal3_multiDiv$median.curve[wh.overlap, 1])
diff.p.68 <- diff(p.68[wh.overlap, 3])
diff.p.2 <- diff(p.2[wh.overlap, 3])
diff.p.31 <- diff(p.31[wh.overlap, 3])
diff.p.80 <- diff(p.80[wh.overlap, 3])
round(cor(diff.cal3, diff.p.68), 3)  #0.927 *** most correlated
round(cor(diff.cal3, diff.p.2), 3)   #0.916
round(cor(diff.cal3, diff.p.31), 3)  #0.923
round(cor(diff.cal3, diff.p.80), 3)  #0.916

# Confirm that it is a close match
round(phangorn::treedist(tree1 = consensus.cal3.tree, 
                         tree2 = par.cal3trees[[68]]), 6)

# While at it, what's the distribution of estimated root ages?
roots <- sapply(sq, function(sq) par.cal3trees[[sq]]$root.time)
summary(roots) # Mean age is very late Ediacaran, but 25% within Cambrian
hist(roots)

# Save cal3 tree for Claddis analyses
cal3.tree <- par.cal3trees[[68]]
phyloDiv(cal3.tree)
# save(cal3.tree, file = "cal3.tree")
# load("cal3.tree")
# pdf(file = "EchinoPhylogeny_cal3.pdf", height = 20)
# plot(cal3.tree, cex = 0.4); axisPhylo()
# dev.off()



## EXPLORE POPULATION OF TIME TREES & CHOOSE MOST TYPICAL TREE FOR ANALYSES ####
# How different are the 100 cal3 trees?
class(par.cal3trees) <- "multiPhylo"
# Force to treat as 'multiPhylo' b/c built built in-parallel earlier
dist1 <- path.dist(par.cal3trees)
summary(as.vector(dist1))
hist(dist1)

# How different are the 100 "equal" trees?
dist2 <- path.dist(equal.trees)
summary(as.vector(dist2))
hist(dist2)

# How different are the "equal" trees from the "cal3" trees (here using all
# pairs)?
lcal3 <- length(par.cal3trees)
lequal <- length(equal.trees)
dist3 <- matrix(nrow = lcal3, ncol = lequal)
for(r in 1:lcal3){
  for(c in 1:lequal){
  dist3[r, c] <- path.dist(tree1 = par.cal3trees[[r]], tree2 = equal.trees[[c]])  
  }
}
hist(dist3)

# How different are the "consensus" trees from the "equal" and cal3 algorithms?
treedist(tree1 = consensus.equal.tree, tree2 = consensus.cal3.tree)

# How different are the "consensus" trees from the "equal" and cal3 algorithms?
treedist(tree1 = equal.tree, tree2 = cal3.tree)
# Quite similar (slightly less than the consensus trees overall)


# Plot in joint histograms
comb <- c(dist1, dist2, dist3)
brks <- pretty(comb, n = 20)
hist(dist1, breaks = brks, prob = TRUE, border = "white", col = "darkgray")
hist(dist2, breaks = brks, prob = TRUE, add = TRUE, border = "black", 
     col = "transparent")
hist(dist3, breaks = brks, prob = TRUE, add = TRUE, border = "blue", 
     col = "transparent")
abline(v = path.dist(tree1 = consensus.cal3.tree, tree2 = consensus.equal.tree))
abline(v = path.dist(tree1 = cal3.tree, tree2 = equal.tree), lty = 2)
legend("topright", inset = .1, c("cal3-cal3", "equal-equal", "cal3-equal"), 
       pch = c(22, 22, 22), pt.bg = c("darkgray", "transparent", "transparent"), 
       col = c("white", "black", "blue"), cex = 1.5, pt.cex = 2)

# Analysis of variance to test for differences
group <- factor(c(rep("cal3", length(dist1)), rep("equal", length(dist2)), 
                  rep("both", length(dist3))))
boxplot(log(comb) ~ group)
aov <- aov(log(comb) ~ group)
summary(aov)

# CONCLUSION: The cal3 trees are quite self-similar, as are the equal trees
# (although slightly more variable). The typical equal tree is slightly more
# distinct than the typical cal3, but they all overlap significantly. However,
# each population is statistically different. The average tree in each sample is
# somewhat dissimilar. However, the selected equal and cal3 trees chosen for
# subsequent disparity analyses are quite similar to one another and to their
# respective time-tree sample.

# What most matters for our hypothesis is whether the lineage richness trends
# are distinct. Let's compare the median lineage richness curves from the
# samples of time trees
plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))
int.start <- bin_multiDiv$int.times[, 1]
int.end <- bin_multiDiv$int.times[, 2]
times <- apply(rbind(int.start, int.end), 2, mean)
lines(x = times, y = bin_multiDiv$median.curve[, 1], col = "red", lwd = 2)
abline(v = l5s$max_ma) # Interval boundaries
legend("topleft", inset = .1, c("cal3", "equal"), pch = c("-", "-"), 
       col = c("black", "red"), cex = 1.5, pt.cex = 5)

# Confirm statistically that they're not different
wh.overlap <-
  which(round(cal3_multiDiv$int.times[, 1]) <= round(max(bin_multiDiv$int.times[, 1])) &
        round(cal3_multiDiv$int.times[, 1]) >= round(min(bin_multiDiv$int.times[, 1])))
# First difference correlation coefficient (simplified because equal bins)
diff.cal3 <- diff(cal3_multiDiv$median.curve[wh.overlap, 1])
diff.equal <- diff(bin_multiDiv$median.curve[, 1])
cor(diff.cal3, diff.equal) # 0.91
lm <- lm(diff.equal ~ diff.cal3)
plot(diff.cal3, diff.equal, xlab = "change in cal3 diversity", 
     ylab = "change in cal3 diversity")
abline(lm)
summary(lm) # p < 2.2e-16 *** significantly correlated

# Absolute difference?
abs.diff <- bin_multiDiv$median.curve[, 1] - cal3_multiDiv$median.curve[wh.overlap, 1]
summary(abs.diff) # ~ 10% greater, on average

# RESULTS: Both time-tree types are statistically very similar in all ways
# (topologically, phylogenetic diversity, and phylogenetic distance). However,
# the resulting diversity curve for the equal treatment produces a slightly
# (~10%) greater net lineage richness during most time intervals because of the
# more drawn-out diversification interval. Because subsequent analyses do not
# work with zero-length branches and unresolved polytomies, we use the 'equal'
# tree for all subsequent analyses.
