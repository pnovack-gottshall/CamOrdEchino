## NO LONGER USED R SCRIPTS
# See original files in "Old R script" for original workflow


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





