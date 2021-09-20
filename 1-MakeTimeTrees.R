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
setwd("C:/Users/pnovack-gottshall/OneDrive - Benedictine University/Documents/Manuscripts/CamOrdEchinos")
# setwd("[filepath to folder containing data files on your personal machine]")

library(paleotree)  # v. 3.3.25
library(phytools)   # v. 0.7-83
library(phangorn)   # v. 2.7.0
library(beepr)      # v. 1.3
library(doParallel) # v. 1.0.16

## Import cladogram in Nexus format
tree <- read.nexus("echinodermtree.nex")
if (class(tree) != "phylo") stop("Nexus file is not class 'phylo")
# pdf(file = "EchinoPhylogeny_supertree.pdf", height = 20)
# plot(tree, cex = 0.4); axisPhylo()
# dev.off()


## Process stratigraphic ranges

# Import strat ranges (previously downloaded fom PBDB June 23, 2020)
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
# likely caused by presence of zero-length branches (ZLBs).

# These methods assume the duration of intervals is constant. Is that
# approximately true? (Limiting to intervals with sampled genera)
Int.durations <-
  -apply(timeList[[1]][seq.int(max(timeList[[2]]$min_bin)), ], 1, diff)
hist(Int.durations, n = 20)
meanInt <- mean(Int.durations)
meanInt
# 6.22 Myr / epoch for Phanerozoic; approximately the same (5.74 Myr/epoch) if
# restrict to Cambrian-Ordovician

# Convert sampling probability to mean sampling rate:
sRate <- sProb2sRate(pres2$par[2], int.length = meanInt)
sRate
# Infinite sampling rate per lineage-million years??? Bapst and Hopkins (2017)
# note that this can occur as by-product of ZLBs. Because this is clearly
# unrealistic, setting sRate as the value of 0.10 calculated by Foote and Raup
# (1996) for early Paleozoic crinoids (as converted to same units in Bapst and
# Hopkins 2017: p. 54, table 2). Note that the cal3 time trees produced below
# using sRate = Inf are not substantially different in their output.
sRate <- 0.10

# To get the extinction rate (also the branching rate), divide the extinction
# rate by the interval length.
divRate <- pres2$par[1] / meanInt
divRate # 0.047 genus extinctions per million years




## PREPARING THE CAL3 TIME-SCALED TREES ########################################

# The 'cal3' time-scaling algorithm was developed by David Bapst (Bapst, 2013,
# 2014; Bapst, et al., 2012; Bapst and Hopkins, 2016). It takes into account not
# just stratigraphic ranges of tip taxa but also rates of origination,
# extinction, and sampling (preservation). Sensitivity analyses (e.g., Bapst,
# 2014; Bapst and Hopkins, 2016) demonstrate that the 'cal3' consistently ranks
# as good as or better than other available time-scaling methods.

# Other time-scaling methods (e.g., the 'equal' or 'basic' algorithms) are
# inappropriate here because the presence of taxon-rich and phylogenetically
# not-yet-studied higher taxa artificially inflates per-origination diversity by
# pushing their subclade roots back in time by arbitrary, user-defined steps (as
# opposed to using the actual stratigraphic ranges and
# sampling/origination/extinction probabilities drawn from the data themselves).

# Because the stratigraphic ranges we have available are based on discrete bins,
# we use the paleotree::bin_cal3TimePaleoPhy() function to build cal3 trees.
# Some subsequent analyses are intractable with zero-length branches (ZLBs) and
# polytomies. The 'cal3' algorithm inherently resolves polytomies during the
# time-scaling algorithm. ZLBs are only present in the 'cal3' tree within nodes;
# no terminal branch tips have ZLBs. For analyses that can not handle ZLBs,
# these are removed as a secondary step, by replacing them with a negligible
# (0.001 Myr) length, and the root of the tree is adjusted accordingly.

# We implement dataTreatment = 'firstLast' (the default treatment) because genus
# ranges are resolved to epochs. (David Bapst, e-mail 8/3/2021, confirms this is
# the appropriate choice for our data.) Sensitivity analyses comparing the
# 'firstLast' algorithm to the 'randObs' show negligible differences in terms of
# resulting phylogenetic diversity, ZLBs, and time-scaled ranges. FAD.only =
# FALSE (default) is set so that FADs are used for rootward node ages and the
# tip ages are set as LADs. In this way, the 'cal3' functions always add
# terminal ranges to taxa, allowing the time-scaled ranges to be used to
# estimate phylogenetic diversity.

# Because the 'cal3' algorithm involves stochastic processes, we use 50
# stochastic trees to evaluate variability of subsequent results to the tree
# structure. Because we have no a priori knowledge of whether any taxa represent
# the ancestors for tips, we also evaluate time-scaled trees built using anc.wt
# = 1 and anc.wt = 0, to allow to both possibilities, building 25 of each tree
# type.

# The 'cal3' algorithm requires probabilities of per-capita origination
# (=branching), extinction, and preservation. We estimate per-interval taxonomic
# origination and extinction rates using the built-in functions within
# 'paleoTree' (as coded above, using the raw stratigraphic ranges and
# non-time-scaled phylogeny). These functions estimate preservation probability
# for these Cambrian-Ordovician echinoderms to equal 0.805 using the FreqRat
# method (Foote and Raup, 1996), a high value among fossil invertebrates, and
# implies a rather complete fossil record. The improved maximum-likelihood
# method (Foote, 1997) yields unrealistically high sampling estimates because of
# zero-length branches in the raw form, so a crinoid-wide sampling rate from
# (Foote and Raup, 1996), converted to appropriate units (c.f., Bapst and
# Hopkins, 2016), was used for the 'cal3' trees. Using either sampling rate
# value yielded similar tree structures.

# For visualizations of ordination spaces through time, convergence dynamics,
# and other visual examples, we use on a single time-scaled tree that most
# represents the typical cal3 tree, using phytools::ls.consensus to build a
# consensus tree and then using various tests to identify the single tree among
# the 60 that is most similar to this consensus. (See below for details.)


# Function to replace ZLBs with negligible branch length. Also accordingly fixes
# the root time of the tree.
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


# Function to calculate time-scaled ranges for tips and nodes. Useful for
# downstream analyses where require identifying which lineages are present in a
# given time interval. Also useful for calculating phylogenetic lineage richess.
# Sets the root FAD = LAD as the tree root age. Modified from code written by
# David Bapst. Thanks, Dave!
strat.ranges <- function(tree) {
  ts.ranges <- tree$edge
  ntime <- tree$root.time - ape::node.depth.edgelength(tree)
  ts.ranges[, 1] <- ntime[ts.ranges[, 1]]
  ts.ranges[, 2] <- ntime[ts.ranges[, 2]]
  terminalEdges <- which(tree$edge[, 2] <= Ntip(tree))
  row.names(ts.ranges) <- as.character(tree$edge[, 2])
  row.names(ts.ranges)[terminalEdges] <-
    tree$tip.label[tree$edge[terminalEdges, 2]]
  # Change order so matches input data rownames
  br.order <-
    as.character(c(row.names(tree$ranges.used), (1 + Ntip(tree)):(Ntip(tree) + Nnode(tree))))
  ts.ranges <- ts.ranges[match(br.order, row.names(ts.ranges)),]
  # Manually set root FAD = LAD
  ts.ranges[(1 + Ntip(tree)),] <- rep(tree$root.time, 2)
  row.names(ts.ranges)[1 + Ntip(tree)] <-
    as.character((1 + Ntip(tree)))
  return(ts.ranges)
}

## BUILD A SINGLE CAL3 TREE ####################################################
set.seed(17)
(t.start <- Sys.time())
cal3tree <- cal3tree.noZLBs <- NULL
cal3tree <- bin_cal3TimePaleoPhy(tree, timeList, brRate = divRate, 
                                 extRate = divRate, sampRate = sRate, 
                                 dateTreatment = "firstLast", FAD.only = FALSE, 
                                 anc.wt = 0, randres = FALSE, ntrees = 1, 
                                 plot = FALSE)
Sys.time() - t.start # 1.73 minutes
beep(3)

# Visualizations and tests

# 1. Are polytomies present?
!is.binary.phylo(cal3tree) # No polytomies

# 2. How many ZLBs?
length(which(cal3tree$edge.length == 0)) # 57 ZLBs

# 3. Compare pre and post FAD/LAD ranges.
tail(ranges)
(nr <- strat.ranges(cal3tree))[361:370,]
par(mfrow = c(1, 2))
hist(ranges[, 1] - nr[1:366, 1], main = "FADs")
hist(ranges[, 2] - nr[1:366, 2], main = "LADs")

# 4. Phylogenetic diversity
phyloDiv(cal3tree)

# 5. Effect of replacing ZLBs with negligible branch lengths (and adjusting root
# accordingly)

# Shortest branch length
min(cal3tree$edge.length[cal3tree$edge.length > 0]) # minimum is 0.1

# Replace ZLBs with 0.001
cal3tree.noZLBs <- replace.ZLBs(cal3tree)
table(compareTermBranches(cal3tree, cal3tree.noZLBs))
table(round(compareNodeAges(cal3tree, cal3tree.noZLBs), 2))
# No changes to tips and negligible changes to nodes

# 6. Effect on root age adjustment
cal3tree$root.time - cal3tree.noZLBs$root.time # root moved backward negligible 0.003 Myr



## BUILD 50 CAL3 TREES IN PARALLEL #############################################

# First make 25 with anc.wt = 0
cl <- makeCluster(detectCores())
registerDoParallel(cl)
# Use load-balancing because of different run times for the MCMC optimizations
opts <- list(preschedule = FALSE)
clusterSetRNGStream(cl, 3142)  # Set L-Ecuyer RNG seed
(t.start <- Sys.time())
nreps <- 25
cal3trees.0anc <- foreach(i = 1:nreps, .options.snow = opts, 
                              .packages = "paleotree") %dopar% {
  bin_cal3TimePaleoPhy(tree, timeList, brRate = divRate, extRate = divRate, 
                       sampRate = sRate, dateTreatment = "firstLast", 
                       FAD.only = FALSE, anc.wt = 0, randres = FALSE, 
                       ntrees = 1, plot = FALSE)
                              }
Sys.time() - t.start # 6.7 minutes on 8-core laptop
stopCluster(cl)
beep(3)

# Second make 25 with anc.wt = 1
cl <- makeCluster(detectCores())
registerDoParallel(cl)
# Use load-balancing because of different run times for the MCMC optimizations
opts <- list(preschedule = FALSE)
clusterSetRNGStream(cl, 1234)  # Set L-Ecuyer RNG seed
(t.start <- Sys.time())
nreps <- 25
cal3trees.1anc <- foreach(i = 1:nreps, .options.snow = opts, 
                              .packages = "paleotree") %dopar% {
  bin_cal3TimePaleoPhy(tree, timeList, brRate = divRate, extRate = divRate, 
                       sampRate = sRate, dateTreatment = "firstLast", 
                       FAD.only = FALSE, anc.wt = 0, randres = FALSE, 
                       ntrees = 1, plot = FALSE)
  }
Sys.time() - t.start # 6.5 minutes on 8-core laptop
stopCluster(cl)
beep(3)

# Combine into a single forest of frees
cal3trees <- c(cal3trees.0anc, cal3trees.1anc)

# Save trees
# save(cal3trees, file = "cal3trees")
# load("cal3trees")



## EXPLORE POPULATION OF TIME TREES ############################################
# Plot first tree
phyloDiv(cal3trees[[1]])

# Compare first 5 and last 4
par(mfrow = c(3, 3), mar = c(1, 1, 1, 1))
for (i in c(1:5, 47:50)) {
  plot(ladderize(cal3trees[[i]]), show.tip.label = FALSE, no.margin = TRUE)
}
par(op)

# Plot median diversity curve with 95%iles
cal3_multiDiv <- multiDiv(cal3trees, plot = FALSE)
plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))
abline(v = l5s$max_ma, col = "darkgray") # Interval boundaries

# Note the relatively little variation in diversity curves for these 50 cal3
# trees.

# Are polytomies present?
table(!is.binary.multiPhylo(cal3trees)) # No polytomies!

# How many zero-length branches (ZLBs)?
sq <- 1:length(cal3trees)
summary(sapply(sq, function(sq) 
  length(which(cal3trees[[sq]]$edge.length == 0)))) # 43 - 78 ZLBs, mean = 61

# Compare pre and post FAD/LAD ranges.
par(mfrow = c(1, 2))
tail(ranges)
(nr <- strat.ranges(cal3trees[[1]]))[361:370,]
hist(ranges[, 1] - nr[1:366, 1], main = "FAD")
hist(ranges[, 2] - nr[1:366, 2], main = "LAD")
par(op)

# Effect of removing ZLBs on terminal and node branch lengths
cal3trees.noZLBs <- lapply(cal3trees, replace.ZLBs)
summary(br <- as.vector(sapply(sq, function(sq) 
  compareTermBranches(cal3trees[[sq]], cal3trees.noZLBs[[sq]]))))
summary(nd <- as.vector(sapply(sq, function(sq) 
  compareNodeAges(cal3trees[[sq]], cal3trees.noZLBs[[sq]]))))
100 * round(table(round(nd, 3)) / length(nd), 3)
hist(nd)
# 0 branch adjustments; -0.010 - + 0.017 (median = 0.001) node adjustments

# Effect on root age adjustment
summary(root.adj <- sapply(sq, function(sq) 
  cal3trees[[sq]]$root.time - cal3trees.noZLBs[[sq]]$root.time))
hist(root.adj)
# root adjusted -0.010 - 0 Myr, median = -0.003 Myr



## CHOOSE MOST TYPICAL TREE FOR ANALYSES #######################################
# See comments above for rationale.

# Create "consensus" cal3 tree. (See comments above.) 
(t.start <- Sys.time())
consensus.cal3.tree <- phytools::ls.consensus(trees = cal3trees)
(Sys.time() - t.start) # 15 minutes
beep(3)
# save(consensus.cal3.tree, file = "consensus.cal3.tree")
# load("consensus.cal3.tree")

# This consensus tree lacks time-calibrated branch lengths (and other output
# from paleoTS). Use minTreeDist() to find the tree that is most similar to use
# as the consensus. (See ?treedist for descriptions of each method.) All 4
# available methods choose the same tree, so using the default quadratic
# (=weighted) path difference method of Steel and Penny (1993) that uses branch
# lengths. The method is relatively slow (symmetric = 20 min., path = 23 min.,
# quadratic = 51 min., and branch = 1.5 hrs). The resulting chosen trees are
# identical in all cases, except for the branch.score.difference.
(t.start <- Sys.time())
closest.cal3 <- phytools::minTreeDist(tree = consensus.cal3.tree, 
                                      trees = cal3trees, 
                                      method = "quadratic.path.difference")
(Sys.time() - t.start)
beep(3)
# save(closest.cal3, file = "closest.cal3")
# load("closest.cal3")
# Confirm that it is a close match (they're identical):
round(phangorn::treedist(tree1 = consensus.cal3.tree, tree2 = closest.cal3), 6)


# Which tree in the cal3trees sample is this most like? Here we identify
# close matches, then pick the match with the closest lineage-richness trend
# line for analyses.
sq <- 1:length(cal3trees)
treedists <- matrix(ncol = 4, nrow = max(sq))
for (s in sq) {
  treedists[s, ] <- treedist(cal3trees[[s]], closest.cal3)
}
bests <- apply(treedists, 2, which.min)
bests
# RF/sym = 5, BSD/KF = 43, path = 50, QPD = 10
treedists[bests, ]

# Use the one with most similar diversity curve:
p.5 <- phyloDiv(cal3trees[[5]], int.times = cal3_multiDiv$int.times)
p.43 <- phyloDiv(cal3trees[[43]], int.times = cal3_multiDiv$int.times)
p.50 <- phyloDiv(cal3trees[[50]], int.times = cal3_multiDiv$int.times)
p.10 <- phyloDiv(cal3trees[[10]], int.times = cal3_multiDiv$int.times)

# Plot all on same graph, limiting to Cambro-Ordovician
plotMultiDiv(cal3_multiDiv, timelims = c(550, 440))
lines(p.5[, c(1, 3)], col = "yellow", lwd = 2)
lines(p.43[, c(1, 3)], col = "red", lwd = 2)
lines(p.50[, c(1, 3)], col = "blue", lwd = 2)
lines(p.10[, c(1, 3)], col = "orange", lwd = 2)

# Limit correlation analysis to Cambro-Ordovician
wh.overlap <- which(cal3_multiDiv$int.times[, 1] <= 541 & 
                      cal3_multiDiv$int.times[, 1] >= 443.4)
# First difference correlation coefficient (simplified because equal bins)
diff.cal3 <- diff(cal3_multiDiv$median.curve[wh.overlap, 1])
diff.p.5 <- diff(p.5[wh.overlap, 3])
diff.p.43 <- diff(p.43[wh.overlap, 3])
diff.p.50 <- diff(p.50[wh.overlap, 3])
diff.p.10 <- diff(p.10[wh.overlap, 3])
round(cor(diff.cal3, diff.p.5), 3)  #0.931
round(cor(diff.cal3, diff.p.43), 3)   #0.927
round(cor(diff.cal3, diff.p.50), 3)  #0.936 *** most correlated
round(cor(diff.cal3, diff.p.10), 3)  #0.894

# Confirm that it is a close match
round(phangorn::treedist(tree1 = consensus.cal3.tree, 
                         tree2 = cal3trees[[50]]), 6)

# While at it, what's the distribution of estimated root ages?
roots <- sapply(sq, function(sq) cal3trees[[sq]]$root.time)
summary(roots) # Mean age is very late Ediacaran, but 25% within Cambrian
hist(roots)

# Save cal3 tree for later visualizations
cal3.tree <- cal3trees[[50]]
phyloDiv(cal3.tree)
# save(cal3.tree, file = "cal3.tree")
# load("cal3.tree")
# pdf(file = "EchinoPhylogeny_cal3.pdf", height = 20)
# plot(cal3.tree, cex = 0.4); axisPhylo()
# dev.off()



# How different are the 50 cal3 trees?
class(cal3trees) <- "multiPhylo"
# Force to treat as 'multiPhylo' b/c built built in-parallel earlier
dist1 <- path.dist(cal3trees)
summary(as.vector(dist1))
hist(dist1)

# CONCLUSION: The cal3 trees are quite self-similar, despite using two ancestral
# models.
