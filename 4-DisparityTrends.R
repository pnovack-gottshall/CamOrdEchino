## CALCULATE DISPARITY TRENDS ##################################################

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' and run 3-DisparityDistances.R to calculate disparity distance
# matrices. The code below is modified from that used in GSA 2019 analyses (ref
# below), and uses a modification of calc_metrics() in 'ecospace' package.

# Novack-Gottshall, P. M., A. Sultan, J. N. Purcell, I. Ranjha, and B. Deline.
# 2019. Phylogenetic constraint and ecological opportunity in the Cambrian and
# Ordovician radiation of echinoderms. Geological Society of America Annual
# Meeting, GSA Abstracts with Programs 51.

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 7/29/2021, v.
# 0.6.3. Users accustomed with earlier versions will need to either download the
# archived version of the package from GitHub or alter the code accordingly to
# use appropriate argument and function names.


## 1 - PREPARATIONS ############################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")
setwd("~/Manuscripts/CamOrdEchinos/Data files/NA Reformatted")

# Load packages
library(beepr)      # v. 1.3
library(paleotree)  # v. 3.3.25
library(ape)        # v. 5.5
library(geoscale)   # v. 2.0
library(ade4)       # v. 1.7-17
library(geometry)   # v. 0.4.5
library(FD)         # v. 1.0-12
library(doParallel) # v. 1.0.16
library(phytools)   # v. 0.7-83
library(vegan)      # v. 2.5-7
library(plotrix)    # v. 3.8-1
library(viridisLite)# v. 0.4.0
library(Claddis)    # v. 0.6.3 - Check SI for Lloyd 2018 for walk-through on code functions
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub or CRAN.\n")

# Modification of geoscale::geoscalePlot to allow ICS 2020 timescale 
source("~/Manuscripts/CamOrdEchinos/geoscalePlot2.R")

# Modification of FD::dbFD() to simplify calculation of FD metrics, to allow
# standardization of PCoA vector coordinates by relative magnitude of
# eigenvalues, and to directly use pre-calculated distance matrices and PCoA
# output.
source("~/Manuscripts/CamOrdEchinos/simple.dbFD.R")

# Modification of ecospace::calc_metrics. See source file for modifications.
source("~/Manuscripts/CamOrdEchinos/calc_metrics2.R")




## 2- IMPORT FILES #############################################################

# Note because the cal3trees were modified after creation (to remove zero-length
# branches), using the tree objects appended to the Claddis objects below.

# Import ancestral states from 2-InferancestralStates.R
load("mode.anc")
load("constant.anc")
load("raw.anc")
load("morph.anc")

# Import distance matrices from 3-DisparityDistances.R
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")



## 3 - GET EPOCH TIME BINS #####################################################
# strat_names <- read.csv("https://www.paleobiodb.org/data1.2/intervals/list.csv?all_records&vocab=pbdb")
strat_names <- read.csv("~/Manuscripts/CamOrdEchinos/strat_names.csv")
# Eons are level 1, eras = level 2, periods = 3, epochs = 4, ages = 5
ages <- strat_names[which(strat_names$scale_level == 5), ]
# Limit to Cambrian and Ordovician
ages <- ages[which(ages$max_ma > 444), ]
# Add Ediacaran for any pre-Cambrian nodes
ages <-
  rbind(ages, strat_names[which(strat_names$interval_name == "Ediacaran"),])
int.times <- ages[ ,9:10]  # Time bins used for disparity analyses below.
mids <- apply(ages[ ,9:10], 1, mean)
summary(int.times$max_ma - int.times$min_ma)
# Replace outdated Series 3 with Miaolingian (if using epochs)
ages$interval_name <-
  replace(ages$interval_name, which(ages$interval_name == "Series 3"), 
          "Miaolingian")
ages[, c(5, 9:10)]

# Import ICS 2020 timescale to use in plotting
ICS2020 <- read.csv("~/Manuscripts/CamOrdEchinos/timescales2020.csv", 
                    stringsAsFactors =  TRUE)
head(ICS2020)
# Note that only using the 2020 ICS timescale for plotting and not analyses,
# because the strat ranges in PBDB (and 'geoscale') are still using the 2015 ICS
# ages. The primary change is the naming of Series 3 as Miaolingian.






## GET STRAT RANGES (FADs AND LADs) FOR TIPS AND NODES #########################

# Based on code written by Dave Bapst. Root node set as root age, with FAD =
# LAD.

# All the 'X.anc' objects contain the same time-scaled trees as appended
# objects. Using the ecological objects because they have smaller object sizes
# (which enhances memory allocation).

ranges <- vector("list", length(mode.anc))
for(t in 1:length(mode.anc)) {
  tree <- mode.anc[[t]]$topper$tree
  tree.ranges <- tree$edge
  # Use time-scaled tree branch lengths to get relative node ages, and use root
  # time to scale to absolute time
  ntime <- tree$root.time - ape::node.depth.edgelength(tree)
  tree.ranges[, 1] <- ntime[tree.ranges[, 1]]
  tree.ranges[, 2] <- ntime[tree.ranges[, 2]]
  # Now assign names for each lineage as row names
  terminalEdges <- which(tree$edge[, 2] <= Ntip(tree))
  colnames(tree.ranges) <- c("FAD", "LAD")
  rownames(tree.ranges) <- as.character(tree$edge[, 2])
  rownames(tree.ranges)[terminalEdges] <- tree$tip.label[tree$edge[terminalEdges, 2]]
  # Add root node with FAD = LAD = root.time
  tree.ranges <- rbind(tree.ranges, rep(tree$root.time, 2))
  rownames(tree.ranges)[Nnode(tree) + Ntip(tree)] <- as.character(Ntip(tree) + 1)
  # Reorder so matches order of the cladistic matrix
  tree.ranges <-
    tree.ranges[match(rownames(mode.anc[[t]]$matrix_1$matrix), rownames(tree.ranges)),]
  ranges[[t]] <- tree.ranges
}
# Save for later use
# save(ranges, file = "ranges")

# Observe mean ranges for tips and nodes
sq <- 1:length(ranges)
mean.FADs <-
  apply(simplify2array(lapply(sq, function(sq) ranges[[sq]][, "FAD"])), 1, mean)
mean.LADs <-
  apply(simplify2array(lapply(sq, function(sq) ranges[[sq]][, "LAD"])), 1, mean)
mean.ranges <- cbind(FAD = mean.FADs, LAD = mean.LADs)

head(mean.ranges)
tail(mean.ranges)

# Export time-scaled ranges, converting to .csv
# write.csv(mean.ranges, file = "TSRanges.csv", row.names = TRUE)

# Which are the 30 oldest taxa?
mean.ranges[tail(order(mean.ranges[, "FAD"]), 30), ]
# On average, the 10 oldest are all ancestral nodes. Helicoplacoids Polyplacus,
# Helicoplacus, and Waucobella are the oldest tips, followed by eocrinoids
# Felbabkacystis and Rhopalocystis, Helicocystis of uncertain affinity, stem
# solute Coleicarpus, and edrioasteroid Camptostroma.

# Stratigraphic range lengths
diffs <- mean.ranges[, "FAD"] - mean.ranges[, "LAD"]
hist(diffs, 50)
summary(diffs)     # median = 8 Myr, min = ~0 Kyr, max = 130 Myr
which(diffs == 0)  # Basal echinoderm (node 367)
which(diffs < .5)  # 13 stem lineages
which.max(diffs)   # Parisocrinus, a Late Ord-Mississippian crinoid

# Explore ancestral reconstructions (using tree #50)
LHs <- apply(mode.anc[[50]]$matrix_1$matrix, 1, paste, collapse="")
mode.anc[[50]]$matrix_1$matrix[which(LHs == LHs[367]), ]
# Stem echinoderm reconstructed as relatively small (0.1 - 1 cubic cm),
# epifaunal, unattached, intermittently mobile, microbivorous filter feeder with
# low filtration density feeding elements, atop a soft substrate. (Same life
# habit as found in stem-echinoderm Ctenoimbricata, cinctans Davidocinctus,
# Elliptocinctus, Graciacystis, Gyrocystis, Lignanicystis, Ludwigicinctus,
# Nelegerocystis, Protocinctus, Rozanovicystis, Sucocystis, Trochocystites,
# Trochocystoides, and Undatacinctus, solutes Castericystis, Pahvanticystis, and
# Scalenocystis, and mitrates Aspidocarpus, Balanocystites, and Lagynocystis.)

# For the morphological dataset, there is no match. It is reconstructed as
# following:
morph.anc[[50]]$matrix_1$matrix[367, ]
which(morph.distances.GED.5[[50]]$distance_matrix[367, ] < 1.1)
# The closest morphological tip is the cinctan Asturicystis followed by cinctant
# Protocinctus (although with so many missing states, probably not a good
# comparitor.




## BUILD PHYLOGENETIC DIVERSITY CURVES #########################################

## Use to build diversity curve (using 'Total Diversity' of Foote 2000)
# X-FL, number of taxa w/ 1st and last appearance in interval (incl. singletons)
#  [including those that range through the entire interval but do not pass into
#  adjacent intervals],
# X-bL, number of taxa crossing bottom boundary & last appear in interval,
# X-Ft, number of taxa crossing top boundary and first appear in interval,
# X-bt, number of taxa crossing both bottom and top interval boundaries,
#
## Total D, total number of taxa in interval (X-FL + X-bL + X-Ft + X-bt)

# Switch between using 'mean.ranges' and 'ranges[[x]]' depending on what want
# diversity curve you wish to calculate and plot

# Genus-level diversity curve (using database occurrences above)
mids <- apply(ages[ ,9:10], 1, mean)
divs <- data.frame(interval = ages$interval_name, base = ages$max_ma,
                   top = ages$min_ma, midpt = mids, div = NA)
# Use if wish to restrict to tips (i.e., to not run as a phylogenetic lineage
# richness)
incl.nodes <- TRUE
inc.rows <- if(incl.nodes) 1:nrow(mean.ranges) else 1:Ntip(tree)
for(t in 1:nrow(divs)) {
  FL <- length(which(mean.ranges[inc.rows, 1] <= divs$base[t] &
                       mean.ranges[inc.rows, 2] >= divs$top[t]))
  bL <- length(which(mean.ranges[inc.rows, 1] > divs$base[t] &
                       mean.ranges[inc.rows, 2] < divs$base[t] &
                       mean.ranges[inc.rows, 2] >= divs$top[t]))
  Ft <- length(which(mean.ranges[inc.rows, 2] < divs$top[t] &
                       mean.ranges[inc.rows, 1] <= divs$base[t] &
                       mean.ranges[inc.rows, 1] > divs$top[t]))
  bt <- length(which(mean.ranges[inc.rows, 1] > divs$base[t] &
                       mean.ranges[inc.rows, 2] < divs$top[t]))
  divs$div[t] <- FL + bL + Ft + bt
}
divs
summary(divs$div)


# Plot diversity curve
geoscalePlot2(divs$midpt, divs$div, units = c("Epoch", "Period"), 
              tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
              cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), 
              data.lim = c(0, 310), ts.col = TRUE, label = "lineage richness", 
              timescale = ICS2020, type = "n", abbrev = "Period")
mtext(text = "binned richness", side = 3, cex = 1.25)
lines(divs$midpt, divs$div, lwd=3)


# David Bapst's simple paleotree way to calculate phylogenetic diversity (+ 1
# added for root). Use tips for LADs and nodes for FADs
LAD <- ntime[1:Ntip(tree)]
FAD <- ntime[(Ntip(tree) + 1):length(ntime)]
int.start <- int.times[, 1]
int.end <- int.times[, 2]
div <- sapply(1:length(int.start), function(x)
  1 + sum(FAD >= int.end[x]) - sum(LAD > int.start[x]))

# Compare binned and "continuous" richness trends
# Binned version using paleotree method
mids <- (int.start + int.end) / 2
geoscalePlot2(mids, div, units = c("Epoch", "Period"), tick.scale = "Period",
              boxes = "Age", cex.age = 0.65, cex.ts = 0.7, cex.pt = 1,
              age.lim = c(540, 445), data.lim = c(0, 310), ts.col = TRUE, 
              label = "lineage richness", timescale = ICS2020, type = "n", 
              abbrev = "Period")
mtext(text = "binned lineage richness", side = 3, cex = 1.25)
lines(mids, div, lwd=3)


# Continuous version
div.cont <- phyloDiv(tree = tree, plot = FALSE)
mids.cont <- (div.cont[ ,1] + div.cont[ ,2]) / 2
geoscalePlot2(mids.cont, div.cont[, 3], units = c("Epoch", "Period"), 
             tick.scale = "Period", boxes = "Age", cex.age = 0.65, cex.ts = 0.7, 
             cex.pt = 1, age.lim = c(540, 445), data.lim = c(0, 310), 
             ts.col = TRUE, label = "lineage richness", timescale = ICS2020, 
             type = "n", abbrev = "Period")
mtext(text = "continuous lineage richness", side = 3, cex = 1.25)
lines(mids.cont, div.cont[ ,3], lwd=3)





## CREATE TAXON - TIME BIN MATRIX ##############################################
# Used to sample taxa in each bin for disparity analyses below.
taxon.bins <- vector("list", length(mode.anc))
for(t in 1:length(mode.anc)) {
  tree <- mode.anc[[t]]$topper$tree
  t.taxon.bins <- matrix(FALSE, nrow = (ape::Ntip(tree) + ape::Nnode(tree)),
                         ncol = nrow(ages))
  colnames(t.taxon.bins) <- ages$interval_name
  rownames(t.taxon.bins) <- rownames(morph.anc[[t]]$matrix_1$matrix)
  for (i in 1:ncol(t.taxon.bins)) {
    wh.row <- which(ranges[[t]][, 1] >= ages$min_ma[i] &
                      ranges[[t]][, 2] <= ages$max_ma[i])
    t.taxon.bins[wh.row, i] <- TRUE
  }
  taxon.bins[[t]] <- t.taxon.bins
}

# View output
taxon.bins[[1]][360:370,]

# Save object
# save(taxon.bins, file = "taxon.bins")
# load("taxon.bins")

# Confirm matches the diversity curve
new.div <- apply(taxon.bins[[1]], 2, sum)
data.frame(divs$div, new.div)
identical(divs$div, as.vector(new.div)) # TRUE (if use 'ranges[[1]]' above when calc divs)





## ASSIGN TAXONOMIC CLASS AND SUBPHYLUM FROM TREE ##############################

# Logic: Move from tips to nodes, classifying a node to a taxonomic class and
# subphylum (only) if both descendants are in same class/subphylum, and leaving
# unclassified otherwise. Also appends genus names to ease in future handling,
# and more intuitive confirmation working as intended. Because the order of tip
# labels varies among trees, need to re-order each set so output in the same
# order as the input data sets (which are roughly alphabetical).

# Import ecological data set:
data <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", header = TRUE, stringsAsFactors = FALSE)

# Confirm list of genera in spreadsheet is same as in tree
setdiff(mode.anc[[1]]$topper$tree$tip.label, data$Genus)
setdiff(data$Genus, mode.anc[[1]]$topper$tree$tip.label)
# Anatifopsis is false-positive because replaced with two subgenera
# (Anatiferocystis & Guichenocarpos in the tree files. These require special
# handling below.

# Change names to match tree tip labels
data$Genus[which(data$Genus == "Anatifopsis")] <- c("Anatiferocystis", "Guichenocarpos")

# Two passes. Start at tips and move down list, then reverse, moving from last
# node to root. Note the order of the list is determined by the tree$tip.labels,
# which differs from that in the raw input data (e.g.,
# EchinoLHData_Mode_NAreformatted.csv).
taxon.list <- vector("list", length(mode.anc))
for(t in 1:length(taxon.list)) {
  tree <- mode.anc[[t]]$topper$tree
  t.class.list <- t.subphylum.list <- character(Ntip(tree) + Nnode(tree))
  genus.match <- match(tree$tip.label, data$Genus)
  genus.list <- data$Genus[genus.match]
  t.class.list[seq.int(Ntip(tree))] <- data$Class[genus.match]
  t.subphylum.list[seq.int(Ntip(tree))] <- data$Subphylum[genus.match]
  orig.classes <- t.class.list
  orig.subphyla <- t.subphylum.list
  # First pass along tips
  for (n in 1:Ntip(tree)) {
    wh.node <-
      tree$edge[which(tree$edge[, 2] == n), ][1]   # Find node for this tip
    tip.classes <- t.class.list[tree$edge[which(tree$edge[,1] == wh.node), 2]]
    if (identical(tip.classes[1], tip.classes[2]))
      t.class.list[wh.node] <- tip.classes[1]
    tip.subphyla <- t.subphylum.list[tree$edge[which(tree$edge[,1] == wh.node), 2]]
    if (identical(tip.subphyla[1], tip.subphyla[2]))
      t.subphylum.list[wh.node] <- tip.subphyla[1]
    }
  # Second pass moving nodes down to root (incl. some redundancy from above)
  start <- Ntip(tree) + Nnode(tree)
  end <- Ntip(tree) + 1               # End at root (= first node)
  for (n in start:end) {
    wh.node <-
      tree$edge[which(tree$edge[, 2] == n), ][1]   # Find node for this tip
    tip.classes <- t.class.list[tree$edge[which(tree$edge[,1] == wh.node), 2]]
    if (identical(tip.classes[1], tip.classes[2]))
      t.class.list[wh.node] <- tip.classes[1]
    tip.subphyla <- t.subphylum.list[tree$edge[which(tree$edge[,1] == wh.node), 2]]
    if (identical(tip.subphyla[1], tip.subphyla[2]))
      t.subphylum.list[wh.node] <- tip.subphyla[1]
  }
  # Replace class/subphylum = "" with class = "UNCERTAIN"
  t.class.list <-
    replace(t.class.list, which(t.class.list == ""), "UNCERTAIN")
  t.subphylum.list <-
    replace(t.subphylum.list, which(t.subphylum.list == ""), "UNCERTAIN")
  # Change back to original input data order
  tip.match <- match(data$Genus, tree$tip.label)
  reordered.classes <-
    c(t.class.list[tip.match], t.class.list[(Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))])
  reordered.subphyla <-
    c(t.subphylum.list[tip.match], t.subphylum.list[(Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))])
  genus.list <-
    c(data$Genus, as.character((Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))))
  taxon.list[[t]] <- cbind(subphylum = reordered.subphyla, 
                           class = reordered.classes, genus = genus.list)
}

# save(taxon.list, file = "taxon.list")

# View output to confirm worked as intended (If monophyletic, nodes should be
# adjacent and ne less than the tips)
head(taxon.list[[1]])
taxon.list[[1]][360:370, ]
taxon.list[[1]][which(taxon.list[[1]][, "class"] == "Helicoplacoidea"), ]
taxon.list[[1]][which(taxon.list[[1]][, "class"] == "Ctenocystoidea"), ]
taxon.list[[1]][which(taxon.list[[1]][, "subphylum"] == "Echinozoa"), ]
table(taxon.list[[1]][, "class"], taxon.list[[1]][, "subphylum"])
# Note all Subphylum Uncertain also in Class Uncertain, but not the other way
# around (as expected if tips are in different classes in same subphylum)

# If each class is monophyletic, the final tallies should be 2 * Ntip(class i) -
# 1. Let's see:
orig.classes <- replace(orig.classes, which(orig.classes == ""), "UNCERTAIN") 
sort((table(orig.classes) * 2 - 1) - table(taxon.list[[10]][, "class"]), 
     decreasing = FALSE)
# Results: Most are monophyletic, but rhombiferans and eocrinoids have
# substantial non-monophyly, with smaller amounts among stenuroids,
# diploporitans, somasteroids, edrioasteroids, paracrinoids, asteroids, and the
# known paraphyletic 'diploporitans'.

orig.subphyla <- replace(orig.subphyla, which(orig.subphyla == ""), "UNCERTAIN") 
sort((table(orig.subphyla) * 2 - 1) - table(taxon.list[[10]][, "subphylum"]), 
     decreasing = FALSE)
# Results: Blastozoans have decent amount of non-monophyly and a few crinozoans
# do, too. But asterozoans and echinozoans are monophyletic.

sort(table(orig.classes), decreasing = FALSE)
sort(table(taxon.list[[50]][, "class"]), decreasing = FALSE)

sort(table(orig.subphyla), decreasing = FALSE)
sort(table(taxon.list[[50]][, "subphylum"]), decreasing = FALSE)


## CLASS DIVERSITY TRENDS ######################################################
unique.classes <- sort(unique(taxon.list[[1]][, "class"]))
class.bins <- matrix(NA, nrow = length(unique.classes), ncol = nrow(ages))
colnames(class.bins) <- ages$interval_name
rownames(class.bins) <- unique.classes
for (c in 1:nrow(class.bins)) {
  wh.class <- which(taxon.list[[1]][, "class"] == rownames(class.bins)[c])
  if (length(wh.class) == 1L)
    class.bins[c, ] <- as.numeric(taxon.bins[[1]][wh.class, ])
  else
    class.bins[c, ] <- apply(taxon.bins[[1]][wh.class, ], 2, sum)
}
head(class.bins)

# Plot diversity curve
cl.div <- sapply(1:ncol(class.bins), function(x) sum(class.bins[, x] > 0))
geoscalePlot2(mids, cl.div, units = c("Epoch", "Period"), 
              tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
              cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), 
              data.lim = c(0, max(cl.div) + 1), ts.col = TRUE, label = "class richness", 
              timescale = ICS2020, type = "n", abbrev = "Period")
mtext(text = "class lineage richness", side = 3, cex = 1.25)
lines(mids, cl.div, lwd=3)

# Plot stacked plot version
# pdf(file = "Stacked class richness.pdf")
par(mar = c(5, 4, 2, 2))
# Put most diverse at bottom
class.order <- order(table(taxon.list[[1]][, "class"]), decreasing = TRUE)
stack.color <- rev(viridisLite::inferno(nrow(class.bins)))
# Use next if prefer non-adjacents
# set.seed(1234); stack.color <- sample(rev(viridisLite::inferno(nrow(class.bins))))
stackpoly(x = -mids, y = t(class.bins[class.order, ]), col = stack.color,
          xlab = "time", ylab = "number of genera", stack = TRUE, 
          xlim = c(-541, -444), main = "genus richness (by class)")
legend("topleft", legend = rev(unique.classes[class.order]), pch = 15, ncol = 3, 
       pt.cex = 1, cex = .55, col = rev(stack.color), bty = "n")
par(op)
# dev.off()







## CALCULATE PCoA ORDINATIONS AND MORPHO/ECOSPACES #############################

# Given high correlations among distance matrices, using GED.5. (The distance
# metric is Will's 1988 generalized Euclidean distance, which replaces missing
# character states using the mean pairwise dissimilarity for each pairwise
# comparison. Hierarchical character dependencies are weighted using Hopkins and
# St. John 2018 alpha = 0.5, which in which the primary character contributes
# weight according to the fraction of shared secondary [and tertiary, etc.]
# characters.


## Pick desired distance matrix
dist.matrix <- mode.distances.GED.5
# dist.matrix <- constant.distances.GED.5
# dist.matrix <- raw.distances.GED.5
# dist.matrix <- morph.distances.GED.5

# View sample
dist.matrix$distance_matrix[1:10, 1:4]

# Observe data structure of the data matrix:
hist(c(dist.matrix$distance_matrix))
summary(c(dist.matrix$distance_matrix))

# Any missing distances? (Should only occur in pre-trimmed raw treatment)
anyNA(dist.matrix$distance_matrix)

# Any non-diagonal zeros (potentially problematic for some ordination analyses)
length(dist.matrix$distance_matrix == 0) - nrow(dist.matrix$distance_matrix)

# Is it Euclidean?
ade4::is.euclid(as.dist(dist.matrix$distance_matrix))




## Perform a phylogenetic Principal Coordinates Analysis:

# Skip 'MorphMatrix2PCoA' [which bundles ancestral state reconstruction and
# calculation of distance matrices] by using direct calculation of PCoA so can
# directly import in previously calculated distance matrices.

# For consistency with other Claddis functions, using ape::pcoa() instead of
# 'ecospace' and 'FD's use of ade4::dudi.pco(), which produces identical output
# (in different formats), and are essentially equally fast. Variation typically
# has to do with how many non-negative eigenvectors are returned, and with
# arbitrary reversing of axes. Eigenvectors are identical for at least first
# 11-15 axes. Using Claddis' recommended "Cailliez" correction here and below
# for consistency in how all morphospace/ecospace disparity metrics are
# calculated. Tests (not demonstrated here, see "X_Diff pcoa trials.R")
# demonstrate that the "Lingoes" correction works less often with the data
# structure of some data treatments, especially the high proportion of zero
# lengths for many data pairs. In the case of the constant treatment, a Cailliez
# correction is not possible, so "none" is used for this single treatment.

# Want to remove missing taxa? (Only used for 'raw' data treatment)
# trim.matrix <- TrimMorphDistMatrix(dist.matrix$distance_matrix, Tree = tree)
# dist.matrix$distance_matrix <- trim.matrix$DistMatrix
# dist.matrix$Tree <- trim.matrix$Tree
# dist.matrix$RemovedTaxa <- trim.matrix$RemovedTaxa

pcoa_results <- ape::pcoa(dist.matrix$distance_matrix,
                          correction = "cailliez", 
                          rn = rownames(dist.matrix$distance_matrix))
# Append Tree, required for Claddis function
if(is.null(pcoa_results$Tree)) pcoa_results$Tree <- tree
# Use following ONLY for trimmed 'raw' treatment:
# if(is.null(pcoa_results$Tree)) pcoa_results$Tree <- trim.matrix$Tree
beep()

# Save PCoA output
# mode.pcoa <- pcoa_results; save(mode.pcoa, file = "mode.pcoa")
# constant.pcoa <- pcoa_results; save(constant.pcoa, file = "constant.pcoa")
# raw.pcoa <- pcoa_results; save(raw.pcoa, file = "raw.pcoa")
# morph.pcoa <- pcoa_results; save(morph.pcoa, file = "morph.pcoa")
# load("morph.pcoa"); pcoa_results <- morph.pcoa
# load("mode.pcoa"); pcoa_results <- mode.pcoa
# load("constant.pcoa"); pcoa_results <- constant.pcoa
# load("raw.pcoa"); pcoa_results <- raw.pcoa


# Plot scree plot and biplot
barplot(100 * pcoa_results$values$Rel_corr_eig[1:15], main = "first 10 axes", names.arg = 1:15,
        xlab = "Axes", ylab = "relative % explained")
barplot(100 * cumsum(pcoa_results$values$Rel_corr_eig)[1:15], names.arg = 1:15,
        xlab = "Axes", ylab = "relative % explained", 
        main = "cum var explained for first 10 axes")
plot(pcoa_results$vectors.cor[, 1], pcoa_results$vectors.cor[, 2], xlab = "PCoA axis 1", 
     ylab = "PCoA axis 2")

round(100 * pcoa_results$values$Rel_corr_eig[1:10], 2)
round(100 * cumsum(pcoa_results$values$Rel_corr_eig)[1:30], 1)


# RESULTS (% explained; cumulative % explained):

# Morph:    1- 1.50%, 2-0.94% (2.4%), 3-0.66% (3.1%), 4-0.42% (3.5%), 
#                     5-0.32% (3.8%), 6-0.27% (4.1%),   
# Mode:     1- 1.54%, 2-0.62% (2.2%), 3-0.47% (2.6%), 4-0.37% (3.0%), 
#                     5-0.32% (3.3%), 6-0.27% (3.6%),   
# Constant: 1- 0.81%, 2-0.29% (1.1%), 3-0.29% (1.3%), 4-0.24% (1.5%), 
#                     5-0.19% (1.7%), 6-0.17% (1.9%),   
# Raw:      1- 0.63%, 2-0.43% (1.1%), 3-0.33% (1.4%), 4-0.30% (1.7%), 
#                     5-0.26% (1.9%), 6-0.25% (2.2%),   

# Plot all scree plots on one figure
# pdf(file = "PCoAScreePlots.pdf")
par(mfrow = c(2,2), mar = c(4, 4, 2, .1))
load("morph.pcoa")
barplot(100 * morph.pcoa$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "morphology", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
load("mode.pcoa")
barplot(100 * mode.pcoa$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (mode)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
load("constant.pcoa")
# note no correction used, but pcoa() still calculates 'corrected' eigenvalues
barplot(100 * constant.pcoa$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (constant)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
load("raw.pcoa")
barplot(100 * raw.pcoa$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (raw)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
par(op)
# dev.off()



## Calculate (mock) "factor loadings"
# Because PCoA uses the distance matrix instead of observed variables, it is not
# possible to relate the original variables to the PCoA axes. However, the basic
# premise of a correlation between PCoA space and original variables can be
# informative.
#   orig.vars = data with original variables
#   ord.coord = principal coordinate vectors output from ape::pcoa
#   vars      = how many principal coordinate axes do you want?
#   cutoff    = what absolute value do you want to mask uncorrelated variables? 
mock.loadings <- function(orig.vars = NULL, ord.coord = NULL, vars = 6, 
                          cutoff = 0.5) {
  out <- matrix(NA, nrow = ncol(orig.vars), ncol = vars)
  for (r in 1:nrow(out)) {
    # exclude any characters that are all NAs
    for (c in 1:vars) {
      if (all(is.na(as.numeric(orig.vars[, r]))))
        next
      # exclude any characters in which only two two are coded (b/c r is
      # uninformative because must equal 1)
      if (sum(!is.na(as.numeric(orig.vars[, r]))) < 3L)
        next
      out[r, c] <- cor(as.numeric(orig.vars[, r]), ord.coord[, c],
                       use = "complete.obs")
    }
  }
  out <- round(out, 3)
  out <- as.data.frame(replace(out, which(abs(out) < cutoff), "-"))
  return(out)
}

loadings.morph <- mock.loadings(orig.vars = morph.anc$matrix_1$matrix,
                                ord.coord = morph.pcoa$vectors.cor, vars = 6,
                                cutoff = 0.8)
na.omit(loadings.morph)

loadings.mode <- mock.loadings(orig.vars = mode.anc$matrix_1$matrix, 
                               ord.coord = mode.pcoa$vectors.cor, vars = 6, 
                               cutoff = 0.4)
na.omit(loadings.mode)

# Note using uncorrected eigenvectors for 'constant' treatment:
loadings.constant <- mock.loadings(orig.vars = constant.anc$matrix_1$matrix,
                                   ord.coord = constant.pcoa$vectors, vars = 6,
                                   cutoff = 0.4)
na.omit(loadings.constant)

# For the raw treatment, need to remove genera with missing data excluded from
# PCoA (and just focusing on tips to keep things simpler because nodes were
# renumbered when making distance matrix)
nt <- Ntip(raw.pcoa$Tree)
wh.to.keep <- na.omit(match(rownames(raw.pcoa$vectors.cor[1:nt, ]),
                            rownames(raw.anc$matrix_1$matrix)[1:nt]))
loadings.raw <- mock.loadings(orig.vars = raw.anc$matrix_1$matrix[wh.to.keep, ],
                              ord.coord = raw.pcoa$vectors.cor[1:length(wh.to.keep), ], vars = 6, 
                              cutoff = 0.3)
na.omit(loadings.raw)

# RESULTS:
# MORPHOLOGY DATA SET:     -                               +
#  PCO 1:     32, 135, 210, 227, 325, 411               4, 134, 141, 409
#  PCO 2:    160                                       14, 46, 126, 132
#  PCO 3:     21, 121, 126, 134, 141, 146, 227, 380    46
#  PCO 4:     46, 111, 121, 172
#  PCO 5:    131, 135                                 376
#  PCO 6:     31, 129, 130, 201                        46, 140, 143

# MODE DATA SET:      -                                +
#  PCO 1:     2-5, (12), 15, (28-29)                    6, (13), 16, (31)                 
#   +: Attached (filter feeders) living far from (hard) substrate
#   -: Free-living, mobile (and mass-feeders living on soft substrates)
#  PCO 2:                                              (9, 29)
#   +: High-density filters on biotic substrates
#   -: everything else?
#  PCO 3:     1 (24, 26, 31, 36, 40)                   (9, 23, 25, 28)
#   +: Smaller? (Epifaunal filter feeders on biotic substrates)
#   -: Large-bodied (and carnivorous bulk feeders with infaunal food)
#  PCO 4:     (9, 24, 26, 31, 32, 36, 40)             (23, 25, 28)
#   +: (Epifaunal filter feeders)
#   -: (Scavengers on biotic substrates with infaunal food) 
#  PCO 5:     (10, 24, 26, 36, 40)                    (23, 25)
#   +: (Food above substrate)
#   -: (Carnivorous scavengers on lithic substrates with infaunal food)
#  PCO 6:                                              10
#   +: Lithic substrates
#   -: (Biotic substrates)

# CONSTANT DATA SET:  -                                +
#  PCO 1:        6, (13), 16, (31)                  2-5, 12, 15, (28-29)
#  PCO 2:      (23, 25, 28, 39)                     1 (24, 26, 31, 36, 40)
#  PCO 3:      (29)                                        
#  PCO 4:      (13, 23, 25, 28, 39)               (12, 24, 26, 36, 40)
#  PCO 5:       (9)
#  PCO 6:       (6)                               (13, 29)

#    OVERALL: PCOs basically same as for mode, with 1 reversed, 2 = mode's 3, 3
#    = mode's 2, and 4 reversed. 5 and 6 slightly different, but low loading
#    scores.

# RAW DATA SET:       -                                +
#  PCO 1:            1-5, 12, 15                  6, 13, 16
#  PCO 2:           12, 15                        1, 6, 13, 16
#  PCO 3:                                        (2-4, 29)
#  PCO 4:                                        (1)
#  PCO 5:          (28)                         (31)                      
#  PCO 6:                                         1

#    OVERALL: PCOs basically same as for mode. 1 identical, 2 = mode's 1, 3 =
#    mode's 2, 3-6 combination of axes.


# Effect of applying the Dineen, et al. (2019, Biology Lettters) standardization
# for (later) PCoA convex hull volume by re-scaling PCoA coordinates according
# to the magnitude of the first eigenvalue. (Their eq. 1 in Supplementary
# Information.)
stand.pcoa <- function(vectors = NULL, eigenvalues = NULL) {
  nc <- ncol(vectors)
  first.e <- eigenvalues[1]
  for (c in 1:nc) {
    vectors[, c] <- vectors[, c] * eigenvalues[c] / first.e
  }
  return(vectors)
}

# Compare before and after standardizing
new <- pcoa_results
new$vectors.cor <- stand.pcoa(vectors = pcoa_results$vectors.cor,
                          eigenvalues = pcoa_results$values$Corr_eig)
pcoa_results$vectors.cor[1:5, 1:5]
new$vectors.cor[1:5, 1:5]
plot(pcoa_results$vectors.cor[, 1], pcoa_results$vectors.cor[, 2])
plot(new$vectors.cor[, 1], new$vectors.cor[, 2])

# CONCLUSION: The axes are rescaled correctly. (Note that doing so affects
# [i.e., stabilizes across time intervals] the values returned for FRic
# [especially so], FEve [minimally so], and FDiv [minimally so].)








## PLOT ORDINATIONS USING CLADDIS FUNCTIONS ####################################

# Plot a 2-dimensional morphospace with tips, nodes, and root:
MorphospacePlot(pcoa_results)
tip.seq <- 1:Ntip(pcoa_results$Tree)
node.seq <-
  (Ntip(pcoa_results$Tree) + 1):(Ntip(pcoa_results$Tree) + Nnode(pcoa_results$Tree))
points(x = pcoa_results$vectors.cor[node.seq, 1],
       y = pcoa_results$vectors.cor[node.seq, 2], col = "gray", pch = 16) # nodes
points(x = pcoa_results$vectors.cor[tip.seq, 1], 
       y = pcoa_results$vectors.cor[tip.seq, 2])                          # tips
points(x = pcoa_results$vectors.cor[max(node.seq), 1], 
       y = pcoa_results$vectors.cor[max(node.seq), 2], col = "red", pch = 16, 
       cex = 1.5)                                                     # root

# Plot additional morphospace axes:
MultiMorphospacePlot(pcoa_results, N_axes = 3)

# Plot a chronophylomorphospace:
ChronoPhyloMorphospacePlot(pcoa_results)



## Plot stacked (Footeogram) ordination spaces
# Note that this uses raw strat ranges of tips and not the ancestral nodes

# 1. Import strat ranges (previously downloaded fom PBDB)
strat.ranges <- read.csv("~/Manuscripts/CamOrdEchinos/Data files/GenusStratRanges.csv", 
                 header = TRUE)
strat.ranges <- strat.ranges[, -1]
colnames(strat.ranges) <- c("FAD", "LAD")

# Use series ages for bins
strat_names <- read.csv("~/Manuscripts/CamOrdEchinos/strat_names.csv")
ages2 <- strat_names[which(strat_names$scale_level == 4), ]
ages2 <- ages2[which(ages2$max_ma > 443), ]

StackPlot(ordination_axes = pcoa_results$vectors.cor, ages = strat.ranges, 
          time_slices = ages2$max_ma)

# pdf(file = "FootePlot_morph.pdf")
pcoa_results <- morph.pcoa
StackPlot(ordination_axes = pcoa_results$vectors.cor, ages = strat.ranges, 
          time_slices = ages2$max_ma, axis_label = "PCoA")
mtext(text = "morphology, tips only", side = 3, cex = 1.25)
# dev.off()

# pdf(file = "FootePlot_mode.pdf")
pcoa_results <- mode.pcoa
StackPlot(ordination_axes = pcoa_results$vectors.cor, ages = strat.ranges, 
          time_slices = ages2$max_ma, axis_label = "PCoA")
mtext(text = "ecology (mode), tips only", side = 3, cex = 1.25)
# dev.off()

# pdf(file = "FootePlot_constant.pdf")
# Using uncorrected eigenvectors here
pcoa_results <- constant.pcoa
StackPlot(ordination_axes = pcoa_results$vectors, ages = strat.ranges, 
          time_slices = ages2$max_ma, axis_label = "PCoA")
mtext(text = "ecology (constant), tips only", side = 3, cex = 1.25)
# dev.off()

# pdf(file = "FootePlot_raw.pdf")
pcoa_results <- raw.pcoa
StackPlot(ordination_axes = pcoa_results$vectors.cor, ages = strat.ranges, 
          time_slices = ages2$max_ma, axis_label = "PCoA")
mtext(text = "ecology (raw), tips only", side = 3, cex = 1.25)
# dev.off()





## DISPARITY / FUNCTIONAL DIVERSITY TRENDS #####################################

# Load PCoA output from ape::pcoa
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")
load("morph.pcoa")

# Choose ancestral states, Wills GED-0.5 distance matrices, and PCoA output
anc <- mode.anc$matrix_1$matrix; dist.matrix <- mode.distances.GED.5$distance_matrix; pcoa_results <- mode.pcoa
# anc <- constant.anc$matrix_1$matrix; dist.matrix <- constant.distances.GED.5$distance_matrix; pcoa_results <- constant.pcoa
# anc <- raw.anc$matrix_1$matrix; dist.matrix <- raw.distances.GED.5$distance_matrix; pcoa_results <- raw.pcoa
# anc <- morph.anc$matrix_1$matrix; dist.matrix <- morph.distances.GED.5$distance_matrix; pcoa_results <- morph.pcoa

# Load taxon.bins (built above)
load("taxon.bins")
if(nrow(taxon.bins) != nrow(dist.matrix))
  stop("Need to re-load or re-build taxon.bins because some taxa were removed when running the 'raw' treatment.\n")

# Because a correction wasn't able to be applied for the 'constant' treatment,
# need to amend as if it were (for consistency in downstream scripts)
# pcoa_results$vectors.cor <- pcoa_results$vectors
# pcoa_results$values$Corr_eig <- pcoa_results$values$Eigenvalues

## Want to remove missing taxa? (Only used for 'raw' data treatment)
# Note the 'tree' is not included to prevent renumbering! For some reason I
# can't figure out, there is a discrepancy in which taxa are removed when the
# tree is present, with 515 taxa with tree and 511 without. In order for the
# trends code to work, the following is required to get all associated disparity
# objects with same taxa.
# trim.matrix <- TrimMorphDistMatrix(dist.matrix)
# trim.matrix.for.pcoa <- TrimMorphDistMatrix(dist.matrix, Tree = tree)
# dim(trim.matrix.for.pcoa$DistMatrix)
# dim(trim.matrix$DistMatrix)

# first pass:
# wh.to.keep <- match(rownames(trim.matrix$DistMatrix), rownames(anc))
# anc <- anc[wh.to.keep, ]
# taxon.bins <- taxon.bins[wh.to.keep, ]
# dist.matrix <- trim.matrix$DistMatrix
# Note still wrong dims:
# dim(anc); dim(taxon.bins); dim(dist.matrix); dim(pcoa_results$vectors.cor)

# Easy to fix with named tips:
# wh.to.cut <- which(!rownames(dist.matrix)[1:258] %in% rownames(pcoa_results$vectors.cor)[1:256])
# rownames(dist.matrix)[wh.to.cut] # row 202: Protocrinites & row 228: Sphaeronites

# Not so easy with re-named ancestors, so need a work-around: Go through and
# manually find the non-kept row based on dissimilarities :(
# (Ignore error for now)
# cbind(trim.matrix$DistMatrix[-wh.to.cut, 1], trim.matrix.for.pcoa$DistMatrix[, 1])

# Turns out two additional ancestors were removed (row 388: ancestor 530 and
# row: 394: ancestor 536). (Technically, ancestor 530 or 531 was cut, as they
# have identical distances to all other taxa, but keeping 531 as it has the
# longer stratigraphic range.)
# wh.to.cut <- c(202, 228, 388, 394)
# anc <- anc[-wh.to.cut, ]
# taxon.bins <- taxon.bins[-wh.to.cut, ]
# dist.matrix <- dist.matrix[-wh.to.cut, -wh.to.cut]
# dim(anc); dim(taxon.bins); dim(dist.matrix); dim(pcoa_results$vectors.cor)
# Yay, all now match! Use these for disparity analyses.


# View to confirm 
taxon.bins[1:5, 1:5]
anc[1:10, 1:10]
dist.matrix[1:4, 1:4]
pcoa_results$vectors.cor[1:4, 1:4]

## How many dimensions in PCoA for FRic and FDiv?

# The optimal choice (e.g., Villeger, et al., 2011; Lloyd,2018) is the highest
# value that is computationally tractable; MAX possible = H - 1, but anything
# greater than 8 is much too slow (especially when use resampling below to
# assess variability). Using m = 6 because all PCoA results show that the first
# 6 PC axes account for at least 80% of explained variance.)

# Timing (and qual.FRic) using largest 'mode' and morph' data sets:
#             m = 3:    30 s (82%)      26 s (78%)
#             m = 4:    25 s (88%)      24 s (87%)
#             m = 5:    22 s (92%)      33 s (90%)
#             m = 6:    22 s (95%)      25 s (92%) *** best trade-off ***
#             m = 7:    27 s (97%)      36 s (95%)
#             m = 8:    46 s (99%)      79 s (97%)
#             m = 9:   164 s (~100%)   337 s (98%)


# See Novack-Gottshall (2015, 2016) for details of each statistic
#  S  =  phylogenetic lineage richness (using 'anc')
#  H  =  number of unique life habits or morphotypes (using 'dist.matrix')
#  D  =  mean pairwise Wills generalized Euclidean distance (GED) (using
#        'dist.matrix')
#  M  =  maximum range in GED (using 'dist.matrix')
#  V  =  total variance (using 'anc')
# FRic = functional richness / convex hull volume (using 'pcoa' eigenvectors)
# FEve = functional evenness (evenness in minimum spanning tree) (using
#        'pcoa' eigenvectors)
# FDiv = functional divergence (mean distance of genera from PCoA centroid)
#        (using 'pcoa' eigenvectors)
# FDis = functional dispersion (total deviance from circle with radius equal to
#        mean distance from PCoA centroid) (using 'dist.matrix')

# Note that because of problematic taxa, unable to calculate FRic and FDiv with
# the 'raw' treatment. So set calc.FRic.and.FDiv = FALSE
m <- 6
(start <- Sys.time())
nc <- ncol(taxon.bins)
metrics <- data.frame(Age = as.numeric(mids), S = NA, H = NA, D = NA, M = NA, 
                      V = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, 
                      qual.FRic = NA)
for(t in 1:nc) {
  wh.gr <- unname(which(taxon.bins[, t]))
  S <- length(wh.gr)
  cat("age =", mids[t], ", ", S, "genera", "\n")
  anc.sam <- anc[wh.gr, ]
  dist.anc.sam <- dist.matrix[wh.gr, wh.gr]
  pcoa.sam <- pcoa_results
  pcoa.sam$vectors.cor <- pcoa_results$vectors.cor[wh.gr, ]
  H <- nrow(unique(dist.anc.sam)) # NEW CHANGE: BETTER NA HANDLING
  if (any(is.nan(dist.anc.sam)) | length(dist.anc.sam) == 0) next
  if (S <= m | H <= m) next
  FD <- calc_metrics2(sample = anc.sam, dist.sam = dist.anc.sam, 
                      pcoa = pcoa.sam, m = m, stand.pcoa = TRUE, 
                      calc.FRic.and.FDiv = TRUE)
  metrics[t, 1 + seq.int(length(FD))] <- FD
}
(Sys.time() - start)

# What proportion of PCoA eigenvalues included in ordination-based metrics?
# (Note that because calculated on the entire PCoA eigenvectors, will not change
# when dealing with subsets.)
round(summary(metrics$qual.FRic), 3)

beep(3)
head(metrics)


## Save metrics
# write.csv(metrics, file = "metrics_morph.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_LH_mode.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_LH_constant.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_LH_raw.csv", row.names = FALSE)
# metrics <- read.csv(file = "metrics_morph.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_LH_mode.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_LH_constant.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_LH_raw.csv", header = TRUE)

par(mfrow = c(2, 4), mar = c(4, 4, 1, 0.2))
for (c in 3:10) {
  if (sum(is.na(metrics[, c])) == length(metrics[, c]))
    plot( 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  else
    plot(metrics$S, metrics[, c], main = colnames(metrics)[c], 
         ylab = colnames(metrics)[c], xlab = "S")
}
par(op)


# Plot trends
par(mar = c(0, 4, 2, 2))
for (c in 2:10) {
  var <- metrics[, c]
  if (sum(is.na(var)) == length(var))
    next
  lim <- range(var, na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
               tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
               cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
               ts.col = TRUE, label = colnames(metrics)[c], timescale = ICS2020, 
               type = "n", abbrev = "Period")
  mtext(text = colnames(metrics)[c], side = 3, cex = 1.25)
  lines(mids, var, lwd = 3)
}
par(op)



## How correlated are the metrics?
morph <- read.csv(file = "metrics_morph.csv", header = TRUE)
mode <- read.csv(file = "metrics_LH_mode.csv", header = TRUE)
constant <- read.csv(file = "metrics_LH_constant.csv", header = TRUE)
raw <- read.csv(file = "metrics_LH_raw.csv", header = TRUE)

# Calculate first-differences (change in slope across time):
diff.morph <- as.data.frame(apply(morph, 2, diff))
diff.morph <- diff.morph / diff.morph$Age

diff.mode <- as.data.frame(apply(mode, 2, diff))
diff.mode <- diff.mode / diff.mode$Age

diff.constant <- as.data.frame(apply(constant, 2, diff))
diff.constant <- diff.constant / diff.constant$Age

diff.raw <- as.data.frame(apply(raw, 2, diff))
diff.raw <- diff.raw / diff.raw$Age

# morphological vs. LH-mode
round(diag(cor(diff.morph[, 3:10], diff.mode[, 3:10], use = "pairwise.complete.obs")), 4)
summary(lm(diff.morph$H ~ diff.mode$H))

# morphological vs. LH-constant
round(diag(cor(diff.morph[, 3:10], diff.constant[, 3:10], use = "pairwise.complete.obs")), 4)
summary(lm(diff.morph$H ~ diff.constant$H))

# morphological vs. LH-raw
round(diag(cor(diff.morph[, 3:10], diff.raw[, 3:10], use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-constant
round(diag(cor(diff.mode[, 3:10], diff.constant[, 3:10], use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-raw
round(diag(cor(diff.mode[, 3:10], diff.raw[, 3:10], use = "pairwise.complete.obs")), 4)

# LH-constant vs. LH-raw
round(diag(cor(diff.constant[, 3:10], diff.raw[, 3:10], use = "pairwise.complete.obs")), 4)





# Plot raw H for morphology and mode
lim <- c(0, max(c(morph$H, mode$H), na.rm = TRUE))
geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
              tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
              cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
              ts.col = TRUE, label = "#", timescale = ICS2020, 
              type = "n", abbrev = "Period")
mtext(text = "# morphotypes / life habits", side = 3, cex = 1.25)
cols <- viridisLite::plasma(3)         # Uses RGB-sensitive red and blue
lines(mids, mode$H, lwd = 4, lty = 1, col = cols[1])
lines(mids, morph$H, lwd = 4, lty = 6, col = cols[2])
legend("topleft", legend = c("ecology", "morphology"), col = cols, 
       bty = "n", lty = c(1, 6), lwd = 4, inset = 0.02, cex = 1.5)





## DISPARITY / FD TRENDS (AT STANDARD SAMPLE SIZE) #############################

# corrected resampling function
sample2 <- function(x, ...) x[sample.int(length(x), ...)]

# Load PCoA output from ape::pcoa
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")
load("morph.pcoa")

# Import ancestral states from 2-InferancestralStates.R
load("mode.anc")
load("constant.anc")
load("raw.anc")
load("morph.anc")

# Import distance matrices from 3-DisparityDistances.R
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")

# Choose ancestral states, Wills GED-0.5 distance matrices, and PCoA output
anc <- mode.anc$matrix_1$matrix; dist.matrix <- mode.distances.GED.5$distance_matrix; pcoa_results <- mode.pcoa
# anc <- constant.anc$matrix_1$matrix; dist.matrix <- constant.distances.GED.5$distance_matrix; pcoa_results <- constant.pcoa
# anc <- raw.anc$matrix_1$matrix; dist.matrix <- raw.distances.GED.5$distance_matrix; pcoa_results <- raw.pcoa
# anc <- morph.anc$matrix_1$matrix; dist.matrix <- morph.distances.GED.5$distance_matrix; pcoa_results <- morph.pcoa

# Load taxon.bins and taxon.list (built above)
load("taxon.bins")
load("taxon.list")
if(nrow(taxon.bins) != nrow(dist.matrix))
  stop("Need to re-load or re-build taxon.bins because some taxa were removed when running the 'raw' treatment.\n")
if(nrow(taxon.bins) != length(taxon.list))
  stop("Need to re-build taxon.list because some taxa were removed when running the 'raw' treatment.\n")

# Because a correction wasn't able to be applied for the 'constant' treatment,
# need to amend as if it were (for consistency in downstream scripts)
# pcoa_results$vectors.cor <- pcoa_results$vectors
# pcoa_results$values$Corr_eig <- pcoa_results$values$Eigenvalues

## Want to remove missing taxa? (Only used for 'raw' data treatment)
# Note the 'tree' is not included to prevent renumbering! See above for more
# explanation.
# trim.matrix <- TrimMorphDistMatrix(dist.matrix)
# wh.to.keep <- match(rownames(trim.matrix$DistMatrix), rownames(anc))
# anc <- anc[wh.to.keep, ]
# taxon.bins <- taxon.bins[wh.to.keep, ]
# dist.matrix <- trim.matrix$DistMatrix
# dim(anc); dim(taxon.bins); dim(dist.matrix); dim(pcoa_results$vectors.cor)
# wh.to.cut <- c(202, 228, 388, 394)
# anc <- anc[-wh.to.cut, ]
# taxon.bins <- taxon.bins[-wh.to.cut, ]
# dist.matrix <- dist.matrix[-wh.to.cut, -wh.to.cut]
# Confirm all these 4 objects have the same number of rows!
# dim(anc); dim(taxon.bins); dim(dist.matrix); dim(pcoa_results$vectors.cor)

# View to confirm
anc[1:5, 1:10]
dist.matrix[1:5, 1:4]
taxon.bins[1:5, 1:4]
pcoa_results$vectors.cor[1:5, 1:5]

## Number of genera per interval

# Use 50 for echinoderm-wide
std.g <- 50
# Use TRUE if want to first check sample sizes to obtain std.g
check.std.g <- FALSE

# Initialize parallel cluster:
# See numreps_confirmation.xlsx for confirmation for why 3,000 replicates is
# sufficient to get relative error (CV = sd/mean) within 0.1% for all time
# intervals and (phylum-wide) sample sizes. Tested using the two time intervals
# with largest (= Floian) and smallest (= Paibian) species richness, using the
# fast-to-calculate D metric. (As expected, intervals with greater richness have
# greater sampling variability, which decreases with smaller resampling pools
# and greater numbers of replicates.)
numrep <- 3000
cl <- makeCluster(detectCores())
registerDoParallel(cl)
# Use load-balancing because of different run times for the MCMC optimizations
opts <- list(preschedule = FALSE)
clusterSetRNGStream(cl, 3142)  # Set L-Ecuyer RNG seed

# How many dimensions in PCoA for FRic and FDiv? Use 6 for phylum-level
m <- 6
nc <- ncol(taxon.bins)
metrics <- data.frame(Age = as.numeric(mids), S = std.g, H = NA, SE.H = NA,
                      D = NA, SE.D = NA, M = NA, SE.M = NA, V = NA, SE.V = NA, FRic = NA, 
                      SE.FRic = NA, FEve = NA, SE.FEve = NA, FDiv = NA, SE.FDiv = NA, 
                      FDis = NA, SE.FDis = NA, qual.FRic = NA, SE.qual.FRic = NA)
(start <- Sys.time())
# Loop through each time interval, running the sample-standardization replicates
# in parallel

# (Was done in following manner in summer 2020, but if need to run again, refer
# to 6-Phylogenetic inertia.R for better way to implement load-balancing where
# the L-Ecuyer RNG seed can allow replication.)

for(t in 1:nc) {
  wh.gr <- unname(which(taxon.bins[, t]))
  S <- length(wh.gr)
  cat(round(Sys.time() - start, 0), "min., age =", mids[t], ", ", 
      S, "genera", "\n")
  if (check.std.g) next
  if(S < std.g) metrics[t, 3:20] <- NA
  if(S < std.g) next
  anc.sam <- anc[wh.gr, ]
  dist.anc.sam <- dist.matrix[wh.gr, wh.gr]
  pcoa.sam <- pcoa_results
  pcoa.sam$vectors.cor <- pcoa.sam$vectors.cor[wh.gr, ]
  H <- nrow(unique(dist.anc.sam))
  if (any(is.nan(dist.anc.sam)) | length(dist.anc.sam) == 0) next
  if (S <= m | H <= m) next
  max.seq <- seq.int(nrow(anc.sam))
  
  # Build numrep in parallel
  bin.metrics <- foreach(i = 1:numrep, .options.snow = opts, .combine = "rbind", 
                         .inorder = FALSE, .packages = "FD") %dopar% {
                           
    # Create a new subsample each replicate
    sampled <- sample2(max.seq, std.g, replace = FALSE)
    sub.sam <- anc.sam[sampled, ]
    sub.dist <- dist.anc.sam[sampled, sampled]
    sub.pcoa <- pcoa.sam
    sub.pcoa$vectors.cor <- sub.pcoa$vectors.cor[sampled, ]

    # Calculate metrics for each subsample
    FD <- calc_metrics2(sample = sub.sam, dist.sam = sub.dist, 
                        pcoa = sub.pcoa, m = m, stand.pcoa = TRUE, 
                        calc.FRic.and.FDiv = FALSE)
    }

  # Summarize bin metrics
  mean.cols <- c(2, 3, 5, 7, 9, 11, 13, 15, 17, 19)
  sd.cols <- c(4, 6, 8, 10, 12, 14, 16, 18, 20)
  metrics[t, mean.cols] <- apply(bin.metrics, 2, mean)
  metrics[t, sd.cols] <- apply(bin.metrics[,-1], 2, sd)
}
stopCluster(cl)
(Sys.time() - start)
beep(3)
# 28 min for 8 cores, 3000 replicates using mode treatment
# 29 min for 8 cores, 3000 replicates using constant treatment
# 28 min for 8 cores, 3000 replicates using raw treatment (b/c no FRic & DViv)
# 33 min for 8 cores, 3000 replicates using morph treatment

## Save metrics
# write.csv(metrics, file = "metrics_StdG50_morph.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_mode.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_constant.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_raw.csv", row.names = FALSE)
# metrics <- read.csv(file = "metrics_StdG50_morph.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_mode.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_constant.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_raw.csv", header = TRUE)
head(metrics, 3)

# With m = 6 PCoA axes used for calculation of FRic and FDiv, what proportion of
# the variability is explained in different time bins? (I.e., how much reduction
# has occurred by the ordination?)
summary(metrics$qual.FRic)
# Morphology: 4.1% used
# Ecology: mode: 3.6% used, constant: 24-45% used, raw: 2.2% used (but irrelevant)

# Plot trends
means <- c(3, 5, 7, 9, 11, 13, 15, 17) # Columns with mean values
SEs <- c(4, 6, 8, 10, 12, 14, 16, 18)  # Columns with SE values
par(mar = c(0, 4, 4, 2))
for (c in 1:8) {
  var.means <- metrics[ ,means[c]]
  var.SEs <- metrics[ ,SEs[c]]
  if (sum(is.na(var.means)) == length(var.means))
    next
  lim <- range(c(var.means - var.SEs, var.means + var.SEs), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
                tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
                cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), 
                data.lim = lim, ts.col = TRUE, 
                label = colnames(metrics)[means[c]], 
                timescale = ICS2020, type = "n", abbrev = "Period")
  mtext(text = paste(colnames(metrics)[means[c]], " (G = ", std.g, ")", sep = ""), 
        side = 3, cex = 1.25)
  lines(mids, var.means, lwd = 3)
  lines(mids, var.means + var.SEs, lwd = 2, lty = 2)
  lines(mids, var.means - var.SEs, lwd = 2, lty = 2)
}
par(op)









## How correlated are the standardized metrics?
morph <- read.csv(file = "metrics_StdG50_morph.csv", header = TRUE)
mode <- read.csv(file = "metrics_StdG50_LH_mode.csv", header = TRUE)
constant <- read.csv(file = "metrics_StdG50_LH_constant.csv", header = TRUE)
raw <- read.csv(file = "metrics_StdG50_LH_raw.csv", header = TRUE)

# Which columns have statistical means?
means <- c(3, 5, 7, 9, 11, 13, 15, 17)

# Calculate first-differences (change in slope across time):
diff.morph <- as.data.frame(apply(morph, 2, diff))
diff.morph <- diff.morph / diff.morph$Age

diff.mode <- as.data.frame(apply(mode, 2, diff))
diff.mode <- diff.mode / diff.mode$Age

diff.constant <- as.data.frame(apply(constant, 2, diff))
diff.constant <- diff.constant / diff.constant$Age

diff.raw <- as.data.frame(apply(raw, 2, diff))
diff.raw <- diff.raw / diff.raw$Age


# Because there are no duplicated morphotypes, H = 50 for all intervals. Ignore
# "H" (and warning) in comparisons with the morphology data set.

# morphological vs. LH-mode: not very correlated (D highest)
round(diag(cor(diff.morph[, means], diff.mode[, means], 
               use = "pairwise.complete.obs")), 4)

# Are the D and FRic statistically correlated?
summary(lm(diff.morph$D ~ diff.mode$D))
summary(lm(diff.morph$FRic ~ diff.mode$FRic))

# morphological vs. LH-constant: not very correlated (D highest)
round(diag(cor(diff.morph[, means], diff.constant[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-constant: Moderately positively, except FRic & FEve; mean r = 0.70
round(diag(cor(diff.mode[, means], diff.constant[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-raw: Moderately: Moderately positively, except for D and FDis; mean r = 0.55
round(diag(cor(diff.mode[, means], diff.raw[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-raw vs. LH-constant: Moderately positively: Except for FEve; mean r = 0.49
round(diag(cor(diff.raw[, means], diff.constant[, means], 
               use = "pairwise.complete.obs")), 4)





## COMPARING LIFE-HABIT PROPOGATION ALGORITHMS #################################

# Compare three algorithms (not sample-standardized)
metrics_mode <- read.csv(file = "metrics_LH_mode.csv", header = TRUE)
metrics_constant <- read.csv(file = "metrics_LH_constant.csv", header = TRUE)
metrics_raw <- read.csv(file = "metrics_LH_raw.csv", header = TRUE)

# Plot trends
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
for (c in c(2:10)) {
  var_m <- metrics_mode[ ,c]
  var_c <- metrics_constant[ ,c]
  var_r <- metrics_raw[ ,c]
  lim <- range(c(var_m, var_c, var_r), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
               tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
               cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
               ts.col = TRUE, label = colnames(metrics_mode)[c],
               timescale = ICS2020, type = "n", abbrev = "Period")
  mtext(text = colnames(metrics_mode)[c], side = 3, cex = 1.25)
  lines(mids, var_r, lwd = 3, lty = 1)
  lines(mids, var_m, lwd = 3, lty = 2, col = "blue")
  lines(mids, var_c, lwd = 3, lty = 3, col = "red")
  legend("bottomright", legend = c("raw", "mode", "constant"), 
         col = c("black", "blue", "red"), bty = "n", lty = c(1, 2, 3), 
         lwd = 2, inset = 0.05)
}



# Compare three algorithms (%-standardized)
metrics_mode <- read.csv(file = "metrics_StdG50_LH_mode.csv", header = TRUE)
metrics_constant <- read.csv(file = "metrics_StdG50_LH_constant.csv", header = TRUE)
metrics_raw <- read.csv(file = "metrics_StdG50_LH_raw.csv", header = TRUE)

# Plot trends
# pdf(file = "3EcoTrends_StdG.pdf"); 
par(mar = c(0, 4, 4, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17)
colnames(metrics_mode)[var.cols]       # Confirm plotting correct columns
cols <- viridisLite::plasma(4)
trans.cols <- viridisLite::plasma(4, alpha = 0.5)
for (c in var.cols) {
  var_m <- metrics_mode[ ,c]
  var_c <- metrics_constant[ ,c]
  var_r <- metrics_raw[ ,c]
  var_m_top <- var_m + metrics_mode[ ,c + 1]
  var_m_bottom <- var_m - metrics_mode[ ,c + 1]
  var_m_top <- (var_m_top - min(var_m, na.rm = TRUE)) /
    (max(var_m, na.rm = TRUE) - min(var_m, na.rm = TRUE))
  var_m_bottom <- (var_m_bottom - min(var_m, na.rm = TRUE)) / 
    (max(var_m, na.rm = TRUE) - min(var_m, na.rm = TRUE))
  var_m <- (var_m - min(var_m, na.rm = TRUE)) / 
    (max(var_m, na.rm = TRUE) - min(var_m, na.rm = TRUE))
  var_r <- (var_r - min(var_r, na.rm = TRUE)) / 
    (max(var_r, na.rm = TRUE) - min(var_r, na.rm = TRUE))
  var_c <- (var_c - min(var_c, na.rm = TRUE)) / 
    (max(var_c, na.rm = TRUE) - min(var_c, na.rm = TRUE))
  lim <- range(c(var_m, var_c, var_r), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
                tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
                cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
                ts.col = TRUE, label = colnames(metrics_mode)[c],
                timescale = ICS2020, type = "n", abbrev = "Period")
  mtext(text = paste(colnames(metrics_mode)[c], " (G = ", std.g, ")", sep = ""), 
        side = 3, cex = 1.25)
  column <- cbind(c(mids, rev(mids)), c(var_m_bottom, rev(var_m_top)))
  column <- na.omit(column)
  polygon(column[ ,1], column[ ,2], col = trans.cols[2], lwd = 2, border = NA)
  lines(mids, var_r, lwd = 4, lty = 1, col = cols[1])
  lines(mids, var_m, lwd = 4, lty = 2, col = cols[2])
  lines(mids, var_c, lwd = 4, lty = 3, col = cols[3])
  legend("bottomright", legend = c("raw", "mode", "constant"), col = cols, 
         bty = "n", lty = c(1, 2, 3), lwd = 2, inset = 0.05)
}
par(op)
# dev.off()







## COMPARING MORPHOLOGICAL & ECOLOGICAL TRENDS #################################

# Compare four data sets (not sample-standardized)
metrics_morph <- read.csv(file="metrics_morph.csv", header=TRUE)
metrics_mode <- read.csv(file = "metrics_LH_mode.csv", header = TRUE)
metrics_constant <- read.csv(file = "metrics_LH_constant.csv", header = TRUE)
metrics_raw <- read.csv(file = "metrics_LH_raw.csv", header = TRUE)

# Plot trends (% transformed)
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
cols <- viridisLite::plasma(4)
for (c in c(2:10)) {
  var_md <- metrics_mode[ ,c]
  var_c <- metrics_constant[ ,c]
  var_r <- metrics_raw[ ,c]
  var_mr <- metrics_morph[ ,c]
  var_md <- (var_md - min(var_md, na.rm = TRUE)) / 
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_c <- (var_c - min(var_c, na.rm = TRUE)) / 
    (max(var_c, na.rm = TRUE) - min(var_c, na.rm = TRUE))
  var_r <- (var_r - min(var_r, na.rm = TRUE)) / 
    (max(var_r, na.rm = TRUE) - min(var_r, na.rm = TRUE))
  var_mr <- (var_mr - min(var_mr, na.rm = TRUE)) / 
    (max(var_mr, na.rm = TRUE) - min(var_mr, na.rm = TRUE))
  lim <- range(c(var_md, var_c, var_r, var_mr), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
                tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
                cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
                ts.col = TRUE, timescale = ICS2020,  
                label = paste(colnames(metrics_mode)[c], "(% transformed)"),
                type = "n", abbrev = "Period")
  mtext(text = colnames(metrics_mode)[c], side = 3, cex = 1.25)
  lines(mids, var_mr, lwd = 3, lty = 1, col = cols[1])
  lines(mids, var_md, lwd = 3, lty = 2, col = cols[2])
  lines(mids, var_c, lwd = 3, lty = 3, col = cols[3])
  lines(mids, var_r, lwd = 3, lty = 4, col = cols[4])
  legend("bottomright", legend = c("morph", "eco-mode", "eco-constant", "eco-ra"), 
         col = cols, bty = "n", lty = c(1, 2, 3, 4), lwd = 2, inset = 0.05)
}
par(op)




# Compare three algorithms (sample-standardized) (% transformed)
metrics_mode <- read.csv(file = "metrics_StdG50_LH_mode.csv", header = TRUE)
metrics_constant <- read.csv(file = "metrics_StdG50_LH_constant.csv", header = TRUE)
metrics_raw <- read.csv(file = "metrics_StdG50_LH_raw.csv", header = TRUE)
metrics_morph <- read.csv(file = "metrics_StdG50_morph.csv", header = TRUE)
std.g <- 50

# Plot trends
# pdf(file = "Morph&3EcoTrends_StdG.pdf")
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17)
colnames(metrics_mode)[var.cols]       # Confirm plotting correct columns
cols <- viridisLite::plasma(4)         # Uses RGB-sensitive red and blue
trans.cols <- viridisLite::plasma(4, alpha = 0.5)
for (c in var.cols) {
  var_md <- metrics_mode[ ,c]
  var_c <- metrics_constant[ ,c]
  var_r <- metrics_raw[ ,c]
  var_mr <- metrics_morph[ ,c]
  var_md_top <- var_md + metrics_mode[ ,c + 1]
  var_md_bottom <- var_md - metrics_mode[ ,c + 1]
  var_md_top <- (var_md_top - min(var_md, na.rm = TRUE)) /
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_md_bottom <- (var_md_bottom - min(var_md, na.rm = TRUE)) / 
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_md <- (var_md - min(var_md, na.rm = TRUE)) / 
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_r <- (var_r - min(var_r, na.rm = TRUE)) / 
    (max(var_r, na.rm = TRUE) - min(var_r, na.rm = TRUE))
  var_c <- (var_c - min(var_c, na.rm = TRUE)) / 
    (max(var_c, na.rm = TRUE) - min(var_c, na.rm = TRUE))
  # Sample-standardization made H constant (=StdG) for morphological data set,
  # so treated differently to force to equal S (or 1 when %-max transformed)
  if (length(unique(na.omit(var_mr))) ==  1L)
    var_mr <- replace(na.omit(var_mr), !is.na(var_mr), 1) else
    var_mr <- (var_mr - min(var_mr, na.rm = TRUE)) / 
      (max(var_mr, na.rm = TRUE) - min(var_mr, na.rm = TRUE))
  lim <- range(c(var_md, var_c, var_r, var_mr), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
                tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
                cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
                ts.col = TRUE, timescale = ICS2020, type = "n", abbrev = "Period",
                label = paste(colnames(metrics_mode)[c], "(% transformed)"))
  mtext(text = paste(colnames(metrics_mode)[c], " (G = ", std.g, ")", sep = ""), 
        side = 3, cex = 1.25)
  column <- cbind(c(mids, rev(mids)), c(var_md_bottom, rev(var_md_top)))
  column <- na.omit(column)
  polygon(column[ ,1], column[ ,2], col = trans.cols[2], lwd = 2, border = NA)
  lines(mids, var_mr, lwd = 4, lty = 1, col = cols[1])
  lines(mids, var_md, lwd = 4, lty = 2, col = cols[2])
  lines(mids, var_c, lwd = 4, lty = 3, col = cols[3])
  lines(mids, var_r, lwd = 4, lty = 4, col = cols[4])
  legend("bottomright", legend = c("morph", "eco-mode", "eco-constant", "eco-raw"), 
         col = cols, bty = "n", lty = c(1, 2, 3, 4), lwd = 2, inset = 0.05)
}
par(op)
# dev.off()




# Same, but just mode and morph
# pdf(file = "2Morph&EcoTrends_StdG.pdf")
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17)
colnames(metrics_mode)[var.cols]       # Confirm plotting correct columns
cols <- viridisLite::plasma(3)         # Uses RGB-sensitive red and blue
trans.cols <- viridisLite::plasma(3, alpha = 0.5)
for (c in var.cols) {
  var_md <- metrics_mode[ ,c]
  var_mr <- metrics_morph[ ,c]
  var_md_top <- var_md + metrics_mode[ ,c + 1]
  var_md_bottom <- var_md - metrics_mode[ ,c + 1]
  var_mr_top <- var_mr + metrics_morph[ ,c + 1]
  var_mr_bottom <- var_mr - metrics_morph[ ,c + 1]
  var_mr_top <- (var_mr_top - min(var_mr, na.rm = TRUE)) /
    (max(var_mr, na.rm = TRUE) - min(var_mr, na.rm = TRUE))
  var_mr_bottom <- (var_mr_bottom - min(var_mr, na.rm = TRUE)) / 
    (max(var_mr, na.rm = TRUE) - min(var_mr, na.rm = TRUE))
  var_mr <- (var_mr - min(var_mr, na.rm = TRUE)) / 
    (max(var_mr, na.rm = TRUE) - min(var_mr, na.rm = TRUE))
  var_md_top <- (var_md_top - min(var_md, na.rm = TRUE)) /
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_md_bottom <- (var_md_bottom - min(var_md, na.rm = TRUE)) / 
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  var_md <- (var_md - min(var_md, na.rm = TRUE)) / 
    (max(var_md, na.rm = TRUE) - min(var_md, na.rm = TRUE))
  # Sample-standardization made H constant (=StdG) for morphological data set,
  # so treated differently to force to equal S (or 1 when %-max transformed).
  # Because constant, no error bar.
  if (length(na.omit(var_mr)) ==  0L)
    var_mr <- replace(var_mr, is.nan(var_mr), 1)
  lim <- range(c(var_md, var_mr), na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
                tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
                cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
                ts.col = TRUE, timescale = ICS2020, type = "n", abbrev = "Period",
                label = paste(colnames(metrics_mode)[c], "(% transformed)"))
  mtext(text = paste(colnames(metrics_mode)[c], " (G = ", std.g, ")", sep = ""), 
        side = 3, cex = 1.25)
  # Add error bars:
  column_md <- cbind(c(mids, rev(mids)), c(var_md_bottom, rev(var_md_top)))
  column_md <- na.omit(column_md)
  polygon(column_md[ ,1], column_md[ ,2], col = trans.cols[1], lwd = 2, border = NA)
  column_mr <- cbind(c(mids, rev(mids)), c(var_mr_bottom, rev(var_mr_top)))
  column_mr <- na.omit(column_mr)
  polygon(column_mr[ ,1], column_mr[ ,2], col = trans.cols[2], lwd = 2, border = NA)
  # Overlay trend lines
  lines(mids, var_md, lwd = 4, lty = 1, col = cols[1])
  lines(mids, var_mr, lwd = 4, lty = 2, col = cols[2])
  # Only add legend to first plot:
  if (length(unique(na.omit(var_mr))) ==  1L)
    legend("bottomright", legend = c("ecology", "morphology"), col = cols, 
           bty = "n", lty = c(1, 2), lwd = 4, inset = 0.02, cex = 1.5)
}
par(op)
# dev.off()








## PLOT PHYLOMORPHOSPACE AND PHYLOECOSPACE #####################################
## (code modified from Brad Deline, bdeline@westga.edu)

# Import time-scaled "equal" phylogeny
load("~/Manuscripts/CamOrdEchinos/equal.tree")
tree <- equal.tree

# Load Wills GED-50 distance matrices
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")

# Load PCoA output from ape::pcoa
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")
load("morph.pcoa")

# Set phylomorphospace/phyloecospace plotting colors
par(mar = c(5, 4, 2, 2))
branch.col <- "gray50"
tip.col <- "#0000FF7F"   # Set blue transparent so overlays as density
node.col <- "#A020F07F"  # Set purple transparent so overlays as density
root.col <- "#FDE725FF"  # viridisLite::viridis(3)[3]

# Set phylomorphospace/phyloecospace plotting variables
Tree <- morph.pcoa$Tree
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))
x.lab1 <- paste0("PCO 1 (", round(100 * morph.pcoa$values$Rel_corr_eig[1], 1), 
                 "% of total variance)")
y.lab2 <- paste0("PCO 2 (", round(100 * morph.pcoa$values$Rel_corr_eig[2], 1), 
                 "% of total variance)")
x.lab3 <- paste0("PCO 3 (", round(100 * morph.pcoa$values$Rel_corr_eig[3], 1), 
                 "% of total variance)")
y.lab4 <- paste0("PCO 4 (", round(100 * morph.pcoa$values$Rel_corr_eig[4], 1), 
                 "% of total variance)")
x.lab5 <- paste0("PCO 5 (", round(100 * morph.pcoa$values$Rel_corr_eig[5], 1), 
                 "% of total variance)")
y.lab6 <- paste0("PCO 6 (", round(100 * morph.pcoa$values$Rel_corr_eig[6], 1), 
                 "% of total variance)")

# pdf(file = "morphospace.pdf")
par(mfrow = c(2, 2))
# PCO 1 vs. PCO 2
phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 1:2], 
                           A = morph.pcoa$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = x.lab1, 
                           ylab = y.lab2, pch = NA)
points(x = morph.pcoa$vectors.cor[node.seq, 1],
       y = morph.pcoa$vectors.cor[node.seq, 2], col = node.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[tip.seq, 1], 
       y = morph.pcoa$vectors.cor[tip.seq, 2], col = tip.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 1], 
       y = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 2], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 3 vs. PCO 4
phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 3:4], 
                           A = morph.pcoa$vectors.cor[node.seq, 3:4], 
                           control = con, label = "off", xlab = x.lab3, 
                           ylab = y.lab4, pch = NA)
points(x = morph.pcoa$vectors.cor[node.seq, 3],
       y = morph.pcoa$vectors.cor[node.seq, 4], col = node.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[tip.seq, 3], 
       y = morph.pcoa$vectors.cor[tip.seq, 4], col = tip.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 3], 
       y = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 4], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 5 vs. PCO 6
phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 5:6], 
                           A = morph.pcoa$vectors.cor[node.seq, 5:6], 
                           control = con, label = "off", xlab = x.lab5, 
                           ylab = y.lab6, pch = NA)
points(x = morph.pcoa$vectors.cor[node.seq, 5],
       y = morph.pcoa$vectors.cor[node.seq, 6], col = node.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[tip.seq, 5], 
       y = morph.pcoa$vectors.cor[tip.seq, 6], col = tip.col, pch = 16, cex = 1.25)
points(x = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 5], 
       y = morph.pcoa$vectors.cor[(max(tip.seq) + 1), 6], col = root.col, pch = 16, 
       cex = 1.25)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "phylomorphospace", legend = c("tips", "nodes", "root"), 
       cex = 2, pch = 16, col = c(tip.col, node.col, root.col), bty = "n", 
       pt.cex = 4)
# dev.off()
par(op)


# First 3 in 3-D (wait a few seconds while renders)
phylomorphospace3d(tree = morph.pcoa$Tree, X = morph.pcoa$vectors.cor[tip.seq, 1:3], 
                 A = morph.pcoa$vectors.cor[node.seq, 1:3], 
                 control = c(con, list(ftype = "off")))
par(op)


## Same, with 'mode' ecological data set
# Set phylomorphospace/phylomorphospace plotting variables
Tree <- mode.pcoa$Tree
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))
x.lab1 <- paste0("PCO 1 (", round(100 * mode.pcoa$values$Rel_corr_eig[1], 1), 
                 "% of total variance)")
y.lab2 <- paste0("PCO 2 (", round(100 * mode.pcoa$values$Rel_corr_eig[2], 1), 
                 "% of total variance)")
x.lab3 <- paste0("PCO 3 (", round(100 * mode.pcoa$values$Rel_corr_eig[3], 1), 
                 "% of total variance)")
y.lab4 <- paste0("PCO 4 (", round(100 * mode.pcoa$values$Rel_corr_eig[4], 1), 
                 "% of total variance)")
x.lab5 <- paste0("PCO 5 (", round(100 * mode.pcoa$values$Rel_corr_eig[5], 1), 
                 "% of total variance)")
y.lab6 <- paste0("PCO 6 (", round(100 * mode.pcoa$values$Rel_corr_eig[6], 1), 
                 "% of total variance)")

# pdf(file = "ecospace_mode.pdf")
par(mfrow = c(2, 2))
# PCO 1 vs. PCO 2
phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 1:2], 
                           A = mode.pcoa$vectors.cor[node.seq, 1:2], 
                           control = con, label = "off", xlab = x.lab1, 
                           ylab = y.lab2, pch = NA)
points(x = mode.pcoa$vectors.cor[node.seq, 1],
       y = mode.pcoa$vectors.cor[node.seq, 2], col = node.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[tip.seq, 1], 
       y = mode.pcoa$vectors.cor[tip.seq, 2], col = tip.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 1], 
       y = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 2], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 3 vs. PCO 4
phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 3:4], 
                           A = mode.pcoa$vectors.cor[node.seq, 3:4], 
                           control = con, label = "off", xlab = x.lab3, 
                           ylab = y.lab4, pch = NA)
points(x = mode.pcoa$vectors.cor[node.seq, 3],
       y = mode.pcoa$vectors.cor[node.seq, 4], col = node.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[tip.seq, 3], 
       y = mode.pcoa$vectors.cor[tip.seq, 4], col = tip.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 3], 
       y = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 4], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 5 vs. PCO 6
phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 5:6], 
                           A = mode.pcoa$vectors.cor[node.seq, 5:6], 
                           control = con, label = "off", xlab = x.lab5, 
                           ylab = y.lab6, pch = NA)
points(x = mode.pcoa$vectors.cor[node.seq, 5],
       y = mode.pcoa$vectors.cor[node.seq, 6], col = node.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[tip.seq, 5], 
       y = mode.pcoa$vectors.cor[tip.seq, 6], col = tip.col, pch = 16, cex = 1.25)
points(x = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 5], 
       y = mode.pcoa$vectors.cor[(max(tip.seq) + 1), 6], col = root.col, pch = 16, 
       cex = 1.25)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "phyloecospace (mode)", legend = c("tips", "nodes", "root"), 
       cex = 2, pch = 16, col = c(tip.col, node.col, root.col), bty = "n", 
       pt.cex = 4)
# dev.off()


# First 3 in 3-D (wait a few seconds while renders)
phylomorphospace3d(tree = mode.pcoa$Tree, X = mode.pcoa$vectors.cor[tip.seq, 1:3], 
                   A = mode.pcoa$vectors.cor[node.seq, 1:3], 
                   control = c(con, list(ftype = "off")))
par(op)





## Same, with 'constant' ecological data set
# Set phylomorphospace/phylomorphospace plotting variables
# Plotting uncorrected eigenvectors instead of corrected ones
Tree <- constant.pcoa$Tree
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))
x.lab1 <- paste0("PCO 1 (", round(100 * constant.pcoa$values$Relative_eig[1], 1), 
                 "% of total variance)")
y.lab2 <- paste0("PCO 2 (", round(100 * constant.pcoa$values$Relative_eig[2], 1), 
                 "% of total variance)")
x.lab3 <- paste0("PCO 3 (", round(100 * constant.pcoa$values$Relative_eig[3], 1), 
                 "% of total variance)")
y.lab4 <- paste0("PCO 4 (", round(100 * constant.pcoa$values$Relative_eig[4], 1), 
                 "% of total variance)")
x.lab5 <- paste0("PCO 5 (", round(100 * constant.pcoa$values$Relative_eig[5], 1), 
                 "% of total variance)")
y.lab6 <- paste0("PCO 6 (", round(100 * constant.pcoa$values$Relative_eig[6], 1), 
                 "% of total variance)")

# pdf(file = "ecospace_constant.pdf")
par(mfrow = c(2, 2))
# PCO 1 vs. PCO 2
phytools::phylomorphospace(tree = Tree, X = constant.pcoa$vectors[tip.seq, 1:2],
                           A = constant.pcoa$vectors[node.seq, 1:2],
                           control = con, label = "off", xlab = x.lab1,
                           ylab = y.lab2, pch = NA)
points(x = constant.pcoa$vectors[node.seq, 1],
       y = constant.pcoa$vectors[node.seq, 2], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[tip.seq, 1], 
       y = constant.pcoa$vectors[tip.seq, 2], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[(max(tip.seq) + 1), 1], 
       y = constant.pcoa$vectors[(max(tip.seq) + 1), 2], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 3 vs. PCO 4
phytools::phylomorphospace(tree = Tree, X = constant.pcoa$vectors[tip.seq, 3:4],
                           A = constant.pcoa$vectors[node.seq, 3:4],
                           control = con, label = "off", xlab = x.lab3,
                           ylab = y.lab4, pch = NA)
points(x = constant.pcoa$vectors[node.seq, 3],
       y = constant.pcoa$vectors[node.seq, 4], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[tip.seq, 3], 
       y = constant.pcoa$vectors[tip.seq, 4], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[(max(tip.seq) + 1), 3], 
       y = constant.pcoa$vectors[(max(tip.seq) + 1), 4], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 5 vs. PCO 6
phytools::phylomorphospace(tree = Tree, X = constant.pcoa$vectors[tip.seq, 5:6],
                           A = constant.pcoa$vectors[node.seq, 5:6],
                           control = con, label = "off", xlab = x.lab5,
                           ylab = y.lab6, pch = NA)
points(x = constant.pcoa$vectors[node.seq, 5],
       y = constant.pcoa$vectors[node.seq, 6], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[tip.seq, 5], 
       y = constant.pcoa$vectors[tip.seq, 6], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors[(max(tip.seq) + 1), 5], 
       y = constant.pcoa$vectors[(max(tip.seq) + 1), 6], col = root.col, pch = 16, 
       cex = 1.25)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "phyloecospace (constant)", legend = c("tips", "nodes", "root"), 
       cex = 2, pch = 16, col = c(tip.col, node.col, root.col), bty = "n", 
       pt.cex = 4)
# dev.off()
par(op)


# First 3 in 3-D (wait a few seconds while renders)
phylomorphospace3d(tree = constant.pcoa$Tree, X = constant.pcoa$vectors.cor[tip.seq, 1:3], 
                A = constant.pcoa$vectors.cor[node.seq, 1:3], 
                control = c(con, list(ftype = "off")))
par(op)





## Same, with 'raw' ecological data set

# Remove missing taxa (Only used for 'raw' data treatment)
trim.matrix <- TrimMorphDistMatrix(raw.distances.GED.5$distance_matrix, 
                                   Tree = tree)
raw.distances.GED.5$distance_matrix <- trim.matrix$DistMatrix
raw.distances.GED.5$Tree <- trim.matrix$Tree
raw.distances.GED.5$RemovedTaxa <- trim.matrix$RemovedTaxa

# Set phylomorphospace/phylomorphospace plotting variables
Tree <- raw.distances.GED.5$Tree  # Note using the trimmed tree here (unlike above)
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))
x.lab1 <- paste0("PCO 1 (", round(100 * raw.pcoa$values$Rel_corr_eig[1], 1), 
                 "% of total variance)")
y.lab2 <- paste0("PCO 2 (", round(100 * raw.pcoa$values$Rel_corr_eig[2], 1), 
                 "% of total variance)")
x.lab3 <- paste0("PCO 3 (", round(100 * raw.pcoa$values$Rel_corr_eig[3], 1), 
                 "% of total variance)")
y.lab4 <- paste0("PCO 4 (", round(100 * raw.pcoa$values$Rel_corr_eig[4], 1), 
                 "% of total variance)")
x.lab5 <- paste0("PCO 5 (", round(100 * raw.pcoa$values$Rel_corr_eig[5], 1), 
                 "% of total variance)")
y.lab6 <- paste0("PCO 6 (", round(100 * raw.pcoa$values$Rel_corr_eig[6], 1), 
                 "% of total variance)")

# pdf(file = "ecospace_raw.pdf")
par(mfrow = c(2, 2))
# PCO 1 vs. PCO 2
phytools::phylomorphospace(tree = Tree, X = raw.pcoa$vectors.cor[tip.seq, 1:2],
                           A = raw.pcoa$vectors.cor[node.seq, 1:2],
                           control = con, label = "off", xlab = x.lab1,
                           ylab = y.lab2, pch = NA)
points(x = raw.pcoa$vectors.cor[node.seq, 1],
       y = raw.pcoa$vectors.cor[node.seq, 2], col = node.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[tip.seq, 1], 
       y = raw.pcoa$vectors.cor[tip.seq, 2], col = tip.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 1], 
       y = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 2], col = root.col, pch = 16, cex = 1.25)

# PCO 3 vs. PCO 4
phytools::phylomorphospace(tree = Tree, X = raw.pcoa$vectors.cor[tip.seq, 3:4], 
                           A = raw.pcoa$vectors.cor[node.seq, 3:4],
                           control = con, label = "off", xlab = x.lab3,
                           ylab = y.lab4, pch = NA)
points(x = raw.pcoa$vectors.cor[node.seq, 3],
       y = raw.pcoa$vectors.cor[node.seq, 4], col = node.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[tip.seq, 3], 
       y = raw.pcoa$vectors.cor[tip.seq, 4], col = tip.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 3], 
       y = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 4], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 5 vs. PCO 6
phytools::phylomorphospace(tree = Tree, X = raw.pcoa$vectors.cor[tip.seq, 5:6], 
                           A = raw.pcoa$vectors.cor[node.seq, 5:6],
                           control = con, label = "off", xlab = x.lab5,
                           ylab = y.lab6, pch = NA)
points(x = raw.pcoa$vectors.cor[node.seq, 5],
       y = raw.pcoa$vectors.cor[node.seq, 6], col = node.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[tip.seq, 5], 
       y = raw.pcoa$vectors.cor[tip.seq, 6], col = tip.col, pch = 16, cex = 1.25)
points(x = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 5], 
       y = raw.pcoa$vectors.cor[(max(tip.seq) + 1), 6], col = root.col, pch = 16, 
       cex = 1.25)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "phyloecospace (raw)", legend = c("tips", "nodes", "root"), 
       cex = 2, pch = 16, col = c(tip.col, node.col, root.col), bty = "n", 
       pt.cex = 4)
# dev.off()
par(op)



# First 3 in 3-D (wait a few seconds while renders)
phylomorphospace3d(tree = raw.distances.GED.5$Tree, 
                   X = raw.pcoa$vectors.cor[tip.seq, 1:3], 
                   A = raw.pcoa$vectors.cor[node.seq, 1:3],
                   control = c(con, list(ftype = "off")))
par(op)




## KMEANS CLUSTERING TO INTERPRET SPACES #######################################

load("taxon.list")
if(nrow(morph.pcoa$vectors.cor) != length(taxon.list))
  stop("Need to re-build taxon.list because some taxa were removed when running the 'raw' treatment.\n")

# Plot within-cluster sum-of-squares vs no. of clusters
# pdf(file = "morph.kmeans.choice.pdf")
no.axes <- 6
ks <- c(1:10, seq(20, 50, by = 10), 75, 100)
sum.sq <- array(dim = length(ks))
for(k in 1:length(ks)){
  sum.sq[k] <- kmeans(morph.pcoa$vectors.cor[, 1:no.axes], centers = ks[k],
                      nstart = 25, iter.max = 1000)$tot.withinss
}
plot(ks, sum.sq, type = "b", pch = 16, xlab = "No. k-means clusters",
     ylab = "Total within-cluster sum of squares",
     main = "Morphological PCoA k-means clusters")
par(op)
# dev.off()

# For morphology, k = 4 & 5 provide clean breaks; 9 divides most classes
# pdf(file = "kmeans_morph.pdf")
k <- 4
set.seed(1234) # To allow replication of order
km <- kmeans(morph.pcoa$vectors.cor[, 1:6], centers = k, nstart = 25, iter.max = 100)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25))
cols <- plasma(k)[km$cluster]
pchs <-as.character(km$cluster)
plot(morph.pcoa$vectors.cor[, 1:2], col = cols, pch = pchs, cex = 0.75)
plot(morph.pcoa$vectors.cor[, 3:4], col = cols, pch = pchs, cex = 0.75)
plot(morph.pcoa$vectors.cor[, 5:6], col = cols, pch = pchs, cex = 0.75)
table(taxon.list, km$cluster) # If phylogenetic structure, lots of zeros
legend.groups <- c("edrio, aster, ech, oph & cyclo", 
                   "stylo, homo, solut, cteno, helico",
                   "crin, eocr, rhomb, diplo, paracr", "more crinoids")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "morphological PCoA", legend = legend.groups, cex = 1.25,
       pch = as.character(1:k), col = plasma(k)[1:k], bty = "n", pt.cex = 1.5)
# dev.off()
par(op)

# Plot within-cluster sum-of-squares vs no. of clusters
# pdf(file = "mode.kmeans.choice.pdf")
no.axes <- 4
ks <- c(1:10, seq(20, 50, by = 10), 75, 100)
sum.sq <- array(dim = length(ks))
for(k in 1:length(ks)){
  sum.sq[k] <- kmeans(mode.pcoa$vectors.cor[, 1:no.axes], centers = ks[k],
                      nstart = 25, iter.max = 1000)$tot.withinss
}
plot(ks, sum.sq, type = "b", pch = 16, xlab = "No. k-means clusters",
     ylab = "Total within-cluster sum of squares",
     main = "Ecological (mode) PCoA k-means clusters")
par(op)
# dev.off()

# For mode, k = 3-5 provides clean breaks; 9 divides most classes
# pdf(file = "kmeans_eco_mode.pdf")
k <- 4
set.seed(1) # To allow replication of order
km <- kmeans(mode.pcoa$vectors.cor[, 1:6], centers = k, nstart = 25, iter.max = 100)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25))
cols <- plasma(k)[km$cluster]
pchs <-as.character(km$cluster)
plot(mode.pcoa$vectors.cor[, 1:2], col = cols, pch = pchs, cex = 0.75)
plot(mode.pcoa$vectors.cor[, 3:4], col = cols, pch = pchs, cex = 0.75)
plot(mode.pcoa$vectors.cor[, 5:6], col = cols, pch = pchs, cex = 0.75)
table(taxon.list, km$cluster) # If phylogenetic structure, lots of zeros
legend.groups <- c("crinoids (& some eocr, edrio & rhomb)",
                   "edr, eocr, rh, dipl, paracr, hel & more cri",
                   "stylo, homo, solut, ctenoc, cyclo, some rhomb",
                   "aster, ech, oph, holo, somas & stenur")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "Ecological PCoA", legend = legend.groups, cex = 1,
       pch = as.character(1:k), col = plasma(k)[1:k], bty = "n", pt.cex = 1.5)
# dev.off()
par(op)



# Plot within-cluster sum-of-squares vs no. of clusters
# pdf(file = "constant.kmeans.choice.pdf")
# Note using uncorrected eigenvectors
no.axes <- 6
ks <- c(1:10, seq(20, 50, by = 10), 75, 100)
sum.sq <- array(dim = length(ks))
for(k in 1:length(ks)){
  sum.sq[k] <- kmeans(constant.pcoa$vectors[, 1:no.axes], centers = ks[k],
                      nstart = 25, iter.max = 1000)$tot.withinss
}
plot(ks, sum.sq, type = "b", pch = 16, xlab = "No. k-means clusters",
     ylab = "Total within-cluster sum of squares",
     main = "Ecological (constant) PCoA k-means clusters")
par(op)
# dev.off()

# For constant, k = 3 (or 4) provides clean breaks; 9 divides most classes
# pdf(file = "kmeans_eco_constant.pdf")
# Note using uncorrected eigenvectors
k <- 4
set.seed(3) # To allow replication of order
km <- kmeans(constant.pcoa$vectors[, 1:6], centers = k, nstart = 25, iter.max = 100)
# Modify so matches 'mode' (switch 2 and 3)
k.3 <- which(km$cluster == 3)
k.4 <- which(km$cluster == 4)
km$cluster[k.3] <- 4
km$cluster[k.4] <- 3
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25))
cols <- plasma(k)[km$cluster]
pchs <-as.character(km$cluster)
plot(constant.pcoa$vectors[, 1:2], col = cols, pch = pchs, cex = 0.75)
plot(constant.pcoa$vectors[, 3:4], col = cols, pch = pchs, cex = 0.75)
plot(constant.pcoa$vectors[, 5:6], col = cols, pch = pchs, cex = 0.75)
table(taxon.list, km$cluster) # If phylogenetic structure, lots of zeros
legend.groups <- c("crinoids (& some eocr, edrio & rhomb)",
                   "edr, eocr, rh, dipl, paracr, hel & more cri",
                   "stylo, homo, solut, ctenoc, cyclo, some rhomb",
                   "aster, ech, oph, holo, somas & stenur")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "Ecological PCoA", legend = legend.groups, cex = 1,
       pch = as.character(1:k), col = plasma(k)[1:k], bty = "n", pt.cex = 1.5)
# dev.off()
par(op)






## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF CLASSES ################
# Question: Within taxonomic classes, is the overall occupation of space greater
# in morphospace or ecospace?

# No need to sample standardize here, as the sample size is the same for
# ecological and morphological subsets. Simply need to observe the differences
# in convex hull hypervolume (FRic), using 'm' PCoA axes like above. And this
# analysis is not observing the trends through time. May need to standardize
# each space to the overall FRic of all echinoderms, as the length of
# eigenvalues may be different in each space.

# Select 'm' (number of PCoA axes to use when calculating convex hull volume)
m <- 2

# Process the classes:
load("taxon.list")
class.tab <- sort(table(taxon.list), decreasing = TRUE)
classes <- names(class.tab)
# Remove genera with UNCERTAIN class affiliation
classes <- classes[-which(classes == "UNCERTAIN")]
class.space <- data.frame(class = classes, morph = NA, eco = NA)

# Process PCoA output
load("morph.pcoa")
load("mode.pcoa")
eco.pcoa <- mode.pcoa
# load("constant.pcoa"); eco.pcoa <- constant.pcoa
# load("raw.pcoa"); eco.pcoa <- raw.pcoa
# Standardize the PCoA eigenvectors
morph.pcoa$vectors.cor <- stand.pcoa(vectors = morph.pcoa$vectors.cor, 
                                     eigenvalues = morph.pcoa$values$Corr_eig)
eco.pcoa$vectors.cor <- stand.pcoa(vectors = eco.pcoa$vectors.cor, 
                                   eigenvalues = eco.pcoa$values$Corr_eig)

for(cl in 1:nrow(class.space)) {
  wh.cl <- which(taxon.list == classes[cl])
  eco.coords <- eco.pcoa$vectors.cor[wh.cl, 1:m]
  morph.coords <- morph.pcoa$vectors.cor[wh.cl, 1:m]
  if (length(eco.coords) == 2L) next
  # Can not calculate if there are less than 'm' distinct points
  H.eco <- nrow(unique(round(eco.coords, 5)))
  H.morph <- nrow(unique(round(morph.coords, 5)))
  if (H.eco <= m | H.morph <= m)
    next
  class.space$morph[cl] <- geometry::convhulln(morph.coords, "FA")$vol
  class.space$eco[cl] <- geometry::convhulln(eco.coords, "FA")$vol
}

# Standardize by entire morphospace / ecospace
morph.total <- geometry::convhulln(morph.pcoa$vectors.cor[, 1:m], "FA")$vol
eco.total <- geometry::convhulln(eco.pcoa$vectors.cor[, 1:m], "FA")$vol
class.space$morph <- class.space$morph / morph.total
class.space$eco <- class.space$eco / eco.total

# Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
class.space$diff <- class.space$morph - class.space$eco

class.space
summary(class.space$diff)

# Statistical test
par(mfrow = c(2, 1))
hist(class.space$diff, 15, main = "differences in class convex hull volume", 
     xlab = "morph volume - eco volume, by class")
abline(v = 0, lwd = 2, lty = 2)
plot(density(na.omit(class.space$diff)), main = "differences in class convex hull volume", 
     xlab = "morph volume - eco volume, by class")
abline(v = 0, lwd = 2, lty = 2)

# Although visually almost normally distributed, not technically normally
# distributed (b/c of overly long tails)
shapiro.test(class.space$diff) # p < 0.05 so non-normally distributed
 
wilcox.test(class.space$diff, mu = 0)
t.test(class.space$diff, mu = 0)

# Confirm with paired t-test
wilcox.test(class.space$morph, class.space$eco, paired = TRUE)
t.test(class.space$morph, class.space$eco, paired = TRUE)





## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF SUBPHYLA ###############
# Question: Within taxonomic subphyla, is the overall occupation of space
# greater in morphospace or ecospace?

# Select 'm' (number of PCoA axes to use when calculating convex hull volume)
m <- 2

# Process the subphyla:
load("subphylum.list")
subphylum.tab <- sort(table(subphylum.list), decreasing = TRUE)
subphyla <- names(subphylum.tab)
# Remove genera with UNCERTAIN subphylum affiliation
subphyla <- subphyla[-which(subphyla == "UNCERTAIN")]
subphylum.space <- data.frame(subphylum = subphyla, morph = NA, eco = NA)

# Process PCoA output
load("morph.pcoa")
load("mode.pcoa")
eco.pcoa <- mode.pcoa
# load("constant.pcoa"); eco.pcoa <- constant.pcoa
# load("raw.pcoa"); eco.pcoa <- raw.pcoa
# Standardize the PCoA eigenvectors
morph.pcoa$vectors.cor <- stand.pcoa(vectors = morph.pcoa$vectors.cor, 
                                     eigenvalues = morph.pcoa$values$Corr_eig)
eco.pcoa$vectors.cor <- stand.pcoa(vectors = eco.pcoa$vectors.cor, 
                                   eigenvalues = eco.pcoa$values$Corr_eig)

for(subph in 1:nrow(subphylum.space)) {
  wh.subph <- which(subphylum.list == subphyla[subph])
  eco.coords <- eco.pcoa$vectors.cor[wh.subph, 1:m]
  morph.coords <- morph.pcoa$vectors.cor[wh.subph, 1:m]
  if (length(eco.coords) == 2L) next
  # Can not calculate if there are less than 'm' distinct points
  H.eco <- nrow(unique(round(eco.coords, 5)))
  H.morph <- nrow(unique(round(morph.coords, 5)))
  if (H.eco <= m | H.morph <= m)
    next
  subphylum.space$morph[subph] <- geometry::convhulln(morph.coords, "FA")$vol
  subphylum.space$eco[subph] <- geometry::convhulln(eco.coords, "FA")$vol
}

# Standardize by entire morphospace / ecospace
morph.total <- geometry::convhulln(morph.pcoa$vectors.cor[, 1:m], "FA")$vol
eco.total <- geometry::convhulln(eco.pcoa$vectors.cor[, 1:m], "FA")$vol
subphylum.space$morph <- subphylum.space$morph / morph.total
subphylum.space$eco <- subphylum.space$eco / eco.total

# Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
subphylum.space$diff <- subphylum.space$morph - subphylum.space$eco

subphylum.space
summary(subphylum.space$diff)

# Statistical test
par(mfrow = c(2, 1))
breaks <- seq(from = -.35, to = 0.05, by = 0.05)
hist(subphylum.space$diff, breaks = breaks, 
     main = "differences in subphylum convex hull volume", 
     xlab = "morph volume - eco volume, by subphylum")
abline(v = 0, lwd = 2, lty = 2)
plot(density(na.omit(subphylum.space$diff)), 
     main = "differences in subphylum convex hull volume", 
     xlab = "morph volume - eco volume, by subphylum")
abline(v = 0, lwd = 2, lty = 2)

# Confirm not normally distributed,
shapiro.test(subphylum.space$diff) # p < 0.05 so non-normally distributed

wilcox.test(subphylum.space$diff, mu = 0)
t.test(subphylum.space$diff, mu = 0)

# Confirm with paired t-test
wilcox.test(subphylum.space$morph, subphylum.space$eco, paired = TRUE)
t.test(subphylum.space$morph, subphylum.space$eco, paired = TRUE)




## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF PHYLUM ###############
# Question: Within all echinoderms, is the OVERALL occupation of space greater
# in morphospace or ecospace?

# Select 'm' (number of PCoA axes to use when calculating convex hull volume)
m <- 2

# Process PCoA output
load("morph.pcoa")
load("mode.pcoa")
eco.pcoa <- mode.pcoa
# load("constant.pcoa"); eco.pcoa <- constant.pcoa
# load("raw.pcoa"); eco.pcoa <- raw.pcoa
# Standardize the PCoA eigenvectors
morph.pcoa$vectors.cor <- stand.pcoa(vectors = morph.pcoa$vectors.cor, 
                                     eigenvalues = morph.pcoa$values$Corr_eig)
eco.pcoa$vectors.cor <- stand.pcoa(vectors = eco.pcoa$vectors.cor, 
                                   eigenvalues = eco.pcoa$values$Corr_eig)
phylum.space <- data.frame(morph = NA, eco = NA)

eco.coords <- eco.pcoa$vectors.cor[, 1:m]
morph.coords <- morph.pcoa$vectors.cor[, 1:m]
# Can not calculate if there are less than 'm' distinct points
H.eco <- nrow(unique(round(eco.coords, 5)))
H.morph <- nrow(unique(round(morph.coords, 5)))
if (H.eco <= m | H.morph <= m)
  stop("Convex hull can not be calculated for this 'm'\n")
phylum.space$morph <- geometry::convhulln(morph.coords, "FA")$vol
phylum.space$eco <- geometry::convhulln(eco.coords, "FA")$vol

# Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
phylum.space$diff <- phylum.space$morph - phylum.space$eco

phylum.space
# Overall morphospace is 2.95-times greater than ecospace





## PLOT PHYLOMORPHO/ECOSPACES THROUGH TIME #####################################

# Plot PCO through time, plotting both spaces to examine convergence. Only
# plotting the 'mode' ecology treatment.

# Note that for SI in manuscript, need to reconstruct a different 'taxon.bin'
# object above, but using strat_names$scale_level == 4 so illustrated for
# epoch-level.

load("taxon.bins")

# Set phylomorphospace/phyloecospace plotting variables
# pdf(file = "Phylospaces_through_time.pdf")
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25), pty = "m")
branch.col <- "gray75"
tip.col <- "#0000FF7F"   # Set blue transparent so overlays as density
node.col <- "#A020F07F"  # Set purple transparent so overlays as density
Tree <- morph.pcoa$Tree  # Same tree topology for both morphology and ecology
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))

nt <- ncol(taxon.bins)
for(t in nt:1) {
  wh.gr <- unname(which(taxon.bins[, t]))
  wh.tip <- wh.gr[which(wh.gr <= Ntip(Tree))]
  wh.node <- wh.gr[which(wh.gr > Ntip(Tree))]

  phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 1:2], 
                             A = morph.pcoa$vectors.cor[node.seq, 1:2], 
                             control = con, label = "off", xlab = "PCO 1", 
                             ylab = "PCO 2", pch = NA)
  mtext(paste(colnames(taxon.bins)[t], "morphospace"))
  points(x = morph.pcoa$vectors.cor[wh.node, 1],
         y = morph.pcoa$vectors.cor[wh.node, 2], col = node.col, pch = 16)
  points(x = morph.pcoa$vectors.cor[wh.tip, 1], 
         y = morph.pcoa$vectors.cor[wh.tip, 2], col = tip.col, pch = 16)

  phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 1:2], 
                             A = mode.pcoa$vectors.cor[node.seq, 1:2], 
                             control = con, label = "off", xlab = "PCO 1", 
                             ylab = "PCO 2", pch = NA)
  mtext(paste(colnames(taxon.bins)[t], "ecospace"))
  points(x = mode.pcoa$vectors.cor[wh.node, 1],
         y = mode.pcoa$vectors.cor[wh.node, 2], col = node.col, pch = 16)
  points(x = mode.pcoa$vectors.cor[wh.tip, 1], 
         y = mode.pcoa$vectors.cor[wh.tip, 2], col = tip.col, pch = 16)
  
  phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 3:4], 
                             A = morph.pcoa$vectors.cor[node.seq, 3:4], 
                             control = con, label = "off", xlab = "PCO 3", 
                             ylab = "PCO 4", pch = NA)
  mtext(paste(colnames(taxon.bins)[t], "morphospace"))
  points(x = morph.pcoa$vectors.cor[wh.node, 3],
         y = morph.pcoa$vectors.cor[wh.node, 4], col = node.col, pch = 16)
  points(x = morph.pcoa$vectors.cor[wh.tip, 3], 
         y = morph.pcoa$vectors.cor[wh.tip, 4], col = tip.col, pch = 16)
  
  phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 3:4], 
                             A = mode.pcoa$vectors.cor[node.seq, 3:4], 
                             control = con, label = "off", xlab = "PCO 3", 
                             ylab = "PCO 4", pch = NA)
  mtext(paste(colnames(taxon.bins)[t], "ecospace"))
  points(x = mode.pcoa$vectors.cor[wh.node, 3],
         y = mode.pcoa$vectors.cor[wh.node, 4], col = node.col, pch = 16)
  points(x = mode.pcoa$vectors.cor[wh.tip, 3], 
         y = mode.pcoa$vectors.cor[wh.tip, 4], col = tip.col, pch = 16)

  }
# dev.off()
par(op)





## PLOT PHYLOMORPHO/ECOSPACES FOR INDIVIDUAL CLASSES ###########################

# Set phylomorphospace/phyloecospace plotting variables (same as above)
# pdf(file = "Class_phylospaces.pdf")
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25), pty = "m")
branch.col <- "gray75"
tip.col <- "#0000FF7F"   # Set blue transparent so overlays as density
node.col <- "#A020F07F"  # Set purple transparent so overlays as density
cam.col <- "#87CEEB7F"   # For Camerata, overlay sky blue 
Tree <- morph.pcoa$Tree  # Same tree topology for both morphology and ecology
tip.seq <- 1:Ntip(Tree)
node.seq <- (Ntip(Tree) + 1):(Ntip(Tree) + Nnode(Tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(Tree$edge)), 
                                as.character(Tree$edge[, 2])))
# Process the classes (removing UNCERTAINs):
load("taxon.list")
# Replace Class Homostelea with (here synonymous) Cincta instead
taxon.list[which(taxon.list == "Homostelea")] <- "Cincta"
(class.tab <- sort(table(taxon.list), decreasing = TRUE))
classes <- names(class.tab)
classes <- classes[which(classes != "UNCERTAIN")]

# Camerata (Crinoidea) are plotted with separate color (but only for tips).
# Import following because lists order names
data <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", 
                 header = TRUE, stringsAsFactors = FALSE)
wh.cam <- which(data$Subclass == "Camerata")

ncl <- length(classes)
for(cl in 1:ncl) {
  wh.gr <- which(taxon.list == classes[cl])
  
  # Combine Diploporita with paraphyletic 'diploporitan'
  if (classes[cl] == "'diploporitan'")
    next
  if (classes[cl] == "Diploporita") {
    wh.gr <- which(taxon.list == "Diploporita" | taxon.list == "'diploporitan'")
  }
  
  wh.tip <- wh.gr[which(wh.gr <= Ntip(Tree))]
  wh.node <- wh.gr[which(wh.gr > Ntip(Tree))]
  
  phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 1:2], 
                             A = morph.pcoa$vectors.cor[node.seq, 1:2], 
                             control = con, label = "off", xlab = "PCO 1", 
                             ylab = "PCO 2", pch = NA)
  mtext(paste(classes[cl], "morphospace"))
  points(x = morph.pcoa$vectors.cor[wh.node, 1],
         y = morph.pcoa$vectors.cor[wh.node, 2], col = node.col, pch = 16)
  points(x = morph.pcoa$vectors.cor[wh.tip, 1], 
         y = morph.pcoa$vectors.cor[wh.tip, 2], col = tip.col, pch = 16)
  # Overlay camerate crinoids
  if (classes[cl] == "Crinoidea")
    points(x = morph.pcoa$vectors.cor[wh.cam, 1],
           y = morph.pcoa$vectors.cor[wh.cam, 2], col = cam.col, pch = 16)

  phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 1:2], 
                             A = mode.pcoa$vectors.cor[node.seq, 1:2], 
                             control = con, label = "off", xlab = "PCO 1", 
                             ylab = "PCO 2", pch = NA)
  mtext(paste(classes[cl], "ecospace"))
  points(x = mode.pcoa$vectors.cor[wh.node, 1],
         y = mode.pcoa$vectors.cor[wh.node, 2], col = node.col, pch = 16)
  points(x = mode.pcoa$vectors.cor[wh.tip, 1], 
         y = mode.pcoa$vectors.cor[wh.tip, 2], col = tip.col, pch = 16)
  if (classes[cl] == "Crinoidea")
    points(x = mode.pcoa$vectors.cor[wh.cam, 1],
           y = mode.pcoa$vectors.cor[wh.cam, 2], col = cam.col, pch = 16)
  
  phytools::phylomorphospace(tree = Tree, X = morph.pcoa$vectors.cor[tip.seq, 3:4], 
                             A = morph.pcoa$vectors.cor[node.seq, 3:4], 
                             control = con, label = "off", xlab = "PCO 3", 
                             ylab = "PCO 4", pch = NA)
  mtext(paste(classes[cl], "morphospace"))
  points(x = morph.pcoa$vectors.cor[wh.node, 3],
         y = morph.pcoa$vectors.cor[wh.node, 4], col = node.col, pch = 16)
  points(x = morph.pcoa$vectors.cor[wh.tip, 3], 
         y = morph.pcoa$vectors.cor[wh.tip, 4], col = tip.col, pch = 16)
  if (classes[cl] == "Crinoidea")
    points(x = morph.pcoa$vectors.cor[wh.cam, 3],
           y = morph.pcoa$vectors.cor[wh.cam, 4], col = cam.col, pch = 16)
  
  phytools::phylomorphospace(tree = Tree, X = mode.pcoa$vectors.cor[tip.seq, 3:4], 
                             A = mode.pcoa$vectors.cor[node.seq, 3:4], 
                             control = con, label = "off", xlab = "PCO 3", 
                             ylab = "PCO 4", pch = NA)
  mtext(paste(classes[cl], "ecospace"))
  points(x = mode.pcoa$vectors.cor[wh.node, 3],
         y = mode.pcoa$vectors.cor[wh.node, 4], col = node.col, pch = 16)
  points(x = mode.pcoa$vectors.cor[wh.tip, 3], 
         y = mode.pcoa$vectors.cor[wh.tip, 4], col = tip.col, pch = 16)
  if (classes[cl] == "Crinoidea")
    points(x = mode.pcoa$vectors.cor[wh.cam, 3],
           y = mode.pcoa$vectors.cor[wh.cam, 4], col = cam.col, pch = 16)
  
}
# dev.off()
par(op)










## COPHYLOGRAMS / TANGLEGRAMS ##################################################
# Tanglegrams ( = cophylograms) compare phylogenetic structure to dendrograms
# based solely on morphology or ecology.

# Import distance matrices from 3-DisparityDistances.R
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")

# Convert distance matrixes (restricting to tips) to cluster analysis
morph.hc <-
  hclust(as.dist(morph.distances.GED.5$distance_matrix[1:Ntip(tree), 1:Ntip(tree)]))
tr.morph <- as.phylo(morph.hc)
tr.morph$tip.label <-  tree$tip.label
plot(tr.morph)

# Cophyloplot between phylogeny and morpho-cluster
# 'association' not needed because strictly matches tip labels by default
cophylo.morph <- cophylo(tree, tr.morph, rotate = TRUE, rotate.multi = TRUE, 
                         print = TRUE)
beep(3)

# Confirm taxa are paired correctly
summary(cophylo.morph)

# save(cophylo.morph, file = "cophylo.morph")
# load("cophylo.morph")
plot(cophylo.morph)
# The phylogeny is the left figure

# plot.cophylo() allows following customizable plotting arguments:
# For the linkages: link.col, link.lwd, and link.lty; part controls width of 
#                   the linkages (but default 0.4 is best)
# For the cladogram/dendrogram: edge.col
# For the middle length: edge.col
# For font size: fsize
# For plotting points at tips: pts = T/F
# pdf(file = "MorphPhylogram.pdf")
edge.colors <- list(left = rep("darkgray", nrow(cophylo.morph$trees[[1]]$edge)), 
                    right = rep("black", nrow(cophylo.morph$trees[[2]]$edge)))
plot(cophylo.morph, edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
     fsize = 0.3, pts = FALSE)
# dev.off()

# How much sum-of-square rank differences between the two "phylogenies"? (same
# metric used in 'cophylo' to optimize the tip alignments)
(diffs.morph <- attr(cophylo.morph$trees[[1]], "minRotate"))
sqrt(diffs.morph / Ntip(tree)) # Average rank difference

# What is the distribution of rank differences between the two "phylogenies"?
phylo.order <- cophylo.morph$trees[[1]]$tip.label
morph.order <- cophylo.morph$trees[[2]]$tip.label
match.offset <- match(phylo.order, morph.order)

# Confirm matched correctly:
identical(phylo.order, morph.order[match.offset])

# Confirm sum of squares matches that in the cophylogram:
morph.offset <- match.offset - seq(phylo.order)
identical(diffs.morph, sum(morph.offset ^ 2))

hist(morph.offset, main = "rank mismatch between phylogeny and morphogram", 
     cex.main = .9)


# pdf(file = "EchinoMorphPhylogeny.pdf", height = 200)
# plot(cophylo.morph, edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
#      fsize = 0.65, pts = FALSE)
# dev.off()



## Do the same for the ecology tree (using mode data set)
# Convert distance matrixes (restricting to tips) to cluster analysis
eco.hc <-
  hclust(as.dist(mode.distances.GED.5$distance_matrix[1:Ntip(tree), 1:Ntip(tree)]))
tr.eco <- as.phylo(eco.hc)
tr.eco$tip.label <-  tree$tip.label
plot(tr.eco)

# Cophyloplot between phylogeny and eco-cluster
# 'association' not needed because strictly matches tip labels by default
cophylo.eco <- cophylo(tree, tr.eco, rotate = TRUE, rotate.multi = TRUE, 
                       print = TRUE)
beep(3)

# Confirm taxa are paired correctly
summary(cophylo.eco)

# save(cophylo.eco, file = "cophylo.eco")
# load("cophylo.eco")
plot(cophylo.eco)
# The phylogeny is the left figure
par(op)

# pdf(file = "EcoPhylogram.pdf")
edge.colors <- list(left = rep("darkgray", nrow(cophylo.eco$trees[[1]]$edge)), 
                    right = rep("blue", nrow(cophylo.eco$trees[[2]]$edge)))
plot(cophylo.eco, edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
     fsize = 0.3, pts = FALSE)
# dev.off()
par(op)

# How much sum-of-square rank differences between the two "phylogenies"? (same
# metric used in 'cophylo' to optimize the tip alignments)
(diffs.eco <- attr(cophylo.eco$trees[[1]], "minRotate"))
sqrt(diffs.eco / Ntip(tree)) # Average rank difference

# What is the distribution of rank differences between the two "phylogenies"?
phylo.order <- cophylo.eco$trees[[1]]$tip.label
eco.order <- cophylo.eco$trees[[2]]$tip.label
match.offset <- match(phylo.order, eco.order)

# Confirm matched correctly:
identical(phylo.order, eco.order[match.offset])

# Confirm sum of squares matches that in the cophylogram:
eco.offset <- match.offset - seq(phylo.order)
identical(diffs.eco, sum(eco.offset ^ 2))

hist(eco.offset, main = "rank mismatch between phylogeny and ecogram", 
     cex.main = .9)

# pdf(file = "EchinoEcoPhylogeny.pdf", height = 200)
# plot(cophylo.eco)
# dev.off()


## Are the offsets different in the two dendrograms?
summary(morph.offset)
sd(morph.offset)

summary(eco.offset)
sd(eco.offset)

# Because both distributions sum to zero, a t-test is not appropriate

# Non-parametric Mann-Whitney U-test to test whether different
wilcox.test(eco.offset, morph.offset)          # Not different in raw offsets
wilcox.test(eco.offset ^ 2, morph.offset ^ 2)  # But sig. diff in sum of square offsets
# CONCLUSION: They overlap significantly, but there are differences in the
# tails.

# Seems the major difference is the much longer tails of the ecological
# distribution. Are they different distributional shapes?
ks.test(eco.offset, morph.offset)
ks.test(eco.offset ^ 2, morph.offset ^ 2)
# CONCLUSION: Distributions are different (p < 0.001)

# Let's focus instead on a test of different variances using the two-sided
# F-test:
var.test(eco.offset, morph.offset)
# CONCLUSION: The ecological data set has substantially greater variance than
# the morphological data set, implying a greater degree of phylogenetic
# structure in the morphological data set (and less in the ecological data set).

boxplot(morph.offset, eco.offset, col = "gray", yaxs = "i",
        names = c("morphology", "ecology"), 
        main = "paired offsets between phylogeny and dendrogram")

# pdf(file = "cophylogram_offsets.pdf")
breaks <- seq(-350, 350, 25)
hist(c(eco.offset, morph.offset), main = "Offsets between phylogeny and dendrogram", 
     xlab = "paired offsets", ylab = "#", breaks = breaks, col = "transparent", 
     border = "transparent")
hist(eco.offset, add = T, border = "white", col = "darkgray", breaks = breaks)
hist(morph.offset, add = T, border = "black", col = "transparent", breaks = breaks)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
# dev.off()



