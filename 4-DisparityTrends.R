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
library(ape)        # v. 5.5.2
if(packageVersion("ape") < "5.5.2")
  stop("outdated version of 'ape'. Get updated version from GitHub or CRAN.\n")
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
  stop("outdated version of 'Claddis'. Get updated version from GitHub or CRAN.\n")

# Modification of geoscale::geoscalePlot to allow ICS 2020 timescale 
source("~/Manuscripts/CamOrdEchinos/geoscalePlot2.R")

# Modification of FD::dbFD() to simplify calculation of FD metrics, to allow
# standardization of PCoA vector coordinates by relative magnitude of
# eigenvalues, and to directly use pre-calculated distance matrices and PCoA
# output.
source("~/Manuscripts/CamOrdEchinos/simple.dbFD.R")

# Modification of ecospace::calc_metrics. See source file for modifications.
source("~/Manuscripts/CamOrdEchinos/calc_metrics2.R")

# Corrected resampling function
sample2 <- function(x, ...) x[sample.int(length(x), ...)]

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
series.boundaries <- c(521, 509, 497, 470, 458.4)
pd.boundaries <- c(541, 485.4, 443.8)
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
# load("taxon.list")

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
class.bins <- vector("list", 50)
for (t in 1:length(class.bins)){
  t.class.bins <- matrix(NA, nrow = length(unique.classes), ncol = nrow(ages))
  colnames(t.class.bins) <- ages$interval_name
  rownames(t.class.bins) <- unique.classes
  for (c in 1:nrow(t.class.bins)) {
    wh.class <- which(taxon.list[[t]][, "class"] == rownames(t.class.bins)[c])
    if (length(wh.class) == 1L)
      t.class.bins[c, ] <- as.numeric(taxon.bins[[t]][wh.class, ])
    else
      t.class.bins[c, ] <- apply(taxon.bins[[t]][wh.class, ], 2, sum)
  }
  class.bins[[t]] <- t.class.bins
}

# View output
head(class.bins[[50]])

# Save (and reload)
# save("class.bins", file = "class.bins")
# load("class.bins")

# Calculate median class.bins
mean.class.bins <- apply(simplify2array(class.bins), 1:2, median)
head(mean.class.bins)

# Plot diversity curve
cl.div <-
  sapply(1:ncol(mean.class.bins), function(x) sum(mean.class.bins[, x] > 0))
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
# stack.color <- rev(viridisLite::turbo(nrow(mean.class.bins)))
# Use next if prefer non-adjacents
set.seed(4); stack.color <- sample(rev(viridisLite::turbo(nrow(mean.class.bins))))
stackpoly(x = -mids, y = t(mean.class.bins[class.order, ]), col = stack.color,
          xlab = "time", ylab = "number of genera", stack = TRUE, 
          xlim = c(-541, -443.8), main = "genus richness (by class)")
abline(v = -series.boundaries, col = "white", lty = 5)
abline(v = -pd.boundaries, col = "white", lty = 1, lwd = 2)
box()
legend("topleft", legend = rev(unique.classes[class.order]), pch = 15, ncol = 3, 
       pt.cex = 1, cex = .55, col = rev(stack.color), bty = "n")
par(op)
# dev.off()



## TRIM MISSING TAXA FROM RAW DATA SET FOR SUBSEQUENT PROCESSING ###############

## Want to remove missing taxa? (Only used for 'raw' data treatment)

## Function to append the numbered index and original names of removed tip and
## node numbers when running Claddis::trim_matrix().
#
# Claddis::trim_matrix() alters the names [= numbers] of nodes in a tree when
# pruning out tips and nodes with missing data. This version wraps around
# trim_matrix() and appends the original node numbers, so that the trimmed
# output better matches the other same-named objects in the workflow. Note that
# the index order of tip names corresponds to the order of 'tip.label' in the
# provided tree, and not necessarily the same order as in other objects used in
# the workflow later. Note also that the length of $removed_taxa_index may be
# shorter than the length of $removed_taxa because of the iterative and
# multi-step manner that trim_matrix identifies problematic taxa, which can
# remove tip taxa multiple times.
#
#      distance_matrix = A distance matrix in the format created by
#                        Claddis::calculate_morphological_distances(), or the 
#                        code used previously in the workflow.
#      tree            = a tree. Used to match tip labels and node numbers to 
#                        trimmed taxa, so must be included, unlike behavior in
#                        trim_matrix().
trim_to_numbers <- function(distance_matrix, tree = NULL) {
  if (is.null(tree))
    stop("'tree' MUST be included.")
  trimmed <- Claddis::trim_matrix(distance_matrix, tree)
  # Get list of tips
  tips.removed <- trimmed$removed_taxa[-grep("%%", trimmed$removed_taxa)]
  tip.numbers.removed <- match(tips.removed, tree$tip.label)
  # Get list of nodes removed (in original 'named' output)
  nodes.removed <- trimmed$removed_taxa[grep("%%", trimmed$removed_taxa)]
  sq <- 1:length(nodes.removed)
  # ... and convert to original node numbers
  node.numbers.removed <- sort(unique(unlist(lapply(sq, function(sq)
    Claddis::find_mrca(unlist(strsplit(trimmed$removed_taxa[which(trimmed$removed_taxa == nodes.removed[sq])], "%%")), tree)))))
  # Confirm index is logical (tip numbers <= no. of tree times, and nodes >)
  if (max(tip.numbers.removed) > ape::Ntip(tree))
    stop("returns tip numbers greater than allowed by 'tree'")
  if (min(node.numbers.removed) <= ape::Ntip(tree))
    stop("returns node numbers less than allowed by 'tree'")
  if (any(node.numbers.removed %in% tip.numbers.removed))
    stop("node numbers can not match tip numbers")
  # Append to trim_matrix output
  trimmed$removed_taxa_index <- unique(c(tip.numbers.removed, node.numbers.removed))
  trimmed$removed_taxa_by_name <-
    c(sort(unique(tree$tip.label[tip.numbers.removed])),
      as.character(sort(unique(node.numbers.removed))))
  return(trimmed)
}

# Make sure to use the original time-scaled tree instead of the trimmed tree for
# subsequent steps. In subsequent steps, make sure to match the original names
# with saved objects in case the order has changed at some point.
load("raw.anc")
load("raw.distances.GED.5")
ntrees <- length(raw.anc)
raw.trimmed <- vector("list", ntrees)
(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
raw.trimmed <- foreach(i = 1:ntrees, .inorder = TRUE) %dopar% { 
  tree.trimmed <-
    trim_to_numbers(distance_matrix = raw.distances.GED.5[[i]]$distance_matrix,
                    tree = raw.anc[[i]]$topper$tree)
}
stopCluster(cl)
(Sys.time() - start) # 2.6 minutes on 8-core laptop
# Save (and reload)
# save(raw.trimmed, file = "raw.trimmed")
beep(3)

# Compare how many taxa removed from the trees?
sq <- 1:ntrees
table(unlist(lapply(sq, function(sq) raw.trimmed[[sq]]$removed_taxa_by_name)))
table(unlist(lapply(sq, function(sq) length(raw.trimmed[[sq]]$removed_taxa_index))))




## CALCULATE PCoA ORDINATIONS AND MORPHO/ECOSPACES #############################

# Given high correlations among distance matrices, using GED.5 as the distance
# metric. (It is Will's 1988 generalized Euclidean distance, which replaces
# missing character states using the mean pairwise dissimilarity for each
# pairwise comparison. Hierarchical character dependencies are weighted using
# Hopkins and St. John 2018 alpha = 0.5, which in which the primary character
# contributes weight according to the fraction of shared secondary [and
# tertiary, etc.] characters.

## Pick desired distance matrix (and corresponding ancestral state
## reconstruction objects, for the data-specific tree object)
# dist.matrix <- mode.distances.GED.5; anc <- mode.anc
# dist.matrix <- constant.distances.GED.5; anc <- constant.anc
dist.matrix <- raw.distances.GED.5; anc <- raw.anc
# dist.matrix <- morph.distances.GED.5; anc <- morph.anc

# View sample
dist.matrix[[50]]$distance_matrix[1:10, 1:4]

# Observe data structure of the data matrix:
hist(c(dist.matrix[[50]]$distance_matrix))
summary(c(dist.matrix[[50]]$distance_matrix))

# Any missing distances? (Should only be TRUE for pre-trimmed raw treatment, and
# some also occur in trees 31-35 in the morphological data set)
sq <- 1:50
table(sapply(sq, function(sq) anyNA(dist.matrix[[sq]]$distance_matrix)))

# Any non-diagonal zeros (all are, which is potentially problematic for some
# ordination analyses)
table(sapply(sq, function(sq)
  (length(dist.matrix[[sq]]$distance_matrix == 0) - nrow(dist.matrix[[sq]]$distance_matrix)) > 0))

# Are they Euclidean? (all are non-Euclidean, so recommend using the Cailliez or
# Lingoes correction for negative eigenvalues). (Does not compute for
# morphological trees 31-35 nor any 'raw' tree)
table(sapply(sq, function(sq)
  ade4::is.euclid(as.dist(dist.matrix[[sq]]$distance_matrix))))



## Perform a phylogenetic Principal Coordinates Analysis:

# Skip 'Claddis::ordinate_cladistic_matrix' [which bundles ancestral state
# reconstruction and calculation of distance matrices] by using direct
# calculation of PCoA so can directly import in previously calculated distance
# matrices.

# For consistency with other Claddis functions, using ape::pcoa() instead of
# 'ecospace' and 'FD's use of ade4::dudi.pco(), which produces identical output
# (in different formats), and are essentially equally fast. Variation typically
# has to do with how many non-negative eigenvectors are returned, and with
# arbitrary reversing of axes. Eigenvectors are identical for at least first
# 11-15 axes (not demonstrated here, see "X_Diff pcoa trials.R")

# Because many of the distance matrices are non-Euclidean (with non-diagonal
# zeroes) and return negative eigenvalues, we are using the "Lingoes" correction
# here and below for consistency in how all morphospace/ecospace disparity
# metrics are calculated. Tests on earlier trials (not demonstrated here, see
# "X_Diff pcoa trials.R") found little difference between using the Cailliez and
# Lingoes corrections, but the Lingoes correction is mathematically more
# tractable with a greater number of the trees and data sets used here.
# Ultimately, the impact of corrections is minor because comparison of the
# character-space (morphospace/ecospace) structure for the first 6 axes are
# essentially identical across corrections (for those distance matrices that can
# be calculated for them). Second, the use of all eigenvalues is only required
# for the FEve functional diversity metric; all others are calculated directly
# on the distance matrix or on a reduced-dimensional PCoA eigenvalue space, and
# so there is negigible impact on the choise of correction used.

# For 'raw' data treatment, need to use previously built 'raw.trimmed' to remove
# incalculable taxa. Recall the list of names was drawn from
# raw.anc[[x]]$topper$tree$tip.label instead of raw.trimmed[[x]]$tree because
# goal is to identify original, trimmed tips and nodes present in other
# objects.)
load("raw.trimmed")
load("raw.anc")
for (t in 1:length(raw.trimmed)) {
  taxa.to.cut <- raw.trimmed[[t]]$removed_taxa_by_name
  wh.to.cut.dist <-
    match(taxa.to.cut, rownames(dist.matrix[[t]]$distance_matrix))
  dist.matrix[[t]]$distance_matrix <-
    dist.matrix[[t]]$distance_matrix[-wh.to.cut.dist, -wh.to.cut.dist]
  dist.matrix[[t]]$tree <- raw.trimmed[[t]]$tree
  dist.matrix[[t]]$removed_taxa <- raw.trimmed[[t]]$removed_taxa
}
# View sample (note NAs now mostly trimmed out)
dist.matrix[[50]]$distance_matrix[1:10, 1:4]



# Running as a loop to enhance diagnostic abilities because of the idiosyncratic
# nature of pcoa calculations (usually involving presence of negative
# eigenvalues and/or non-diagonal zeroes in distance matrices).
(start <- Sys.time())
pcoa.results <- vector("list", length(dist.matrix))
index <- c(1, seq(5, 50, by = 5))
for(t in 1:length(pcoa.results)) {
# Use for morphological data set to restart when hits problematic trees 31-35
# for(t in 36:length(pcoa.results)) { 
    if (t %in% index)
    cat("processing tree", t, "\n")
  t.pcoa.results <- ape::pcoa(dist.matrix[[t]]$distance_matrix,
                              correction = "lingoes",
                              rn = rownames(dist.matrix[[t]]$distance_matrix))
  # Append tree, required for downstream Claddis functions
  if (is.null(t.pcoa.results$tree)) t.pcoa.results$tree <- anc[[t]]$topper$tree
  # Only use next line instead for the 'raw' treatment (where trimmed tree was
  # added in pre-processing):
  # if (is.null(t.pcoa.results$tree)) t.pcoa.results$tree <- dist.matrix[[t]]$tree
  pcoa.results[[t]] <- t.pcoa.results
}
(Sys.time() - start) #~ 1-10 minutes using 1 core

# Save (and re-load) PCoA output:
# mode.pcoa <- pcoa.results; save(mode.pcoa, file = "mode.pcoa")
# constant.pcoa <- pcoa.results; save(constant.pcoa, file = "constant.pcoa")
# raw.pcoa <- pcoa.results; save(raw.pcoa, file = "raw.pcoa")
# morph.pcoa <- pcoa.results; save(morph.pcoa, file = "morph.pcoa")
# load("morph.pcoa"); pcoa.results <- morph.pcoa
# load("mode.pcoa"); pcoa.results <- mode.pcoa
# load("constant.pcoa"); pcoa.results <- constant.pcoa
# load("raw.pcoa"); pcoa.results <- raw.pcoa
beep(3)

# Using new pcoa() with Cailliez correction, no errors for Mode or constant or
# morph (except trees 31-35)! but Cailliez errors for raw. Seems to be avoided
# when using Lingoes correction.

# No errors (and much faster!, <<6 minutes) in raw and morph when using the
# Lingoes correction!!!

# Plot sample scree plot and biplot
barplot(100 * pcoa.results[[50]]$values$Rel_corr_eig[1:15], main = "first 10 axes", names.arg = 1:15,
        xlab = "Axes", ylab = "relative % explained")
barplot(100 * cumsum(pcoa.results[[50]]$values$Rel_corr_eig)[1:15], names.arg = 1:15,
        xlab = "Axes", ylab = "relative % explained", 
        main = "cum var explained for first 10 axes")
plot(pcoa.results[[50]]$vectors.cor[, 1], pcoa.results[[50]]$vectors.cor[, 2], xlab = "PCoA axis 1", 
     ylab = "PCoA axis 2")

round(100 * pcoa.results[[50]]$values$Rel_corr_eig[1:10], 2)
round(100 * cumsum(pcoa.results[[50]]$values$Rel_corr_eig)[1:30], 1)


# RESULTS FOR TREE 50 (% explained; cumulative % explained):

# Morph:    1- 1.68%, 2-1.02% (2.7%), 3-0.51% (3.2%), 4-0.42% (3.7%), 
#                     5-0.28% (3.9%), 6-0.23% (4.2%),   
# Mode:     1- 1.42%, 2-0.38% (1.8%), 3-0.33% (2.1%), 4-0.26% (2.4%), 
#                     5-0.22% (2.6%), 6-0.20% (2.8%),   
# Constant: 1- 0.81%, 2-0.29% (1.1%), 3-0.25% (1.3%), 4-0.20% (1.6%), 
#                     5-0.19% (1.7%), 6-0.17% (1.9%),   
# Raw:      1- 0.66%, 2-0.39% (1.1%), 3-0.30% (1.4%), 4-0.27% (1.6%), 
#                     5-0.23% (1.9%), 6-0.23% (2.1%),   

# Plot all scree plots on one figure
# pdf(file = "PCoAScreePlots.pdf")
load("morph.pcoa")
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")
par(mfrow = c(2,2), mar = c(4, 4, 2, .1))
barplot(100 * morph.pcoa[[50]]$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "morphology", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
barplot(100 * mode.pcoa[[50]]$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (mode)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
# note no correction used, but pcoa() still calculates 'corrected' eigenvalues
barplot(100 * constant.pcoa[[50]]$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (constant)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
barplot(100 * raw.pcoa[[50]]$values$Rel_corr_eig[1:10], names.arg = 1:10, 
        main = "ecology (raw)", xlab = "PCoA eigenvalues", 
        ylab = "relative % explained", cex.names = 1)
par(op)
# dev.off()



## Calculate (mock) "factor loadings" (for tree #50)
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

loadings.morph <- mock.loadings(orig.vars = morph.anc[[50]]$matrix_1$matrix,
                                ord.coord = morph.pcoa[[50]]$vectors.cor, vars = 6,
                                cutoff = 0.8)
na.omit(loadings.morph)

loadings.mode <- mock.loadings(orig.vars = mode.anc[[50]]$matrix_1$matrix, 
                               ord.coord = mode.pcoa[[50]]$vectors.cor, vars = 6, 
                               cutoff = 0.4)
na.omit(loadings.mode)

loadings.constant <- mock.loadings(orig.vars = constant.anc[[50]]$matrix_1$matrix,
                                   ord.coord = constant.pcoa[[50]]$vectors.cor, vars = 6,
                                   cutoff = 0.4)
na.omit(loadings.constant)

# For the raw treatment, need to remove genera with missing data excluded when
# trimmed
load("raw.trimmed")
taxa.to.cut <- raw.trimmed[[50]]$removed_taxa_by_name
wh.to.cut.anc <- match(taxa.to.cut, rownames(raw.anc[[50]]$matrix_1$matrix))
raw.anc[[50]]$matrix_1$matrix <- raw.anc[[50]]$matrix_1$matrix[-wh.to.cut.anc, ]
loadings.raw <- mock.loadings(orig.vars = raw.anc[[50]]$matrix_1$matrix,
                              ord.coord = raw.pcoa[[50]]$vectors.cor, vars = 6, 
                              cutoff = 0.3)
na.omit(loadings.raw)

# RESULTS (*** UPDATED USING ONLY TREE #50 ***):
# MORPHOLOGY DATA SET:     -                               +
#  PCO 1:     4, 134, 141                             46, 227, 325
#  PCO 2:     41, 46, 126, 140, 143, 181, 274         31, 129, 160, 173, 339
#  PCO 3:     126, 134, 141, 172, 227                409
#  PCO 4:     31, 129                                 41, 46, 126, 140, 143, 
#  PCO 5:     38, 46, 170, 380                        32, 42, 126, 135, 173
#  PCO 6:     31, 129                                 41, 46, 121, 140, 143

# MODE DATA SET:      -                                +
#  PCO 1:     6, 13, 16, (24, 26, 31, 36)            2-5, (9), 12, 15, (23, 25, 28-29)
#   +: Attached (filter feeders) living far from hard (or biotic) substrate
#   -: Free-living, mobile (mass-feeding deposit/carnivore/scavenging feeders) living on soft substrates)
#  PCO 2:   (28)                                       1, (24, 26, 31, 36)
#   +: Large bodied echinoderms (especially mass/deposit/carnivore/scavenging feeders)
#   -: Filter feeders
#  PCO 3:    (1)                                     (29)
#   +: (Smallest with) highest density filter feeders
#   -: Largest (non-filterers?)
#  PCO 4:    (3)                                      (9, 12)
#   +: (Biotic and hard substrates) 
#   -: (Largest RelStrats)
#  PCO 5:   (10, 24, 26, 32, 40)                     (23, 25)
#   +: (Food above substrate)
#   -: (Raptorial/bulk feeders [= scavengers] on lithic substrates with infaunal food)
#  PCO 6:    (5)
#   +: Small RelFoodStrat (= eating in contact with food)
#   -: Large RelFoodStrat (= food at a distance)

# CONSTANT DATA SET:  -                                +
#  PCO 1:        6, 13, 16, (31)                  2-5, 12, 15, (28-29)
#  PCO 2:      (23, 25, 28, 39)                     1 (24, 26, 31, 36, 40)
#  PCO 3:      (12)                               (2-4, 13), 29     
#  PCO 4:       (5, 13, 23, 25, 28, 34, 39)       (12), 24, 26, (31-32, 36, 40)
#  PCO 5:       (9, 29)                           (10)
#  PCO 6:       (6)                               (13, 29)

#    OVERALL: PCOs basically same as for mode, with 4 = mode's #5 reversed. #5
#    and 6 emphasize slightly different habits (#5 for lithic / biotic & dense
#    filtering, and #6 for mobile vs. dense filterers on soft-substrate), but
#    with very low loading scores.

#  PCO 1:        6, 13, 16, (24, 26, 31)          1-5, 12, 15, (28-29)
#  PCO 2:      (23, 25, 28, 39)                     1, (6), 24, 26, 31, (32, 36)
#  PCO 3:     (2-4, 24, 26, 29-32, 36, 40)        (23, 25, 28)
#  PCO 4:       (6, 13, 16, 29)                    12, (15, 36, 40)
#  PCO 5:      (10, 31)                            (9, 28-29)                      
#  PCO 6:      (12)                               (13)

#    OVERALL: PCOs basically same as for mode (but more like constant), with #5
#    reversed from constant.


# Compare before and after standardizing (for tree #50)
new <- pcoa.results[[50]]
new$vectors.cor <- stand.pcoa(vectors = pcoa.results[[50]]$vectors.cor,
                          eigenvalues = pcoa.results[[50]]$values$Corr_eig)
pcoa.results[[50]]$vectors.cor[1:5, 1:5]
new$vectors.cor[1:5, 1:5]
plot(pcoa.results[[50]]$vectors.cor[, 1], pcoa.results[[50]]$vectors.cor[, 2])
plot(new$vectors.cor[, 1], new$vectors.cor[, 2])

# CONCLUSION: The axes are rescaled correctly. (Note that doing so affects
# [i.e., stabilizes across time intervals] the values returned for FRic
# [especially so], FEve [minimally so], and FDiv [minimally so].)








## PLOT ORDINATIONS USING CLADDIS FUNCTIONS ####################################

## *** ONLY USING TREE #50 ***

# Plot a 2-dimensional morphospace with tips, nodes, and root:
plot_morphospace(pcoa.results[[50]])
tip.seq <- 1:Ntip(pcoa.results[[50]]$tree)
node.seq <-
  (Ntip(pcoa.results[[50]]$tree) + 1):(Ntip(pcoa.results[[50]]$tree) + Nnode(pcoa.results[[50]]$tree))
points(x = pcoa.results[[50]]$vectors.cor[node.seq, 1],
       y = pcoa.results[[50]]$vectors.cor[node.seq, 2], col = "gray", pch = 16) # nodes
points(x = pcoa.results[[50]]$vectors.cor[tip.seq, 1], 
       y = pcoa.results[[50]]$vectors.cor[tip.seq, 2])                          # tips
points(x = pcoa.results[[50]]$vectors.cor[max(node.seq), 1], 
       y = pcoa.results[[50]]$vectors.cor[max(node.seq), 2], col = "red", pch = 16, 
       cex = 1.5)                                                     # root





## DISPARITY / FUNCTIONAL DIVERSITY TRENDS #####################################

# Import ancestral states from 2-InferancestralStates.R
load("morph.anc")
load("mode.anc")
load("constant.anc")
load("raw.anc")

# Import distance matrices from 3-DisparityDistances.R
load("morph.distances.GED.5")
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")

# Load PCoA output from ape::pcoa
load("morph.pcoa")
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")

# Choose ancestral states, Wills GED-0.5 distance matrices, and PCoA output
# anc <- morph.anc; dist.matrix <- morph.distances.GED.5; pcoa.results <- morph.pcoa
# anc <- mode.anc; dist.matrix <- mode.distances.GED.5; pcoa.results <- mode.pcoa
# anc <- constant.anc; dist.matrix <- constant.distances.GED.5; pcoa.results <- constant.pcoa

# Load taxon.bins (built above)
load("taxon.bins")
if(nrow(taxon.bins[[1]]) != nrow(dist.matrix[[1]]$distance_matrix))
  stop("Need to re-load or re-build taxon.bins because some taxa were removed when running the 'raw' treatment.\n")

## For 'raw' treatment (ONLY!), need to remove the taxa trimmed out for the PCOA
## because of high density of missing states.
# Use 'raw.trimmed' output to override the raw objects. (Not needed for
# pcoa.results because used earlier trimmed version when building the PCoA
# object.) Because the tip names are based on the raw.anc[[x]]$topper$tree, need
# to match names for rows to correctly identify the trimmed-out tips and nodes
load("raw.trimmed")
for (t in 1:length(raw.trimmed)) {
  taxa.to.cut <- raw.trimmed[[t]]$removed_taxa_by_name
  # Because rownames are identical across 'dist.matrix', 'raw.anc', and
  # 'taxon.bins', using 'taxon.bins' matches here
  wh.to.cut <- match(taxa.to.cut, rownames(taxon.bins[[50]]))
  anc[[t]]$matrix_1$matrix <- raw.anc[[t]]$matrix_1$matrix[-wh.to.cut, ]
  dist.matrix[[t]]$distance_matrix <-
    raw.distances.GED.5[[t]]$distance_matrix[-wh.to.cut, -wh.to.cut]
  taxon.bins[[t]] <- taxon.bins[[t]][-wh.to.cut, ]
}

# Confirm all these 4 objects have the same number of rows
dim(anc[[50]]$matrix_1$matrix)
dim(taxon.bins[[50]])
dim(dist.matrix[[50]]$distance_matrix)
dim(pcoa.results[[50]]$vectors.cor)

# View to confirm (and in case of 'raw', that taxa with numerous NAs are trimmed
# out)
taxon.bins[[50]][1:5, 1:5]
anc[[50]]$matrix_1$matrix[1:10, 1:10]
dist.matrix[[50]]$distance_matrix[1:4, 1:4]
pcoa.results[[50]]$vectors.cor[1:4, 1:4]


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
#  R  =  sum of ranges (using 'pcoa' eigenvectors)
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
ntrees <- length(anc)
(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
# Use load-balancing because of different run times for the MCMC optimizations
opts <- list(preschedule = FALSE)
clusterSetRNGStream(cl, 3142)  # Set L-Ecuyer RNG seed
(start <- Sys.time())
# Run each tree in parallel
metrics <- vector("list", 50)
metrics <- foreach(i = 1:ntrees, .options.snow = opts, .inorder = TRUE, 
                   .packages = "FD") %dopar% {
  nc <- ncol(taxon.bins[[i]])
  t.metrics <- data.frame(Age = as.numeric(mids), S = NA, H = NA, D = NA, 
                          M = NA, V = NA, R = NA, FRic = NA, FEve = NA, 
                          FDiv = NA, FDis = NA, qual.FRic = NA)
  for (t in 1:nc) {
    wh.gr <- unname(which(taxon.bins[[i]][, t]))
    S <- length(wh.gr)
    anc.sam <- anc[[i]]$matrix_1$matrix[wh.gr,]
    dist.anc.sam <- dist.matrix[[i]]$distance_matrix[wh.gr, wh.gr]
    pcoa.sam <- pcoa.results[[i]]
    pcoa.sam$vectors.cor <- pcoa.results[[i]]$vectors.cor[wh.gr, ]
    H <- nrow(unique(dist.anc.sam))
    if (any(is.nan(dist.anc.sam)) | length(dist.anc.sam) == 0) next
    if (S <= m | H <= m) next
    FD <- rep(NA, 11)
    FD <- calc_metrics2(sample = anc.sam, dist.sam = dist.anc.sam, 
                        pcoa = pcoa.sam, m = m, stand.pcoa = TRUE, 
                        calc.FRic.and.FDiv = TRUE)
    t.metrics[t, 1 + seq.int(FD)] <- FD
  }
  return(t.metrics)
}
stopCluster(cl)
(Sys.time() - start) # 4.9-6.8 minutes on 8-core laptop

## Save / reload metrics
# morph.metrics <- metrics; save(morph.metrics, file = "morph.metrics")
# mode.metrics <- metrics; save(mode.metrics, file = "mode.metrics")
# constant.metrics <- metrics; save(constant.metrics, file = "constant.metrics")
# raw.metrics <- metrics; save(raw.metrics, file = "raw.metrics")
beep(3)
# load("morph.metrics"); metrics <- morph.metrics
# load("mode.metrics"); metrics <- mode.metrics
# load("constant.metrics"); metrics <- constant.metrics
# load("raw.metrics"); metrics <- raw.metrics

# Save tree # 50 to .csv
# write.csv(metrics[[50]], file = "metrics_morph.csv", row.names = FALSE)
# write.csv(metrics[[50]], file = "metrics_LH_mode.csv", row.names = FALSE)
# write.csv(metrics[[50]], file = "metrics_LH_constant.csv", row.names = FALSE)
# write.csv(metrics[[50]], file = "metrics_LH_raw.csv", row.names = FALSE)


# Check for any tree intervals skipped to avert error (ignore t = 18 = Ediacaran
# b/c generally skipped due to too low S or H):
sq <- 1:50
sapply(sq, function(sq) which(is.na(metrics[[sq]]$FRic)))

## Error log (usually qhull errors for FRic and FDiv). Following trees unable to
## calculate Terreneuvian FRic/FDiv:
# Morph:    No errors!
# Mode:     Trees 12, 21, 26, and 29 for Terreneuvian
# Constant: Trees 1-3, 5-6, 8-10, 12-13, 18, 23-24, 28-31, 33, 41, 47, and 50 for 
#           varying intervals, usually only 1-2 per tree (2-Katian, 3-Sandbian, 
#           4-Darriwilian, 5-Dapingian, 6-Floian, 7-Tremadocian, 8-Stage 10, 
#           12-Drumian, 13-Stage 5, 14-Stage 4, 15-Stage 3, and 16-Stage 2)
# Raw:      All trees, except #22, 26, 30, and 47, yield FRic/FDiv errors for 
#           several time intervals.

# What proportion of PCoA eigenvalues included in ordination-based metrics?
# (Note that because calculated on the entire PCoA eigenvectors, will not change
# when dealing with subsets.)
summary(unlist(lapply(sq, function(sq) metrics[[sq]]$qual.FRic)))

head(metrics[[50]])

# Sample relationships between metrics and species richness (for tree #50); see
# Novack-Gottshall (2016b) for why informative
par(mfrow = c(2, 4), mar = c(4, 4, 1, 0.2))
for (c in c(3:6, 8:11)) {
  if (sum(is.na(metrics[[50]][, c])) == length(metrics[[50]][, c]))
    plot( 1, type = "n", axes = FALSE, xlab = "", ylab = "")
  else
    plot(metrics[[50]]$S, metrics[[50]][, c], main = colnames(metrics[[50]])[c], 
         ylab = colnames(metrics[[50]])[c], xlab = "S")
}
par(op)


# Plot trends (using tree #50 as an example)
par(mar = c(0, 4, 2, 2))
for (c in 2:11) {
  var <- metrics[[50]][, c]
  if (sum(is.na(var)) == length(var))
    next
  lim <- range(var, na.rm = TRUE)
  geoscalePlot2(mids, rep(lim[1], length(mids)), units = c("Epoch", "Period"), 
               tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
               cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
               ts.col = TRUE, label = colnames(metrics[[50]])[c], timescale = ICS2020, 
               type = "n", abbrev = "Period")
  mtext(text = colnames(metrics[[50]])[c], side = 3, cex = 1.25)
  lines(mids, var, lwd = 3)
}
par(op)




## How correlated are the metrics? (Using just results from tree #50)
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
round(diag(cor(diff.morph[, 3:11], diff.mode[, 3:11], use = "pairwise.complete.obs")), 4)
summary(lm(diff.morph$H ~ diff.mode$H)) # r2 = 0.89 ***

# morphological vs. LH-constant
round(diag(cor(diff.morph[, 3:11], diff.constant[, 3:11], use = "pairwise.complete.obs")), 4)
summary(lm(diff.morph$H ~ diff.constant$H)) # r2 = 0.92 ***

# morphological vs. LH-raw
round(diag(cor(diff.morph[, 3:11], diff.raw[, 3:11], use = "pairwise.complete.obs")), 4)
summary(lm(diff.morph$H ~ diff.raw$H)) # r2 = 0.93 ***

# LH-mode vs. LH-constant
round(diag(cor(diff.mode[, 3:11], diff.constant[, 3:11], use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-raw
round(diag(cor(diff.mode[, 3:11], diff.raw[, 3:11], use = "pairwise.complete.obs")), 4)

# LH-constant vs. LH-raw
round(diag(cor(diff.constant[, 3:11], diff.raw[, 3:11], use = "pairwise.complete.obs")), 4)





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

# Subsampling samples both across time intervals and across trees simultaneously
# to account for variation in the tree structure, using 3000 total replicates
# (with 50 trees, drawing standard number of taxa per time interval, replicated
# 60 times per tree).

# See numreps_confirmation.xlsx for confirmation for why 3,000 replicates is
# sufficient to get relative error (CV = sd/mean) within 0.1% for all time
# intervals and (phylum-wide) sample sizes. Tested using the two time intervals
# with largest (= Floian) and smallest (= Paibian) species richness, using the
# fast-to-calculate D metric. (As expected, intervals with greater richness have
# greater sampling variability, which decreases with smaller resampling pools
# and greater numbers of replicates.)


# Load PCoA output from ape::pcoa
load("morph.pcoa")
load("mode.pcoa")
load("constant.pcoa")
load("raw.pcoa")

# Import ancestral states from 2-InferAncestralStates.R
load("morph.anc")
load("mode.anc")
load("constant.anc")
load("raw.anc")

# Import distance matrices from 3-DisparityDistances.R
load("morph.distances.GED.5")
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")

# Choose ancestral states, Wills GED-0.5 distance matrices, and PCoA output
# anc <- mode.anc; dist.matrix <- mode.distances.GED.5; pcoa.results <- mode.pcoa
# anc <- constant.anc; dist.matrix <- constant.distances.GED.5; pcoa.results <- constant.pcoa
# anc <- raw.anc; dist.matrix <- raw.distances.GED.5; pcoa.results <- raw.pcoa
# anc <- morph.anc; dist.matrix <- morph.distances.GED.5; pcoa.results <- morph.pcoa

# Load taxon.bins (built above)
load("taxon.bins")
if(nrow(taxon.bins[[50]]) != nrow(dist.matrix[[50]]$distance_matrix))
  stop("Need to re-load or re-build taxon.bins because some taxa were removed when running the 'raw' treatment.\n")

## For 'raw' treatment (ONLY!), need to remove the taxa trimmed out for the PCOA
## because of high density of missing states.
# Use 'raw.trimmed' output to override the raw objects. (Not needed for
# pcoa.results because used earlier trimmed version when building the PCoA
# object.) Because the tip names are based on the raw.anc[[x]]$topper$tree, need
# to match names for rows to correctly identify the trimmed-out tips and nodes
load("raw.trimmed")
for (t in 1:length(raw.trimmed)) {
  taxa.to.cut <- raw.trimmed[[t]]$removed_taxa_by_name
  # Because rownames are identical across 'dist.matrix', 'raw.anc', and
  # 'taxon.bins', using 'taxon.bins' matches here
  wh.to.cut <- match(taxa.to.cut, rownames(taxon.bins[[50]]))
  anc[[t]]$matrix_1$matrix <- raw.anc[[t]]$matrix_1$matrix[-wh.to.cut, ]
  dist.matrix[[t]]$distance_matrix <-
    raw.distances.GED.5[[t]]$distance_matrix[-wh.to.cut, -wh.to.cut]
  taxon.bins[[t]] <- taxon.bins[[t]][-wh.to.cut, ]
}

# Confirm all these 4 objects have the same number of rows
dim(anc[[50]]$matrix_1$matrix)
dim(dist.matrix[[50]]$distance_matrix)
dim(pcoa.results[[50]]$vectors.cor)
dim(taxon.bins[[50]])

# View to confirm (and in case of 'raw', that taxa with numerous NAs are trimmed
# out)
taxon.bins[[50]][1:5, 1:5]
anc[[50]]$matrix_1$matrix[1:10, 1:10]
dist.matrix[[50]]$distance_matrix[1:4, 1:4]
pcoa.results[[50]]$vectors.cor[1:4, 1:4]

# How many reps per tree per time interval? (60 provides sufficient control of
# variability per tree-bin)
nreps <- 60
cat("There will be", nreps * length(anc), 
    "total replicates for each time bin, sampling equally across", length(anc), "trees\n")

# Based on sensitivity tests: Higher 'std.g' more important than increasing
# 'nreps.' 60, 90, 100, and 120 nreps all ~ same (i.e., minimally increased
# precision above 60).

# Create index for pre-allocating position in trees and replicate number)
tree.seq <- rep(seq.int(anc), nreps)

# How many dimensions in PCoA for FRic and FDiv? Use 6 for phylum-level
m <- 6

# How many time bins? (Make sure taxon.bin[[1]] returns 18 bins! If 17, choose
# different taxon.bin tree)
cat((nc <- ncol(taxon.bins[[1]])), "time bins used (if sample sizes allow)\n")

# How many genera to sample per interval? (Use 50 for echinoderm-wide; allows
# all 50 trees to be used for all but early Furongian [Paibian], Terreneuvian
# [Fortunian and Stage 2], and Ediacaran. Mid Furongian [Jiangshanian] can be
# included if wish to use 92% of the trees. And technically, disparity for Stage
# 2 and Paibian are calculated, but using only 4% and 26% of trees.)
std.g <- 50

# Use TRUE if want to first check sample sizes to obtain std.g
check.std.g <- FALSE

# If wanting to confirm appropriate sampling quota
if (check.std.g) {
  for (b in 1:nc) {
    sq <- seq.int(tree.seq)
    bin.richness <- 
      unlist(lapply(sq, function(sq) length(which(taxon.bins[[tree.seq[sq]]][, b]))))
    cat("bin richness", b, "across trees: min = ", min(bin.richness), 
        ", median = ", median(bin.richness), ", max = ", max(bin.richness), 
        ", using", length(bin.richness[bin.richness >= std.g]) / 60, "trees\n")
  }
}
# 29 is lowest median richness across time bins (excl. Ediacaran, and excl. the
# 'raw' treatment, which also excludes bins 17-Fortunian, 16-Stage 2, and
# 10-Paibian)


## Loop through each time interval, running the sample-standardization
## replicates in parallel

# Create data frame to store means and SDs for each bin (across trees/replicates)
metrics <- data.frame(Age = as.numeric(mids), S = std.g, H = NA, SE.H = NA,
                      D = NA, SE.D = NA, M = NA, SE.M = NA, V = NA, SE.V = NA, 
                      R = NA, SE.R = NA, FRic = NA, SE.FRic = NA, FEve = NA, 
                      SE.FEve = NA, FDiv = NA, SE.FDiv = NA, FDis = NA, 
                      SE.FDis = NA, qual.FRic = NA, SE.qual.FRic = NA)
mean.cols <- c(2, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21)
sd.cols <- c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22)

(start <- Sys.time())
for(t in 1:nc) {
  cat("resampling bin =", mids[t], "Mya,", round(Sys.time() - start, 0), 
      "min since started.\n")

  # Initialize parallel cluster:
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  # Use load-balancing because of different run times for the MCMC optimizations
  opts <- list(preschedule = FALSE)
  clusterSetRNGStream(cl, 3142)  # Set L-Ecuyer RNG seed
  
  # Build numrep in parallel
  par.bin.metrics <- foreach(i = 1:(nreps * length(anc)), .options.snow = opts, 
                             .combine = rbind, .inorder = FALSE, 
                             .packages = "FD") %dopar% {

    # Create data.frame to store subsampled metrics
    sam.metrics <- data.frame(Age = as.numeric(mids[t]), S = std.g, H = NA, 
                              D = NA, M = NA, V = NA, R = NA, FRic = NA, 
                              FEve = NA, FDiv = NA, FDis = NA, qual.FRic = NA)

    FD <- rep(NA, 11)
    # Restrict to taxa (from relevant 'taxon.bin') occurring in bin 't' from
    # tree (specified in 'tree.seq')
    wh.gr <- unname(which(taxon.bins[[tree.seq[i]]][, t]))
    S <- length(wh.gr)

    # Only calculate next stats if more taxa than sampling quota
    if (S >= std.g) {
      
      # Extract relevant data for selected time-scaled tree
      anc.tree <- anc[[tree.seq[i]]]$matrix_1$matrix[wh.gr, ]
      dist.anc.tree <- dist.matrix[[tree.seq[i]]]$distance_matrix[wh.gr, wh.gr]
      pcoa.tree <- pcoa.results[[tree.seq[i]]]
      H <- nrow(unique(dist.anc.tree))

      # Only calculate next stats if complete distance matrix for tree sample
      if (!any(is.nan(dist.anc.tree)) & !length(dist.anc.tree) == 0) {

        # And only if more taxa or unique morphotypes/life habits than 'm'
        if (S > m | H > m) {
          
          # Create a new subsample each replicate
          max.seq <- seq.int(nrow(anc.tree))
          sampled <- sample2(max.seq, std.g, replace = FALSE)
          sub.tree <- anc.tree[sampled, ]
          sub.dist <- dist.anc.tree[sampled, sampled]
          sub.pcoa <- pcoa.tree
          sub.pcoa$vectors.cor <- sub.pcoa$vectors.cor[wh.gr, ]
          if (!is.null(sub.pcoa))
            sub.pcoa$vectors.cor <- sub.pcoa$vectors.cor[sampled,]
          
          # Calculate metrics for each subsample
          FD <- calc_metrics2(sample = sub.tree, dist.sam = sub.dist, 
                              pcoa = sub.pcoa, m = m, stand.pcoa = TRUE, 
                              calc.FRic.and.FDiv = TRUE)
          }
        }
    }
    sam.metrics[,2:12] <- FD
  }
  stopCluster(cl)
  
  # Summarize bin metrics
  metrics[t, mean.cols] <- apply(par.bin.metrics, 2, mean, na.rm = TRUE)
  metrics[t, sd.cols] <- apply(par.bin.metrics[, -1], 2, sd, na.rm = TRUE)
}
(Sys.time() - start) # 40-46 minutes using 8 core laptop, 3000 replicates

## Save / reload metrics:
# morph.stdG.metrics <- metrics; save(morph.stdG.metrics, file = "morph.stdG.metrics")
# mode.stdG.metrics <- metrics; save(mode.stdG.metrics, file = "mode.stdG.metrics")
# constant.stdG.metrics <- metrics; save(constant.stdG.metrics, file = "constant.stdG.metrics")
# raw.stdG.metrics <- metrics; save(raw.stdG.metrics, file = "raw.stdG.metrics")
# write.csv(metrics, file = "metrics_StdG50_morph.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_mode.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_constant.csv", row.names = FALSE)
# write.csv(metrics, file = "metrics_StdG50_LH_raw.csv", row.names = FALSE)
# metrics <- read.csv(file = "metrics_StdG50_morph.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_mode.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_constant.csv", header = TRUE)
# metrics <- read.csv(file = "metrics_StdG50_LH_raw.csv", header = TRUE)
# load("morph.stdG.metrics"); metrics <- morph.stdG.metrics
# load("mode.stdG.metrics"); metrics <- mode.stdG.metrics
# load("constant.stdG.metrics"); metrics <- constant.stdG.metrics
# load("raw.stdG.metrics"); metrics <- raw.stdG.metrics

beep(3)

head(metrics, 3)

# With m = 6 PCoA axes used for calculation of FRic and FDiv, what proportion of
# the variability is explained in different time bins? (I.e., how much reduction
# has occurred by the ordination?)
summary(metrics$qual.FRic)
# Morphology: 4.2% used, on average
# Ecology: mode: 2.8% used, constant: 1.9% used, raw: 2.1% used (but irrelevant)

# Plot trends
means <- c(3, 5, 7, 9, 11, 13, 15, 17, 19) # Columns with mean values
SEs <- c(4, 6, 8, 10, 12, 14, 16, 18, 20)  # Columns with SE values
par(mar = c(0, 4, 4, 2))
for (c in 1:9) {
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
means <- c(3, 5, 7, 9, 11, 13, 15, 17, 19)

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
# "H" (and warning) in comparisons with the morphological data set.

# morphological vs. LH-mode: not very correlated (R & FDiv highest)
round(diag(cor(diff.morph[, means], diff.mode[, means], 
               use = "pairwise.complete.obs")), 4)

# Are the D and FRic statistically correlated? (D no; FRic yes, but low R)
summary(lm(diff.morph$D ~ diff.mode$D))
summary(lm(diff.morph$FRic ~ diff.mode$FRic))

# morphological vs. LH-constant: not very correlated (R & FRic highest)
round(diag(cor(diff.morph[, means], diff.constant[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-constant: Moderately positively, except FEve; mean r = 0.71
round(diag(cor(diff.mode[, means], diff.constant[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-mode vs. LH-raw: Moderately: Moderately positively, except FEve; mean r = 0.43
round(diag(cor(diff.mode[, means], diff.raw[, means], 
               use = "pairwise.complete.obs")), 4)

# LH-raw vs. LH-constant: Moderately positively: mean r = 0.64
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
for (c in c(2:11)) {
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
std.g <- 50

# Plot trends
# pdf(file = "3EcoTrends_StdG.pdf"); 
par(mar = c(0, 4, 4, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17, 19)
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
  lines(mids, var_r, lwd = 5, lty = 1, col = cols[1])
  lines(mids, var_m, lwd = 5, lty = 2, col = cols[2])
  lines(mids, var_c, lwd = 5, lty = 3, col = cols[3])
  legend("topleft", legend = c("raw", "mode", "constant"), col = cols, 
         bty = "n", lty = c(1, 2, 3), lwd = 4, inset = 0.05)
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
load("morph.stdG.metrics"); metrics_morph <- morph.stdG.metrics
load("mode.stdG.metrics"); metrics_mode <- mode.stdG.metrics
load("constant.stdG.metrics"); metrics_constant <- constant.stdG.metrics
load("raw.stdG.metrics"); metrics_raw <- raw.stdG.metrics
std.g <- 50

# Because of rounding errors, several intervals in the morphological data set
# contain H values that are rounded to 50 life habits, but negligibly different.
# (E.g., 49.999 instead of 50.) The 100% maximum transformation below causes
# this to be treated as a true value less than 50. The following code overrides
# rounded-to-50 values with 50 so that the standardized plots ar sensible.
metrics_morph$H <- replace(metrics_morph$H, metrics_morph$H > 49.99, 50)

# Plot trends
# pdf(file = "Morph&3EcoTrends_StdG.pdf")
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17, 19)
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
         col = cols, bty = "n", lty = c(1, 2, 3, 4), lwd = 4, inset = 0.05)
}
par(op)
# dev.off()




# Same, but just mode and morph
# pdf(file = "2Morph&EcoTrends_StdG.pdf")
par(mar = c(0, 4, 2, 2))
mids <- metrics_mode$Age
var.cols <- c(3, 5, 7, 9, 11, 13, 15, 17, 19)
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

## ONLY PLOTS TREE #50 ***

# Load Wills GED-50 distance matrices
load("mode.distances.GED.5"); mode.distances.GED.5 <- mode.distances.GED.5[[50]]
load("constant.distances.GED.5"); constant.distances.GED.5 <- constant.distances.GED.5[[50]]
load("raw.distances.GED.5"); raw.distances.GED.5 <- raw.distances.GED.5[[50]]
load("morph.distances.GED.5"); morph.distances.GED.5 <- morph.distances.GED.5[[50]]

# Load PCoA output from ape::pcoa
load("mode.pcoa"); mode.pcoa <- mode.pcoa[[50]]
load("constant.pcoa"); constant.pcoa <- constant.pcoa[[50]]
load("raw.pcoa"); raw.pcoa <- raw.pcoa[[50]]
load("morph.pcoa"); morph.pcoa <- morph.pcoa[[50]]

# Set phylomorphospace/phyloecospace plotting colors
par(mar = c(5, 4, 2, 2))
branch.col <- "gray50"
tip.col <- "#0000FF7F"   # Set blue transparent so overlays as density
node.col <- "#A020F07F"  # Set purple transparent so overlays as density
root.col <- "#FDE725FF"  # viridisLite::viridis(3)[3]

# Set phylomorphospace/phyloecospace plotting variables
tree <- morph.pcoa$tree
tip.seq <- 1:Ntip(tree)
node.seq <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(tree$edge)), 
                                as.character(tree$edge[, 2])))
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
phytools::phylomorphospace(tree = tree, X = morph.pcoa$vectors.cor[tip.seq, 1:2], 
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
phytools::phylomorphospace(tree = tree, X = morph.pcoa$vectors.cor[tip.seq, 3:4], 
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
phytools::phylomorphospace(tree = tree, X = morph.pcoa$vectors.cor[tip.seq, 5:6], 
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
phylomorphospace3d(tree = morph.pcoa$tree, X = morph.pcoa$vectors.cor[tip.seq, 1:3], 
                 A = morph.pcoa$vectors.cor[node.seq, 1:3], 
                 control = c(con, list(ftype = "off")))
par(op)


## Same, with 'mode' ecological data set
# Set phylomorphospace/phylomorphospace plotting variables
tree <- mode.pcoa$tree
tip.seq <- 1:Ntip(tree)
node.seq <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(tree$edge)), 
                                as.character(tree$edge[, 2])))
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
phytools::phylomorphospace(tree = tree, X = mode.pcoa$vectors.cor[tip.seq, 1:2], 
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
phytools::phylomorphospace(tree = tree, X = mode.pcoa$vectors.cor[tip.seq, 3:4], 
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
phytools::phylomorphospace(tree = tree, X = mode.pcoa$vectors.cor[tip.seq, 5:6], 
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
phylomorphospace3d(tree = mode.pcoa$tree, X = mode.pcoa$vectors.cor[tip.seq, 1:3], 
                   A = mode.pcoa$vectors.cor[node.seq, 1:3], 
                   control = c(con, list(ftype = "off")))
par(op)





## Same, with 'constant' ecological data set
# Set phylomorphospace/phylomorphospace plotting variables
# Plotting uncorrected eigenvectors instead of corrected ones
tree <- constant.pcoa$tree
tip.seq <- 1:Ntip(tree)
node.seq <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(tree$edge)), 
                                as.character(tree$edge[, 2])))
x.lab1 <- paste0("PCO 1 (", round(100 * constant.pcoa$values$Rel_corr_eig[1], 1), 
                 "% of total variance)")
y.lab2 <- paste0("PCO 2 (", round(100 * constant.pcoa$values$Rel_corr_eig[2], 1), 
                 "% of total variance)")
x.lab3 <- paste0("PCO 3 (", round(100 * constant.pcoa$values$Rel_corr_eig[3], 1), 
                 "% of total variance)")
y.lab4 <- paste0("PCO 4 (", round(100 * constant.pcoa$values$Rel_corr_eig[4], 1), 
                 "% of total variance)")
x.lab5 <- paste0("PCO 5 (", round(100 * constant.pcoa$values$Rel_corr_eig[5], 1), 
                 "% of total variance)")
y.lab6 <- paste0("PCO 6 (", round(100 * constant.pcoa$values$Rel_corr_eig[6], 1), 
                 "% of total variance)")

# pdf(file = "ecospace_constant.pdf")
par(mfrow = c(2, 2))
# PCO 1 vs. PCO 2
phytools::phylomorphospace(tree = tree, X = constant.pcoa$vectors.cor[tip.seq, 1:2],
                           A = constant.pcoa$vectors.cor[node.seq, 1:2],
                           control = con, label = "off", xlab = x.lab1,
                           ylab = y.lab2, pch = NA)
points(x = constant.pcoa$vectors.cor[node.seq, 1],
       y = constant.pcoa$vectors.cor[node.seq, 2], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[tip.seq, 1], 
       y = constant.pcoa$vectors.cor[tip.seq, 2], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 1], 
       y = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 2], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 3 vs. PCO 4
phytools::phylomorphospace(tree = tree, X = constant.pcoa$vectors.cor[tip.seq, 3:4],
                           A = constant.pcoa$vectors.cor[node.seq, 3:4],
                           control = con, label = "off", xlab = x.lab3,
                           ylab = y.lab4, pch = NA)
points(x = constant.pcoa$vectors.cor[node.seq, 3],
       y = constant.pcoa$vectors.cor[node.seq, 4], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[tip.seq, 3], 
       y = constant.pcoa$vectors.cor[tip.seq, 4], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 3], 
       y = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 4], col = root.col, pch = 16, 
       cex = 1.25)

# PCO 5 vs. PCO 6
phytools::phylomorphospace(tree = tree, X = constant.pcoa$vectors.cor[tip.seq, 5:6],
                           A = constant.pcoa$vectors.cor[node.seq, 5:6],
                           control = con, label = "off", xlab = x.lab5,
                           ylab = y.lab6, pch = NA)
points(x = constant.pcoa$vectors.cor[node.seq, 5],
       y = constant.pcoa$vectors.cor[node.seq, 6], col = node.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[tip.seq, 5], 
       y = constant.pcoa$vectors.cor[tip.seq, 6], col = tip.col, pch = 16, cex = 1.25)
points(x = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 5], 
       y = constant.pcoa$vectors.cor[(max(tip.seq) + 1), 6], col = root.col, pch = 16, 
       cex = 1.25)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "phyloecospace (constant)", legend = c("tips", "nodes", "root"), 
       cex = 2, pch = 16, col = c(tip.col, node.col, root.col), bty = "n", 
       pt.cex = 4)
# dev.off()
par(op)


# First 3 in 3-D (wait a few seconds while renders)
phylomorphospace3d(tree = constant.pcoa$tree, X = constant.pcoa$vectors.cor[tip.seq, 1:3], 
                A = constant.pcoa$vectors.cor[node.seq, 1:3], 
                control = c(con, list(ftype = "off")))
par(op)





## Same, with 'raw' ecological data set

## For 'raw' treatment (ONLY!), need to remove the taxa trimmed out for the PCOA
## because of high density of missing states.
# Use 'raw.trimmed' output to override the raw objects. (Not needed for
# pcoa.results because used earlier trimmed version when building the PCoA
# object.) Because the tip names are based on the raw.anc[[x]]$topper$tree, need
# to match names for rows to correctly identify the trimmed-out tips and nodes
load("raw.trimmed"); raw.trimmed <- raw.trimmed[[50]]
load("raw.anc"); raw.anc <- raw.anc[[50]]
load("raw.distances.GED.5"); raw.distances.GED.5 <- raw.distances.GED.5[[50]]
load("raw.pcoa"); raw.pcoa <- raw.pcoa[[50]]
taxa.to.cut <- raw.trimmed$removed_taxa_by_name
wh.to.cut <- match(taxa.to.cut, rownames(raw.anc$matrix_1$matrix))
raw.anc$matrix_1$matrix <- raw.anc$matrix_1$matrix[-wh.to.cut, ]
raw.distances.GED.5$distance_matrix <-
  raw.distances.GED.5$distance_matrix[-wh.to.cut, -wh.to.cut]
# Confirm all same number of rows
nrow(raw.anc$matrix_1$matrix)
nrow(raw.distances.GED.5$distance_matrix)
nrow(raw.pcoa$vectors.cor)


# Set phylomorphospace/phylomorphospace plotting variables
tree <- raw.trimmed$tree  # Note using the trimmed tree here (unlike above)
tip.seq <- 1:Ntip(tree)
node.seq <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
con <- list(col.edge = setNames(rep(branch.col, nrow(tree$edge)), 
                                as.character(tree$edge[, 2])))
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
phytools::phylomorphospace(tree = tree, X = raw.pcoa$vectors.cor[tip.seq, 1:2],
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
phytools::phylomorphospace(tree = tree, X = raw.pcoa$vectors.cor[tip.seq, 3:4], 
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
phytools::phylomorphospace(tree = tree, X = raw.pcoa$vectors.cor[tip.seq, 5:6], 
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
phylomorphospace3d(tree = raw.distances.GED.5$tree, 
                   X = raw.pcoa$vectors.cor[tip.seq, 1:3], 
                   A = raw.pcoa$vectors.cor[node.seq, 1:3],
                   control = c(con, list(ftype = "off")))
par(op)




## KMEANS CLUSTERING TO INTERPRET SPACES #######################################

## ONLY USING TREE #50 ***

# Load PCoA output from ape::pcoa
load("mode.pcoa"); mode.pcoa <- mode.pcoa[[50]]
load("constant.pcoa"); constant.pcoa <- constant.pcoa[[50]]
load("raw.pcoa"); raw.pcoa <- raw.pcoa[[50]]
load("morph.pcoa"); morph.pcoa <- morph.pcoa[[50]]

load("taxon.list"); taxon.list <- taxon.list[[50]]
if(nrow(morph.pcoa$vectors.cor) != nrow(taxon.list))
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
# If phylogenetic structure, lots of zeros:
(cl.table <- table(taxon.list[, "class"], km$cluster))
sort(cl.table[ ,1])
sort(cl.table[ ,2])
sort(cl.table[ ,3])
sort(cl.table[ ,4])
legend.groups <- c("edrio, aster, ech, oph & cyclo", 
                   "stylo, homo, solut, cteno, helico",
                   "most crinoids", 
                   "more crinoids, eocr, rhomb, diplo, paracr")
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
(cl.table <- table(taxon.list[, "class"], km$cluster))
sort(cl.table[ ,1])
sort(cl.table[ ,2])
sort(cl.table[ ,3])
sort(cl.table[ ,4])
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
no.axes <- 6
ks <- c(1:10, seq(20, 50, by = 10), 75, 100)
sum.sq <- array(dim = length(ks))
for(k in 1:length(ks)){
  sum.sq[k] <- kmeans(constant.pcoa$vectors.cor[, 1:no.axes], centers = ks[k],
                      nstart = 25, iter.max = 1000)$tot.withinss
}
plot(ks, sum.sq, type = "b", pch = 16, xlab = "No. k-means clusters",
     ylab = "Total within-cluster sum of squares",
     main = "Ecological (constant) PCoA k-means clusters")
par(op)
# dev.off()

# For constant, k = 3 (or 4) provides clean breaks; 9 divides most classes
# pdf(file = "kmeans_eco_constant.pdf")
k <- 4
set.seed(3) # To allow replication of order
km <- kmeans(constant.pcoa$vectors.cor[, 1:6], centers = k, nstart = 25, iter.max = 100)
# Modify so matches 'mode' (switch 1 and 4)
k.1 <- which(km$cluster == 1)
k.4 <- which(km$cluster == 4)
km$cluster[k.1] <- 4
km$cluster[k.4] <- 1
par(mfrow = c(2, 2), mar = c(4, 4, 1, 0.25))
cols <- plasma(k)[km$cluster]
pchs <-as.character(km$cluster)
plot(constant.pcoa$vectors.cor[, 1:2], col = cols, pch = pchs, cex = 0.75)
plot(constant.pcoa$vectors.cor[, 3:4], col = cols, pch = pchs, cex = 0.75)
plot(constant.pcoa$vectors.cor[, 5:6], col = cols, pch = pchs, cex = 0.75)
(cl.table <- table(taxon.list[, "class"], km$cluster))
sort(cl.table[ ,1])
sort(cl.table[ ,2])
sort(cl.table[ ,3])
sort(cl.table[ ,4])
legend.groups <- c("crinoids (& some eocr, edrio & rhomb)",
                   "edr, eocr, dipl, paracr, rh, hel & more cri",
                   "stylo, homo, solut, ctenoc, cyclo, some rhomb",
                   "aster, ech, oph, stenur, somas & holo")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes = FALSE, xlab="", ylab = "")
legend("left", title = "Ecological PCoA", legend = legend.groups, cex = 1,
       pch = as.character(1:k), col = plasma(k)[1:k], bty = "n", pt.cex = 1.5)
# dev.off()
par(op)






## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF CLASSES ################
# Question: Within taxonomic classes, is the overall occupation of space greater
# in morphospace or ecospace?

# No need to sample standardize here, as the sample size is identical for each
# ecological and morphological subset. Simply need to observe the differences in
# convex hull hypervolume (FRic), using 'm' PCoA axes like above. And this
# analysis is not observing the trends through time. Standardizing each space to
# the overall FRic of all echinoderms because the length of eigenvalues may be
# different in each space.

# Process PCoA output
load("morph.pcoa")
load("mode.pcoa")
eco.pcoa <- mode.pcoa
# load("constant.pcoa"); eco.pcoa <- constant.pcoa
# load("raw.pcoa"); eco.pcoa <- raw.pcoa
load("taxon.list")

# Select 'm' (number of PCoA axes to use when calculating convex hull volume).
# Use m = 2:7
m <- 2

class.space <- vector("list", length(morph.pcoa))

# Identify classes (which are constant across taxon.lists:
classes <- sort(unique(taxon.list[[1]][, "class"]))
# Remove genera with UNCERTAIN class affiliation and 'diploporitan' (which will
# be processed with Diploporita sensu strictu)
classes <- classes[-which(classes == "UNCERTAIN" | classes == "'diploporitan'")]

(clu <- makeCluster(detectCores()))
registerDoParallel(clu)
(start <- Sys.time())
class.space <- foreach(t = 1:length(morph.pcoa), .inorder = TRUE, 
                        .packages = "geometry") %dopar% {
  
  t.class.space <- matrix(NA, nrow = length(classes), ncol = 4)
  colnames(t.class.space) <- c("morph", "eco", "diff", "diff.prop")
  
  
  # Only proceed if tree objects exist (as in trees 31-35)
  if (!is.null(morph.pcoa[[t]])) {
  
    # Standardize the PCoA eigenvectors
    morph.pcoa[[t]]$vectors.cor <- 
      stand.pcoa(vectors = morph.pcoa[[t]]$vectors.cor,
                 eigenvalues = morph.pcoa[[t]]$values$Corr_eig)
    eco.pcoa[[t]]$vectors.cor <- 
      stand.pcoa(vectors = eco.pcoa[[t]]$vectors.cor,
                 eigenvalues = eco.pcoa[[t]]$values$Corr_eig)
    
    for(cl in 1:nrow(t.class.space)) {
  
      # Combine Diploporita with 'diploporitans'
      if (classes[cl] == "Diploporita") {
        wh.cl <- which(taxon.list[[t]][, "class"] == classes[cl] |
                         taxon.list[[t]][, "class"] == "'diploporitan'")
      } else
        wh.cl <- which(taxon.list[[t]][, "class"] == classes[cl])
      # Skip to next class if < 2 genera in class
      if (length(wh.cl) < 2L)
        next
      
      eco.coords <- eco.pcoa[[t]]$vectors.cor[wh.cl, 1:m]
      morph.coords <- morph.pcoa[[t]]$vectors.cor[wh.cl, 1:m]
  
      # Can not calculate if there are less than 'm' distinct points
      H.eco <- nrow(unique(round(eco.coords, 5)))
      H.morph <- nrow(unique(round(morph.coords, 5)))
      if (H.eco <= m | H.morph <= m) {
        warning("Convex hull can not be calculated for 'm' = ", m, 
                " for class ", classes[cl], " in tree ", t, "\n")
        next
        }
        
      # Standardize by entire morphospace / ecospace
      morph.total.attempt <-
        try(geometry::convhulln(morph.pcoa[[t]]$vectors.cor[, 1:m], "FA")$vol)
      eco.total.attempt <- 
        try(geometry::convhulln(eco.pcoa[[t]]$vectors.cor[, 1:m], "FA")$vol)
      morph.cl.attempt <- try(geometry::convhulln(morph.coords, "FA")$vol)
      eco.cl.attempt <- try(geometry::convhulln(eco.coords, "FA")$vol)
      if (!inherits(morph.total.attempt, "try-error") &
          !inherits(eco.total.attempt, "try-error") &
          !inherits(eco.cl.attempt, "try-error") &
          !inherits(eco.cl.attempt, "try-error")) {
        t.class.space[cl, "morph"] <- morph.cl.attempt / morph.total.attempt
        t.class.space[cl, "eco"] <- eco.cl.attempt / eco.total.attempt
      }
    }
    
    # Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
    t.class.space[, "diff"] <- t.class.space[, "morph"] - t.class.space[, "eco"]
    t.class.space[, "diff.prop"] <- 
      t.class.space[, "morph"] / t.class.space[, "eco"]
  }
  return(t.class.space)
}
stopCluster(clu)
(Sys.time() - start) # 22-68 s, depending on 'm'. 5.6 minutes when m = 7
beep()

# Medians across trees. (Most are close to normally distributed, but some
# classes are definitely not, especially for the proportional difference.)
class.results <-
  as.data.frame(apply(simplify2array(class.space), 1:2, median, na.rm = TRUE))
class.results <- cbind(classes, class.results)
nrow(na.omit(class.results))
na.omit(class.results)

summary(class.results$diff)
summary(class.results$diff.prop)

# Which classes have GREATER morphospace occupation compared to ecological?

#   Crinoids (only when m = 2-4), ctenocystoids (only when m = 3 & 4),
#   diploporitans (m < 7), edrioasteroids (m = 2-3 & 6-7), eocrinoids,
#   homosteleans (m = 2 & 3), ophiuroids (m = 2-4), parablastoids (m = 2),
#   paracrinoids (m = 2-5 & 7), rhombiferans (m = 5-7), and solutes (m = 6-7).

# Which classes have LESSER morphospace occupation compared to ecological?

#   Asteroids, crinoids (only when m = 5-7), ctenocyctoids (m = 2),
#   cyclocystoids (m = 2-5), edrioasteroids (m = 4), holothuroids (m = 2),
#   ophiocistioids (m = 2), ophiuroids (m = 2), paracrinoids (m = 5),
#   rhombiferans (m = 2-4), solutes (m = 2-4), stenuroids (m = 2), and
#   stylophorans (m = 2-4 & 6-7)

# Statistical tests
par(mfrow = c(2, 1))
hist(class.results$diff, 15, main = "differences in class convex hull volume", 
     xlab = "morph volume - eco volume, by class")
abline(v = 0, lwd = 2, lty = 2)
plot(density(na.omit(class.results$diff)), main = "differences in class convex hull volume", 
     xlab = "morph volume - eco volume, by class")
abline(v = 0, lwd = 2, lty = 2)

# Distribution is normally distributed
shapiro.test(class.results$diff) # p > 0.05 so normally distributed (except m > 4)
 
t.test(class.results$diff, mu = 0)
# p = 0.5993 - 0.7871 (depending on 'm'), so not different

# Confirm with paired t-test
t.test(class.results$morph, class.results$eco, paired = TRUE)
# p = 0.6072 - 0.7884 (depending on 'm'), so not diff





## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF SUBPHYLA ###############
# Question: Within taxonomic subphyla, is the overall occupation of space
# greater in morphospace or ecospace?

# Select 'm' (number of PCoA axes to use when calculating convex hull volume).
# Use m = 2:7
m <- 2

subphylum.space <- vector("list", length(morph.pcoa))

# Identify subphyla (which are constant across taxon.lists:
subphyla <- sort(unique(taxon.list[[1]][, "subphylum"]))
# Remove genera with UNCERTAIN subphylum affiliation
subphyla <- subphyla[-which(subphyla == "UNCERTAIN")]

(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
subphylum.space <- foreach(t = 1:length(morph.pcoa), .inorder = TRUE, 
                       .packages = "geometry") %dopar% {
                         
   t.subphylum.space <- matrix(NA, nrow = length(subphyla), ncol = 4)
   colnames(t.subphylum.space) <- c("morph", "eco", "diff", "diff.prop")

   # Only proceed if tree objects exist (as in trees 31-35)
   if (!is.null(morph.pcoa[[t]])) {
     
     # Standardize the PCoA eigenvectors
     morph.pcoa[[t]]$vectors.cor <- 
       stand.pcoa(vectors = morph.pcoa[[t]]$vectors.cor,
                  eigenvalues = morph.pcoa[[t]]$values$Corr_eig)
     eco.pcoa[[t]]$vectors.cor <- 
       stand.pcoa(vectors = eco.pcoa[[t]]$vectors.cor,
                  eigenvalues = eco.pcoa[[t]]$values$Corr_eig)
     
     for(sp in 1:nrow(t.subphylum.space)) {
       wh.sp <- which(taxon.list[[t]][, "subphylum"] == subphyla[sp])
       
       # Skip to next subphylum if < 2 genera in subphylum
       if (length(wh.sp) < 2L)
         next
       
       eco.coords <- eco.pcoa[[t]]$vectors.cor[wh.sp, 1:m]
       morph.coords <- morph.pcoa[[t]]$vectors.cor[wh.sp, 1:m]
       
       # Can not calculate if there are less than 'm' distinct points
       H.eco <- nrow(unique(round(eco.coords, 5)))
       H.morph <- nrow(unique(round(morph.coords, 5)))
       if (H.eco <= m | H.morph <= m) {
         warning("Convex hull can not be calculated for 'm' = ", m, 
                 " for subphylum ", subphyla[sp], " in tree ", t, "\n")
         next
       }
       
       # Standardize by entire morphospace / ecospace
       morph.total.attempt <-
         try(geometry::convhulln(morph.pcoa[[t]]$vectors.cor[, 1:m], "FA")$vol)
       eco.total.attempt <- 
         try(geometry::convhulln(eco.pcoa[[t]]$vectors.cor[, 1:m], "FA")$vol)
       morph.cl.attempt <- try(geometry::convhulln(morph.coords, "FA")$vol)
       eco.cl.attempt <- try(geometry::convhulln(eco.coords, "FA")$vol)
       if (!inherits(morph.total.attempt, "try-error") &
           !inherits(eco.total.attempt, "try-error") &
           !inherits(eco.cl.attempt, "try-error") &
           !inherits(eco.cl.attempt, "try-error")) {
         t.subphylum.space[sp, "morph"] <- morph.cl.attempt / morph.total.attempt
         t.subphylum.space[sp, "eco"] <- eco.cl.attempt / eco.total.attempt
       }
     }
     
     # Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
     t.subphylum.space[, "diff"] <- t.subphylum.space[, "morph"] - t.subphylum.space[, "eco"]
     t.subphylum.space[, "diff.prop"] <- 
       t.subphylum.space[, "morph"] / t.subphylum.space[, "eco"]
   }
   return(t.subphylum.space)
}
stopCluster(cl)
(Sys.time() - start) # 18-26 s, depending on 'm'. 52 s when m = 6, 4.1 minutes when m = 7
beep()

# Medians across trees. (Most are close to normally distributed, but some
# subphyla are definitely not, especially for the proportional difference.)
subphylum.results <-
  as.data.frame(apply(simplify2array(subphylum.space), 1:2, median, na.rm = TRUE))
subphylum.results <- cbind(subphyla, subphylum.results)
nrow(na.omit(subphylum.results))
na.omit(subphylum.results)

summary(subphylum.results$diff)
summary(subphylum.results$diff.prop)

# Which subphyla have GREATER morphospace occupation compared to ecological?
#   Blastozoa, Crinozoa, and Echinozoa

# Which subphyla have LESSER morphospace occupation compared to ecological?
#   Asterozoa

# Statistical tests
par(mfrow = c(2, 1))
hist(subphylum.results$diff, 15, main = "differences in subphylum convex hull volume", 
     xlab = "morph volume - eco volume, by subphylum")
abline(v = 0, lwd = 2, lty = 2)
plot(density(na.omit(subphylum.results$diff)), main = "differences in subphylum convex hull volume", 
     xlab = "morph volume - eco volume, by subphylum")
abline(v = 0, lwd = 2, lty = 2)

# Distribution is normally distributed
shapiro.test(subphylum.results$diff) # p > 0.05 so normally distributed (except m > 4)

t.test(subphylum.results$diff, mu = 0)
# p = 0.1914 - 0.5336 (depending on 'm'), so not different

# Confirm with paired t-test
t.test(subphylum.results$morph, subphylum.results$eco, paired = TRUE)
# p = 0.1914 - 0.5310 (depending on 'm'), so not diff






## COMPARE OVERALL MORPHOSPACE / ECOSPACE OCCUPATION OF PHYLUM ###############
# Question: Within all echinoderms, is the OVERALL occupation of space greater
# in morphospace or ecospace?

# Select 'm' (number of PCoA axes to use when calculating convex hull volume)
m <- 7 # Use values of 2-7

phylum.space <- vector("list", length(morph.pcoa))

(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
phylum.space <- foreach(t = 1:length(morph.pcoa), .inorder = TRUE, 
                        .packages = "geometry") %dopar% { 

  # Only proceed if tree objects exist (as in trees 31-35)
  if (!is.null(morph.pcoa[[t]])) {
    t.phylum.space <- rep(NA, 4)
    
    # Standardize the PCoA eigenvectors
    morph.pcoa[[t]]$vectors.cor <- 
      stand.pcoa(vectors = morph.pcoa[[t]]$vectors.cor,
                 eigenvalues = morph.pcoa[[t]]$values$Corr_eig)
    eco.pcoa[[t]]$vectors.cor <- 
      stand.pcoa(vectors = eco.pcoa[[t]]$vectors.cor,
                 eigenvalues = eco.pcoa[[t]]$values$Corr_eig)
  
    eco.coords <- eco.pcoa[[t]]$vectors.cor[, 1:m]
    morph.coords <- morph.pcoa[[t]]$vectors.cor[, 1:m]
    
    # Can not calculate if there are less than 'm' distinct points
    H.eco <- nrow(unique(round(eco.coords, 5)))
    H.morph <- nrow(unique(round(morph.coords, 5)))
    if (H.eco <= m | H.morph <= m) {
      warning("Convex hull can not be calculated for 'm' = ", m, " for tree ", t, "\n")
      stop
    }
    morph.attempt <- try(geometry::convhulln(morph.coords, "FA")$vol)
    eco.attempt <- try(geometry::convhulln(eco.coords, "FA")$vol)
    if (!inherits(morph.attempt, "try-error") & !inherits(eco.attempt, "try-error")) {
      t.phylum.space[1] <- morph.attempt
      t.phylum.space[2] <- eco.attempt
      
      # Calculate difference. If +, morph > eco. If 0, morph = eco, If -, morph < eco
      t.phylum.space[3] <- t.phylum.space[1] - t.phylum.space[2]
      t.phylum.space[4] <- t.phylum.space[1] / t.phylum.space[2]
    }
  }
  return(t.phylum.space)
}
stopCluster(cl)
(Sys.time() - start) # 20-68 s, depending on 'm'
beep()

# Averages across trees:
phylum.results <- matrix(apply(simplify2array(phylum.space), 1, mean, na.rm = TRUE), ncol = 4)
colnames(phylum.results) <- c("morph", "eco", "diff", "prop.diff")
phylum.results

# Overall morphospace is 2.95-times greater than ecospace
# m = 2: morph is 2.28x greater
# m = 3:  3.87x
# m = 4:  8.40x
# m = 5: 13.16x
# m = 6: 14.57x
# m = 7: 20.09x





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
Tree <- morph.pcoa$tree  # Same tree topology for both morphology and ecology
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
Tree <- morph.pcoa$tree  # Same tree topology for both morphology and ecology
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
load("morph.distances.GED.5")
# Load for the $trees
load("morph.anc")
load("mode.anc")

# Convert distance matrixes (restricting to tips) to cluster analysis
nt <- length(morph.distances.GED.5)
tr.morph <- vector("list", nt)
for(t in 1:nt) {
  tree <- morph.anc[[t]]$topper$tree
  morph.hc <-
    hclust(as.dist(morph.distances.GED.5[[t]]$distance_matrix[1:Ntip(tree), 1:Ntip(tree)]))
  tr.morph[[t]] <- as.phylo(morph.hc)
  tr.morph[[t]]$tip.label <-  tree$tip.label
}
plot(tr.morph[[50]])

# Cophyloplot between phylogeny and morpho-cluster
# 'association' not needed because strictly matches tip labels by default
cophylo.morph <- vector("list", length(tr.morph))
(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
cophylo.morph <- foreach(t = 1:nt, .inorder = TRUE, .packages = "phytools") %dopar% {
  tree <- morph.anc[[t]]$topper$tree
  t.cophylo.morph <- phytools::cophylo(tree, tr.morph[[t]], rotate = TRUE,
                                       rotate.multi = TRUE, print = FALSE)
}
stopCluster(cl)
(Sys.time() - start) # 13.1 minutes on 8-core laptop
beep(3)

# Confirm taxa are paired correctly
summary(cophylo.morph[[50]])

# save(cophylo.morph, file = "cophylo.morph")
# load("cophylo.morph")

# Plot sample (using tree 50)
plot(cophylo.morph[[50]])
# The phylogeny is the left figure

# plot.cophylo() allows following customizable plotting arguments:
# For the linkages: link.col, link.lwd, and link.lty; part controls width of 
#                   the linkages (but default 0.4 is best)
# For the cladogram/dendrogram: edge.col
# For the middle length: edge.col
# For font size: fsize
# For plotting points at tips: pts = T/F
# pdf(file = "MorphPhylogram.pdf")
edge.colors <- list(left = rep("darkgray", nrow(cophylo.morph[[50]]$trees[[1]]$edge)), 
                    right = rep("black", nrow(cophylo.morph[[50]]$trees[[2]]$edge)))
plot(cophylo.morph[[50]], edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
     fsize = 0.3, pts = FALSE)
# dev.off()

# How much sum-of-square rank differences between the two "phylogenies"? (same
# metric used in 'cophylo' to optimize the tip alignments)
(diffs.morph <- attr(cophylo.morph[[50]]$trees[[1]], "minRotate"))
sqrt(diffs.morph / Ntip(morph.anc[[50]]$topper$tree)) # Average rank difference

# What is the distribution of rank differences between the two "phylogenies"?
phylo.order <- cophylo.morph[[50]]$trees[[1]]$tip.label
morph.order <- cophylo.morph[[50]]$trees[[2]]$tip.label
match.offset <- match(phylo.order, morph.order)

# Confirm matched correctly:
identical(phylo.order, morph.order[match.offset])

# Confirm sum of squares matches that in the cophylogram:
morph.offset <- match.offset - seq(phylo.order)
identical(diffs.morph, sum(morph.offset ^ 2))

hist(morph.offset, main = "rank mismatch between phylogeny and morphogram", 
     cex.main = .9)

# pdf(file = "EchinoMorphPhylogeny.pdf", height = 200)
# plot(cophylo.morph[[50]], edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
#      fsize = 0.65, pts = FALSE)
# dev.off()

# Calculate offsets across trees
morph.offset <- vector("list", length(cophylo.morph))
for(t in 1:length(cophylo.morph)){
  phylo.order <- cophylo.morph[[t]]$trees[[1]]$tip.label
  morph.order <- cophylo.morph[[t]]$trees[[2]]$tip.label
  match.offset <- match(phylo.order, morph.order)
  if (!identical(phylo.order, morph.order[match.offset]))
    stop(cat("order does not match for tree", t, "\n"))
  morph.offset[[t]] <- match.offset - seq(phylo.order)
  diffs.morph <- attr(cophylo.morph[[t]]$trees[[1]], "minRotate")
  if (!identical(diffs.morph, sum(morph.offset[[t]] ^ 2)))
    stop(cat("sum of square offsets did not calculate correctly for tree", t, "\n"))
}


## Do the same for the ecology tree (using mode data set)
# Convert distance matrixes (restricting to tips) to cluster analysis
nt <- length(mode.distances.GED.5)
tr.eco <- vector("list", nt)
for(t in 1:nt) {
  tree <- mode.anc[[t]]$topper$tree
  eco.hc <-
    hclust(as.dist(mode.distances.GED.5[[t]]$distance_matrix[1:Ntip(tree), 1:Ntip(tree)]))
  tr.eco[[t]] <- as.phylo(eco.hc)
  tr.eco[[t]]$tip.label <-  tree$tip.label
}
plot(tr.eco[[50]])

# Cophyloplot between phylogeny and eco-cluster
# 'association' not needed because strictly matches tip labels by default
cophylo.eco <- vector("list", length(tr.eco))
(cl <- makeCluster(detectCores()))
registerDoParallel(cl)
(start <- Sys.time())
cophylo.eco <- foreach(t = 1:nt, .inorder = TRUE, .packages = "phytools") %dopar% {
  tree <- mode.anc[[t]]$topper$tree
  t.cophylo.eco <- phytools::cophylo(tree, tr.eco[[t]], rotate = TRUE,
                                     rotate.multi = TRUE, print = FALSE)
}
stopCluster(cl)
(Sys.time() - start) # 10.6 minutes on 8-core laptop
beep(3)

# Confirm taxa are paired correctly
summary(cophylo.eco[[50]])

# save(cophylo.eco, file = "cophylo.eco")
# load("cophylo.eco")

plot(cophylo.eco[[50]])
# The phylogeny is the left figure
par(op)

# pdf(file = "EcoPhylogram.pdf")
edge.colors <- list(left = rep("darkgray", nrow(cophylo.eco[[50]]$trees[[1]]$edge)), 
                    right = rep("blue", nrow(cophylo.eco[[50]]$trees[[2]]$edge)))
plot(cophylo.eco[[50]], edge.col = edge.colors, link.lwd = 1, link.col = "darkgray", 
     fsize = 0.3, pts = FALSE)
# dev.off()
par(op)

# How much sum-of-square rank differences between the two "phylogenies"? (same
# metric used in 'cophylo' to optimize the tip alignments)
(diffs.eco <- attr(cophylo.eco[[50]]$trees[[1]], "minRotate"))
sqrt(diffs.eco / Ntip(mode.anc[[50]]$topper$tree)) # Average rank difference

# What is the distribution of rank differences between the two "phylogenies"?
phylo.order <- cophylo.eco[[50]]$trees[[1]]$tip.label
eco.order <- cophylo.eco[[50]]$trees[[2]]$tip.label
match.offset <- match(phylo.order, eco.order)

# Confirm matched correctly:
identical(phylo.order, eco.order[match.offset])

# Confirm sum of squares matches that in the cophylogram:
eco.offset <- match.offset - seq(phylo.order)
identical(diffs.eco, sum(eco.offset ^ 2))

hist(eco.offset, main = "rank mismatch between phylogeny and ecogram", 
     cex.main = .9)

# pdf(file = "EchinoEcoPhylogeny.pdf", height = 200)
# plot(cophylo.eco[[50]], edge.col = edge.colors, link.lwd = 1, link.col = "darkgray",
#      fsize = 0.65, pts = FALSE)
# dev.off()

# Calculate offsets across trees
eco.offset <- vector("list", length(cophylo.eco))
for(t in 1:length(cophylo.eco)){
  phylo.order <- cophylo.eco[[t]]$trees[[1]]$tip.label
  eco.order <- cophylo.eco[[t]]$trees[[2]]$tip.label
  match.offset <- match(phylo.order, eco.order)
  if (!identical(phylo.order, eco.order[match.offset]))
    stop(cat("order does not match for tree", t, "\n"))
  eco.offset[[t]] <- match.offset - seq(phylo.order)
  diffs.eco <- attr(cophylo.eco[[t]]$trees[[1]], "minRotate")
  if (!identical(diffs.eco, sum(eco.offset[[t]] ^ 2)))
    stop(cat("sum of square offsets did not calculate correctly for tree", t, "\n"))
}



## Statistical comparisons: Are the offsets different in the two dendrograms?

summary(unlist(morph.offset))
sd(unlist(morph.offset))

summary(unlist(eco.offset))
sd(unlist(eco.offset))

# Because both distributions sum to zero, a t-test is not appropriate

# Non-parametric Mann-Whitney U-test to test whether different
wilcox.test(unlist(eco.offset), unlist(morph.offset))
# Not different in raw offsets: p = 0.913
wilcox.test(unlist(eco.offset) ^ 2, unlist(morph.offset) ^ 2)
# But sig. diff in sum of square offsets: p = 2.479e-12

# CONCLUSION: They overlap significantly, but there are differences in the
# tails.

# Seems the major difference is the much longer tails of the ecological
# distribution. Are they different distributional shapes?
ks.test(unlist(eco.offset), unlist(morph.offset))
ks.test(unlist(eco.offset) ^ 2, unlist(morph.offset) ^ 2)
# CONCLUSION: Distributions are different (p < 0.001)

# Let's focus instead on a test of different variances using the two-sided
# F-test:
var.test(unlist(eco.offset), unlist(morph.offset))
# F = 1.0925, p = 2.163e-09

# CONCLUSION: The ecological data set has 9.25% substantially greater variance
# than the morphological data set, implying a greater degree of phylogenetic
# structure in the morphological data set (and less in the ecological data set).

boxplot(unlist(morph.offset), unlist(eco.offset), col = "gray", yaxs = "i",
        names = c("morphology", "ecology"), 
        main = "paired offsets between phylogeny and dendrogram")

# pdf(file = "cophylogram_offsets.pdf")
breaks <- seq(-375, 375, 25)
hist(unlist(morph.offset), main = "Offsets between phylogeny and dendrogram", 
     xlab = "paired offsets", ylab = "density", breaks = breaks, 
     col = "transparent", border = "transparent", prob = TRUE)
hist(unlist(eco.offset), add = T, border = "white", col = "darkgray", 
     breaks = breaks, prob = TRUE)
hist(unlist(morph.offset), add = T, border = "black", col = "transparent", 
     breaks = breaks, prob = TRUE)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
# dev.off()
