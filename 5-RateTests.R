## CHARACTER RATE ANALYSES #####################################################

# Note that substantial coding changes accompanied the update to Claddis v.
# 0.6.0 in August 2020. The code here uses the functions as of 7/2/2020, v.
# 0.4.1, but with test_rates() modified from v. 0.6.0. Future users will need to
# either download the archived version of the package from GitHub or alter the
# code accordingly (primarily changes to function and argument names).

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' package.

## PREPARATIONS ################################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")
setwd("~/Manuscripts/CamOrdEchino/Data files/NA Reformatted")

# Load packages
library(beepr)      # v. 1.3
library(ape)        # v. 5.4
library(geiger)     # v. 2.0.7
library(doParallel) # v. 1.0.15
library(Claddis)    # v. 0.4.1 - Check SI for Lloyd 2018 for walk-through on code functions (in R experiments folder)
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub.\n")



## IMPORT FILES ################################################################

# Import ancestral states from 2-InferancestralStates.R
load("mode.anc")
load("constant.anc")
load("raw.anc")
load("morph.anc")
eco.anc <- mode.anc
# eco.anc <- constant.anc
# eco.anc <- raw.anc

# Import time trees saved from 1-MakeTimeTrees.R
load("~/Manuscripts/CamOrdEchino/equal.tree")
tree <- equal.tree

# Get epoch boundaries
strat_names <-
  read.csv("https://www.paleobiodb.org/data1.2/intervals/list.csv?all_records&vocab=pbdb")
# Eons are level 1, eras=level 2, periods=3, subperiods=4, epochs=5
TimeBins <- strat_names[which(strat_names$scale_level == 5),]
# Restrict to Cambrian and Ordovician:
TimeBins <- TimeBins$max_ma[which(TimeBins$max_ma > 443)]

# Using a modification of Claddis::DiscreteCharacterRate, renamed test_rates()
# starting v. 0.6.0, to allow direct input of the data matrix [from
# ReadMorphNexus()] with ancestral states previously reconstructed [from
# AncStateEstMatrix()], and adding both blocks into the same input object. Also
# includes other updated functions needed to run analyses below. Thanks, Graeme!
source("~/Manuscripts/CamOrdEchino/new_Claddis_functions.R")


## Process matrices into the arg names required by aug 2020 Claddis functions:

# Combine into single data matrix with two blocks
data <- list(topper = morph.anc$Topper, matrix_1 = morph.anc$Matrix_1, 
             matrix_2 = eco.anc$Matrix_1)
# Change names to match expected names in new Claddis functions
new.names <- c("block_names", "datatype", "matrix", "ordering", "character_weights", 
               "minimum_values", "maximum_values", "characters")
names(data$matrix_1) <- new.names
names(data$matrix_2) <- new.names

# Specify which characters are from the morphological block. The remaining
# (ecological) characters are automatically placed in the second block by the
# function.
morph.char <- 1:ncol(data$matrix_1$matrix)
eco.char <- max(morph.char) + (1:ncol(data$matrix_2$matrix))
character_partitions <- list(list(morphology = morph.char, ecology = eco.char))





## TEST RATE HETEROGENEITY BETWEEN ECOLOGY AND MORPHOLOGY ######################
# Using discrete rate test of Lloyd, et al. (2012) that uses a likelihood-ratio
# test. We also use the AIC method. Test compares the morphological and
# ecological data sets.

# Likelihood ratio tests:
LRT.rates <- 
  test_rates2(cladistic_matrix = data, time_tree = tree, time_bins = TimeBins,
              character_partitions = character_partitions, 
              polymorphism_state = "missing", uncertainty_state = "missing",
              inapplicable_state = "missing", time_binning_approach = "lloyd",
              test_type = "lrt", alpha = 0.01,
              multiple_comparison_correction = "benjaminihochberg")
beep()
# View results
round(LRT.rates$character_test_results[[1]]$rates, 4)
LRT.rates$character_test_results[[1]]$p_value

# Morphology rate = 1.3766 character changes / lineage Myr
# Ecology    rate = 0.7078 character changes / lineage Myr   p-value < 1e-66 ****
#                   0.6180 using 'constant'                  p < 1e-83 ****
#                   0.7690 using 'raw'                       p < 1e-27 ****



# AIC: 

# Use different partitions b/c comparing one-rate model to two-rate model
character_partitions.aic <-
  list(list(c(morph.char, eco.char)), list(morph.char))
AIC.rates <- 
  test_rates2(cladistic_matrix = data, time_tree = tree, time_bins = TimeBins,
              character_partitions = character_partitions.aic,
              polymorphism_state = "missing", uncertainty_state = "missing",
              inapplicable_state = "missing", time_binning_approach = "lloyd",
              test_type = "aic")
beep()
# View lambdas
round(AIC.rates$character_test_results[[1]]$rates, 4) # One-rate lambda
round(AIC.rates$character_test_results[[2]]$rates, 4) # Two-rate lambdas

# Convert AICc to akaike weights:
geiger::aicw(c(AIC.rates$character_test_results[[1]]$aicc,
               AIC.rates$character_test_results[[2]]$aicc))
# Given large number of characters, fine to use AIC, but reporting AICc anyway

# Mode:
# Global (overall) rate = 1.1831 character changes / lineage Myr
#      Morphology  rate = 1.3766 character changes / lineage Myr
#         Ecology  rate = 0.7078 character changes / lineage Myr
#
# AIC results (mode):       AICc      deltaAIC   Akaike weight      
# Model 1 (Single-rate)   7196.211      298.613   0.000
# Model 2 (Diff-rates)    6897.598        0.000   1.000 ***

# Constant:
# Global (overall) rate = 1.1710 character changes / lineage Myr
#      Morphology  rate = 1.3766 character changes / lineage Myr
#          Ecology rate = 0.6180 character changes / lineage Myr
#
# AIC results:              AICc        deltaAIC   Akaike weight
# Model 1 (Single-rate)   7205.563      374.802    0.000
# Model 2 (Diff-rates)    6830.761        0.000    1.000 ***

# Raw:
# Global (overall) rate = 1.2860 character changes / lineage Myr
#      Morphology  rate = 1.3766 character changes / lineage Myr
#          Ecology rate = 0.7690 character changes / lineage Myr
#
# AIC results:              AICc        deltaAIC   Akaike weight
# Model 1 (Single-rate)   6396.430      119.376    0.000
# Model 2 (Diff-rates)    6277.053        0.000    1.000 ***

# Visualize character partitions
plot_rates_character(test_rates_output = AIC.rates, model_number = 2)






## SUBSAMPLING ALGORITHM FOR CHARACTER PARTITIONS ##############################
# The lambda rates returned above are sensitive to the number of characters.
# (All things equal, more characters in a partition can produce an inflated
# rate.) To avoid this potential bias, here we subsample each partition to a
# common number of characters (and removing any characters with invariant rate).

# Removing following invariant (rate = 0) characters from ecological matrix: 7,
# 11, 14, 33, 35, 37 dealing with non-sexual reproduction, autotrophy,
# herbivory, and swimming.
eco.rates <- read.csv(file = "EcoRates.csv")     # Created below
(invar.eco <- which(eco.rates$rate == 0))
length(invar.eco) / nrow(eco.rates)              # 15% is invariant

# Removing following invariant (rate = 0) characters from morphogical matrix: 
morph.rates <- read.csv(file = "MorphRates.csv") # Created below
(invar.morph <- which(morph.rates$rate == 0))
length(invar.morph) / nrow(morph.rates)          # 38% is invariant

# Need to make sequential (with morphological characters first)
invars <- c(invar.morph, invar.eco + 413)

# Remove the invariants
no.invars <- prune_cladistic_matrix(cladistic_matrix = data, characters2prune = invars)
n.morph.char <- ncol(no.invars$matrix_1$matrix)  # 255 variable characters
n.eco.char <- ncol(no.invars$matrix_2$matrix)    #  34 variable characters
total.chars <- n.morph.char + n.eco.char

n.char <- 30                                     # Subsample to 30 characters
# Specify the subsampled partitions
character_partitions.aic <- list(list(1:(n.char * 2)), list(1:n.char))

## Resample in parallel

# Initialize cluster
nreps <- 1000
cl <- makeCluster(detectCores())
registerDoParallel(cl)
clusterSetRNGStream(cl, 3142)  # Set L-Ecuyer RNG seed
(start <- Sys.time())

# Run in parallel
aic.results <- foreach(i = 1:nreps, .inorder = FALSE, 
                       .packages = "Claddis") %dopar% {

  # Choose subsampled characters and prune cladistic matrix accordingly
  wh.morph.chars <- sample(1:n.morph.char, n.char, replace = FALSE)
  wh.eco.chars <- sample(1:n.eco.char, n.char, replace = FALSE)
  characters2prune <-
    setdiff(1:total.chars, c(wh.morph.chars, wh.eco.chars + n.morph.char))
  subsample <- prune_cladistic_matrix(cladistic_matrix = no.invars, 
                                      characters2prune = characters2prune)
  AIC.rates <- 
    test_rates2(cladistic_matrix = subsample, time_tree = tree, time_bins = TimeBins,
                character_partitions = character_partitions.aic,
                polymorphism_state = "missing", uncertainty_state = "missing",
                inapplicable_state = "missing", time_binning_approach = "lloyd",
                test_type = "aic")
}

stopCluster(cl)
(Sys.time() - start)     # 31.4 minutes on 8-core laptop
beep(3)

# Save output
# save(aic.results, file = "aic.results")
# load("aic.results")

# Convert into clean data table
l.aic <- length(aic.results)
resampled.aic <- 
  data.frame(rate1 = NA, AIC1 = NA, AICc1 = NA, rate.morph = NA, rate.eco = NA, 
             AIC2 = NA, AICc2 = NA, AW1 = NA, AW2 = NA)
for(r in 1:l.aic) {
  resampled.aic[r, 1:7] <-
    as.numeric(unlist(aic.results[[r]]$character_test_results)[c(1:3, 5:8)])
  resampled.aic[r, 8:9] <-
    geiger::aicw(c(resampled.aic$AICc1[r], resampled.aic$AICc2[r]))$w
}
head(resampled.aic)
round(apply(resampled.aic, 2, mean), 3)
round(apply(resampled.aic, 2, sd), 3)
round(apply(resampled.aic, 2, median), 3)

# pdf(file = "SubsampledRates.pdf")
breaks <- pretty(c(resampled.aic$rate.morph, resampled.aic$rate.eco), 40)
hist(resampled.aic$rate.eco, main = "Subsampled rates", 
     xlab = "Rate (character changes / lineage Myr)", ylab = "Density", 
     breaks = breaks, col = "transparent", border = "transparent", prob = TRUE)
hist(resampled.aic$rate.eco, add = TRUE, border = "white", col = "darkgray", 
     breaks = breaks, prob = TRUE)
hist(resampled.aic$rate.morph, add = TRUE, border = "black", col = "transparent", 
     breaks = breaks, prob = TRUE)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
par(op)
# dev.off()


# pdf(file = "TwoModelAkaikeWeight.pdf")
breaks <- seq(0.25, 1, by = 0.01)
hist(resampled.aic$AW2, main = "Akaike weight for two-rate model", 
     xlab = "Akaike weight", ylab = "Density", breaks = breaks, prob = TRUE, 
     col = "black")
abline(v = 0.99, lty = 2, lwd = 2, col = "red")
par(op)
# dev.off()

summary(resampled.aic$AW2) # mean weight = 0.955, median = 1.000

# What proportion >= 0.9?
100 * length(resampled.aic$AW2[resampled.aic$AW2 >= 0.9]) / nrow(resampled.aic)
# 91.2% greater than 0.90

100 * length(resampled.aic$AW2[resampled.aic$AW2 >= 0.999]) / nrow(resampled.aic)
# 85.2% ~ = 1




## TABULATE BRANCHING DISTANCE & NUMBER OF CHARACTER CHANGES PER NODE ##########
# Import distance matrices from 3-DisparityDistances.R
load("mode.distances.GED.5"); dist.matrix <- mode.distances.GED.5$DistanceMatrix; load("mode.anc"); anc.matrix <- mode.anc$Matrix_1$Matrix
# load("constant.distances.GED.5"); dist.matrix <- constant.distances.GED.5$DistanceMatrix; load("constant.anc"); anc.matrix <- constant.anc$Matrix_1$Matrix
# load("raw.distances.GED.5"); dist.matrix <- raw.distances.GED.5$DistanceMatrix; load("raw.anc"); anc.matrix <- raw.anc$Matrix_1$Matrix
load("morph.distances.GED.5"); dist.matrix <- morph.distances.GED.5$DistanceMatrix; load("morph.anc"); anc.matrix <- morph.anc$Matrix_1$Matrix

# Function to identify number of characters that are identical, including proper
# handling of poymorphisms and missing states (where they match to the paired
# states, if matchable). Influenced by code by Graeme Lloyd in 'Claddis'.
better.identical <- function(a, b) {
  if (!identical(dim(a), dim(b)))
    stop("data frames have different sizes\n")
  a <- matrix(a, nrow = 1)
  b <- matrix(b, nrow = 1)
  colnames(a) <- colnames(b) <- 1:length(a)
  # Convert empty cells to NAs
  a[which(a == "" | a == " ")] <- NA
  b[which(b == "" | b == " ")] <- NA
  # Convert polymorphisms and NAs to tip state (if match-able, leaving unchanged
  # if not)
  PolymorphismAndUncertaintyPositions.a <- 
    sort(unique(c(grep("/|&", a), which(is.na(a)))))
  PolymorphismAndUncertaintyPositions.b <- 
    sort(unique(c(grep("/|&", b), which(is.na(b)))))
  if (length(PolymorphismAndUncertaintyPositions.a) > 0L) {
    # polys <- strsplit(b[PolymorphismAndUncertaintyPositions.a], "/")
    for (ch in PolymorphismAndUncertaintyPositions.a) {
      # First select matched polymorphism:
      opts <- strsplit(a[ch], "/|&")
      if (b[ch] %in% opts[[1]])
        a[ch] <- b[ch]
      # Then set NA:
      if (is.na(a[ch]))
        a[ch] <- b[ch]
    }
  }
  if (length(PolymorphismAndUncertaintyPositions.b) > 0L) {
    # polys <- strsplit(b[PolymorphismAndUncertaintyPositions.b], "/")
    for (ch in PolymorphismAndUncertaintyPositions.b) {
      opts <- strsplit(b[ch], "/|&")
      if (a[ch] %in% opts[[1]])
        b[ch] <- a[ch]
      if (is.na(b[ch]))
        b[ch] <- a[ch]
    }
  }
  wh.diff <- which(!mapply(identical, a, b))
  ab <- rbind(a, b)
  return(ab[, wh.diff, drop = FALSE])
}

# Function to calculate branch distances and number of character changes for all
# edges in a tree.
tally.branch.changes <- function(dist.matrix = NULL, anc.matrix = NULL,
                                 tree = NULL) {
  out <- data.frame(Node = tree$edge[, 1], Tip = tree$edge[, 2],
                    Branch.dist = NA, Char.changes = NA)
  for (n in 1:nrow(out)) {
    out$Branch.dist[n] <- dist.matrix[out$Node[n], out$Tip[n]]
    out$Char.changes[n] <- ncol(better.identical(anc.matrix[out$Node[n], ], 
                                                 anc.matrix[out$Tip[n], ]))
  }
  return(out)
}


eco.branch.changes <- 
  tally.branch.changes(tree = tree, anc.matrix = mode.anc$Matrix_1$Matrix,
                       dist.matrix = mode.distances.GED.5$DistanceMatrix)
plot(eco.branch.changes$Branch.dist, eco.branch.changes$Char.changes)
cor(eco.branch.changes$Branch.dist, eco.branch.changes$Char.changes)

morph.branch.changes <- 
  tally.branch.changes(tree = tree, anc.matrix = morph.anc$Matrix_1$Matrix,
                       dist.matrix = morph.distances.GED.5$DistanceMatrix)
plot(morph.branch.changes$Branch.dist, morph.branch.changes$Char.changes)
cor(morph.branch.changes$Branch.dist, morph.branch.changes$Char.changes)

# Save objects
# save(eco.branch.changes, file = "eco.branch.changes")
# save(morph.branch.changes, file = "morph.branch.changes")
# load("eco.branch.changes")
# load("morph.branch.changes")


# pdf(file = "PerBranchChanges.pdf")
par(mfrow = c(2, 2))
hist(morph.branch.changes$Branch.dist, 20, main = "Morphology", xlab = "Branching distance")
hist(eco.branch.changes$Branch.dist, 20, main = "Ecology", xlab = "Branching distance")
hist(100 * morph.branch.changes$Char.changes / 413, 20, main = "Morphology", xlab = "% character changes")
hist(100 * eco.branch.changes$Char.changes / 40, 20, main = "Ecology", xlab = "% character changes")
# dev.off()
par(op)

# pdf(file = "Dist&CharPerBranch.pdf")
par(mfrow = c(1, 2), mar = c(4, 4, 1, 0))
breaks.dist <- 
  pretty(c(morph.branch.changes$Branch.dist, eco.branch.changes$Branch.dist), 10)
breaks.char <- pretty(c(morph.branch.changes$Char.changes, 
                        eco.branch.changes$Char.changes), 10)
hist(eco.branch.changes$Branch.dist, 
     main = "Distance moved / branching event", xlab = "Distance moved", 
     ylab = "# branching events", breaks = breaks.dist, col = "transparent", 
     border = "transparent", cex.lab = 1, cex.main = 0.8)
hist(eco.branch.changes$Branch.dist, add = TRUE, border = "white", 
     col = "darkgray", breaks = breaks.dist)
hist(morph.branch.changes$Branch.dist, add = TRUE, border = "black", 
     col = "transparent", breaks = breaks.dist)
hist(eco.branch.changes$Branch.dist, xlab = "# character changes", 
     main = "# character changes / branching event", 
     ylab = "# branching events", breaks = breaks.char, col = "transparent", 
     border = "transparent", cex.lab = 1, cex.main = 0.8)
hist(eco.branch.changes$Char.changes, add = TRUE, border = "white", 
     col = "darkgray", breaks = breaks.char)
hist(morph.branch.changes$Char.changes, add = TRUE, border = "black", 
     col = "transparent", breaks = breaks.char)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
par(op)
# dev.off()

# Statistical tests
summary(morph.branch.changes$Branch.dist)
summary(eco.branch.changes$Branch.dist)

summary(morph.branch.changes$Char.changes)
summary(eco.branch.changes$Char.changes)

summary(100 * morph.branch.changes$Char.changes / 413)
summary(100 * eco.branch.changes$Char.changes / 40)

wilcox.test(morph.branch.changes$Branch.dist, eco.branch.changes$Branch.dist)
# W = 346033, p < 2.2e-16

wilcox.test(morph.branch.changes$Char.changes, eco.branch.changes$Char.changes)
# W = 417274, p < 2.2e-16

wilcox.test(morph.branch.changes$Char.changes / 413, eco.branch.changes$Char.changes / 40)
# W = 299869, p = 1.794e-5

# Morphological characters change more often than ecological characters (but
# proportionally less), and when they do, the change is a more substantial one.







## WHICH CHARACTERS EVOLVE FASTEST AND SLOWEST? ################################

# Tally ecological characters
morph.char <- 1:ncol(data$matrix_1$matrix)
eco.char <- max(morph.char) + (1:ncol(data$matrix_2$matrix))
nchar <- max(eco.char)
character_partitions <- lapply(as.list(1:nchar), as.list)

all.char.partitions <- 
  test_rates2(cladistic_matrix = data, time_tree = tree, time_bins = TimeBins,
              character_partitions = character_partitions, 
              polymorphism_state = "missing", uncertainty_state = "missing",
              inapplicable_state = "missing", time_binning_approach = "lloyd",
              test_type = "lrt", alpha = 0.01,
              multiple_comparison_correction = "benjaminihochberg")
beep()

# Convert to cleaner output
rates <- data.frame(character = 1:nchar, rate = NA, LRT.pval = NA, 
                        sig.fast.or.slow = NA)
rates$rate <- 
  round(sapply(1:nchar, function(x) all.char.partitions$character_test_results[[x]]$rates[1]), 5)
# prior line same as  > all.char.partitions$character_rates[, 2]
rates$LRT.pval <- 
  sapply(1:nchar, function(x) all.char.partitions$character_test_results[[x]]$p_value)
signs <- 
  sign(sapply(1:nchar, function(x) all.char.partitions$character_test_results[[x]]$rates[1] - all.char.partitions$character_test_results[[x]]$rates[2]))
adj.alphas <- 
  sapply(1:nchar, function(x) all.char.partitions$character_test_results[[x]]$CorrectedAlpha)
for(ch in 1:nchar) {
  if (rates$LRT.pval[ch] < adj.alphas[ch]) {
    if (signs[ch] > 0)
      rates$sig.fast.or.slow[ch] <- "+"
    if (signs[ch] < 0)
      rates$sig.fast.or.slow[ch] <- "-"
  }
  else
    rates$sig.fast.or.slow[ch] <- ""
}

# Append number of character changes, completeness & duration
rates <- cbind(rates, as.data.frame(all.char.partitions$character_rates[, 3:5]))
beep()

# Visualize and save (+ = significantly faster than other characters, 
#                     - = significantly slower)
head(rates)

morph.rates <- rates[morph.char, ]
eco.rates <- rates[eco.char, ]
# Renumber eco.rates for consistency with raw data
eco.rates$character <- 1:length(eco.char)
rownames(eco.rates) <- as.character(1:length(eco.char))

# write.csv(morph.rates, file = "MorphRates.csv", row.names = FALSE)
# write.csv(eco.rates, file = "EcoRates.csv", row.names = FALSE)
# morph.rates <- read.csv(file = "MorphRates.csv")
# eco.rates <- read.csv(file = "EcoRates.csv")

# 10 fastest evolving characters (* = significantly so):
eco.rates[order(eco.rates$rate, decreasing = TRUE)[1:10], ]
# In order, *Body size, *RelFoodStratification, *AbsFoodStratification,
# *FilterDensity, *AbsStratification, *RelStratification, HardSubstratum,
# SoftSubstratum, and BioticSubstrate.

# 10 slowest evolving characters (* = significantly so):
eco.rates[order(eco.rates$rate, decreasing = FALSE)[1:12], ]
# In order, 6 characters with no variation (*sexual reproduction,
# *FluidicSubstrate & *InsubstantialSubstrate (both = swimming), *autotrophy,
# *herbivory, and *incorporeal (= autotrophy). Those that vary include (in
# order): *living above substrate, *ambient feeding (= osmotrophy and feeding on
# dissolved organic matter), *attachment feeding, *asexual reproduction, and
# *solution feeding. Moderately slow characters (also significantly slower)
# include *AboveImmediate, *WithinImmediate, *FilterFeeding, *SelfSupported,
# *Supported, *MassFeeding, and *Microbivory, followed by (still slightly
# faster) *FeedingWithinPrimary (= deposit feeding), *FeedingWithinImmediate,
# *RaptorialFeeding, *ParticleFeeding, *FeedingAboveImm, *Carnivory,
# *BulkFeeding, *Attached, *Free-living, *FeedingAbovePrimary, *LithicSubstrate,
# and *WithnnPrimary (= infaunal). Except for supported/self-suppported (having
# to do with whether relies on another organism to inhabit microhabitat) and
# whether live epifaunally vs. infaunally, remaining 7 significantly slow
# characters all deal with feeding (diet, foraging, and food location).

# In other words, among the 24 significantly slowest (varying) characters, 14
# have to do with diet/foraging/food location. Remainder have to do with
# reproduction, infaunal/epifaunal, and whether the animal lives atop another
# organism to inhabit its microhabitat. Among the fastest evolving, they involve
# body size, tiering, and filter density. (These 6 characters also have more
# states and are ordered factors or numerics.)


# How many significantly slower, faster, and not sig. diff?
table(eco.rates$sig.fast.or.slow)
# 6 significantly fast and 30 significantly slow

hist(eco.rates$rate, 20, main = "Ecological character rate", 
     xlab = "Rate (character changes / Myr)")




# 10 fastest evolving characters (* = significantly so):
morph.rates[order(morph.rates$rate, decreasing = TRUE)[1:10], ]
# The fastest evolving characters include (in order) *261 (rigid oral surface
# [tegmen] height [height/width]), 61 (oral 1 plate shape), *299 (plate size in
# circlet at top of the theca), 237 (periproct size), 92 (width of recumbant
# ambulacra), *341 (number of series of plate in stele), *338 (number of
# marginal appendages), *263 (tegmen plate scupturing), *298 (shape of plates in
# the circlet at the top of the theca, regularity of plate shape), *304 (circlet
# of plates below top circlet of plates, if not bottom circlet), *315 (interrays
# depressed), 270 (), *266 (periproct position, including position of tube), and
# *55 (number of oral plates).

# 10 slowest evolving characters (only listing those that are significantly so):
morph.rates[order(morph.rates$rate, decreasing = FALSE)[1:10], ]
# The slowest evolving characters include (in order) 110 (lancet plates), 159
# (proximal ambulacral cover plates), 231 (genital plates with single
# perforation surrounding periproct), 233 (periproct distinctive opening), 360
# (branching on column proxistele), all characters with very little variation
# among taxa (often only 1 genus with a state). The next 10 significantly slow
# characters include (in order) 7 (direction of growth), 17 (oral region
# surrounded by ctenoid plates), 50 (calcareous ring surrounding pharynx), 19
# (extra row of plates along oral area on ventral surface = lateral plates in
# ctenocystoids), 386 (Aboral covering mimics ambulacral pattern, as in
# asterozoans), 20 (food groove on basal side of cover plates, as in
# Ctenocystis), 75 (ambulacra), 373 (distinctive thecal appendage not bearing
# food groove = stele in Homoiostelea), 201 (appendage bearing peristome), and
# 229 (Calcified hydropore = madreporite).

# How many significantly slower, faster, and not sig. diff?
table(morph.rates$sig.fast.or.slow)
# 65 significantly fast and 54 significantly slow

hist(eco.rates$rate, 20, main = "Ecological character rate", 
     xlab = "Rate (character changes / Myr)")
hist(morph.rates$rate, 20, main = "Morphological character rate", 
     xlab = "Rate (character changes / Myr)")

summary(morph.rates$rate)
summary(eco.rates$rate)

## *** NOTE IT IS INAPPROPRIATE TO DO SIMPLE TWO-SAMPLE TESTS COMPARING RATES OF
## ECOLOGY VS. MORPHOLOGY. See Lloyd, et al. (2012) and Claddis::test_rates()
## for proper tests. ***

wilcox.test(morph.rates$rate, eco.rates$rate)
# W = 8448.5, p-value = 0.8075 # Not sig. diff on a per-character basis

# Graeme Lloyd claims this is not an appropriate way to test this comparison
# because it ignores differences in character branch completeness and branch
# duration. The Claddis::test_rates() partition test above is the preferable way
# to model these additional factors, and strongly supports that the rate of
# morphology is greater than of ecology. Further two-sample tests are not
# conducted here given inappropriateness.

# pdf(file = "CharacterRates.pdf")
breaks <- pretty(c(morph.rates$rate, eco.rates$rate), 20)
hist(eco.rates$rate, main = "Per-character rates", 
     xlab = "Rate (per-character changes / Myr)", 
     ylab = "Proportion of characters", breaks = breaks, col = "transparent", 
     border = "transparent", prob = TRUE)
hist(eco.rates$rate, add = TRUE, border = "white", col = "darkgray", 
     breaks = breaks, prob = TRUE)
hist(morph.rates$rate, add = TRUE, border = "black", col = "transparent", 
     breaks = breaks, prob = TRUE)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
par(op)
# dev.off()





## Do primary vs dependent (secondary, etc.) characters evolve at different
## rates?

# Let's hide creating this ugliness with brackets:
{
  morph.level <- 
    c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 2, 1, 1, 1, 2, 2, 
      2, 2, 2, 1, 2, 2, 3, 3, 4, 3, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 1, 
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 
      2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
      3, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 4, 4, 
      4, 1, 2, 2, 3, 4, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 4, 2, 1, 2, 2, 2, 1, 2, 2, 
      2, 2, 2, 2, 2, 2, 2, 3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 4, 2, 3, 3, 3, 
      3, 3, 3, 3, 3, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 3, 3, 3, 4, 4, 
      2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 1, 2, 2, 3, 3, 3, 
      3, 3, 4, 1, 2, 1, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 2, 3, 
      3, 3, 3, 3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 1, 2, 1, 2, 
      2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
      3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 4, 5, 5, 5, 3, 4, 4, 2, 3, 3, 3, 2, 
      3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 5, 5, 4, 5, 5, 5, 5, 5, 5, 1, 2, 2, 2, 
      2, 2, 3, 3, 3, 4, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 2, 3, 3, 
      3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 
      4, 3, 4, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2)
  morph.nstates <- 
    c(3, 3, 2, 6, 2, 3, 2, 2, 3, 2, 3, 2, 3, 4, 3, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 
      2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 1, 2, 2, 2, 
      2, 2, 2, 2, 6, 2, 2, 2, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 2, 4, 3, 4, 2, 2, 2, 
      2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 4, 4, 2, 3, 3, 2, 2, 2, 2, 
      2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
      2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 2, 3, 2, 2, 4, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1,
      1, 1, 1, 1, 2, 2, 3, 5, 4, 2, 3, 3, 2, 2, 2, 2, 8, 2, 7, 2, 2, 4, 2, 3, 3,
      2, 2, 2, 3, 3, 2, 2, 1, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 5, 4, 2, 2, 2,
      2, 2, 2, 2, 1, 2, 2, 2, 3, 2, 3, 3, 2, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 3,
      1, 2, 5, 2, 2, 2, 2, 2, 2, 2, 4, 2, 5, 2, 1, 3, 2, 1, 2, 3, 3, 2, 2, 2, 2,
      5, 3, 4, 4, 5, 2, 5, 12, 2, 2, 2, 3, 3, 5, 5, 7, 5, 5, 3, 6, 5, 7, 2, 3, 2,
      3, 3, 2, 2, 3, 2, 2, 4, 4, 5, 2, 5, 2, 2, 2, 2, 3, 2, 4, 3, 2, 1, 2, 2, 2,
      2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 3, 2,
      3, 5, 2, 4, 3, 2, 6, 3, 2, 2, 2, 7, 4, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2,
      2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2)
  eco.level <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  eco.nstates <- c(7, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
                   2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
}
table(morph.level, morph.nstates)
table(eco.level, eco.nstates)

boxplot(morph.rates$rate ~ morph.level, 
        main = "Variation in morphological character rates by level")
aov <- aov(morph.rates$rate ~ morph.level)
summary(aov) # F = 11.74, p = 0.00067

# M-W test by primary vs. any dependent
which.prim.morph <- which(morph.level == 1L)
summary(morph.rates$rate[which.prim.morph])
summary(morph.rates$rate[-which.prim.morph])
# Primary characters DO have slightly slower rate of evolution

wilcox.test(morph.rates$rate[which.prim.morph], morph.rates$rate[-which.prim.morph])
# W = 8333, p = 0.769: But primary not significantly different than dependents.

# Confirm it is only the quinary characters that are different:
which.5.morph <- which(morph.level == 5L)
wilcox.test(morph.rates$rate[-which.5.morph], morph.rates$rate[which.5.morph])
# W = 1255, p = 0.0011: Only the quinary characters evolve at faster rate than
# others.



boxplot(eco.rates$rate ~ eco.level, 
        main = "Variation in ecological character rates by level")
aov <- aov(eco.rates$rate ~ eco.level)
summary(aov) # F = 0.172, p = 0.172

# M-W test by primary vs. any dependent
which.prim.eco <- which(eco.level == 1L)
summary(eco.rates$rate[which.prim.eco])
summary(eco.rates$rate[-which.prim.eco])
# Primary characters DO have slightly slower rate of evolution (but only 1
# dependent = filter density)

wilcox.test(eco.rates$rate[which.prim.eco], eco.rates$rate[-which.prim.eco])
# W = 3, p = 0.165: But primary not significantly different than dependents.





# pdf(file = "PrimaryCharacterRates.pdf")
par(mfrow = c(2,2), mar = c(4.5, 4, 2, 0.1))
boxplot(morph.rates$rate ~ morph.level, horizontal = TRUE, main = "morphology", 
        xlab = "Rate (per-char. changes / Myr)", ylab = "Dependency level")
boxplot(eco.rates$rate ~ eco.level, horizontal = TRUE, main = "ecology", 
        xlab = "Rate (per-char. changes / Myr)", ylab = "Dependency level")
breaks <- pretty(morph.rates$rate, 20)
hist(morph.rates$rate, main = "", 
     xlab = "Rate (per-char. changes / Myr)", 
     ylab = "Proportion of characters", breaks = breaks, col = "transparent", 
     border = "transparent", prob = TRUE)
hist(morph.rates$rate[which.prim.morph], add = TRUE, border = "white", 
     col = "darkgray", breaks = breaks, prob = TRUE)
hist(morph.rates$rate[-which.prim.morph], add = TRUE, border = "black", 
     col = "transparent", breaks = breaks, prob = TRUE)

breaks <- pretty(eco.rates$rate, 20)
hist(eco.rates$rate, main = "", 
     xlab = "Rate (per-char. changes / Myr)", ylab = "# of characters", 
     breaks = breaks, col = "transparent", border = "transparent")
hist(eco.rates$rate[which.prim.eco], add = TRUE, border = "white", 
     col = "darkgray", breaks = breaks)
hist(eco.rates$rate[-which.prim.eco], add = TRUE, border = "black", 
     col = "transparent", breaks = breaks)
legend("topright", inset = .05, c("primary", "dependent"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)
par(op)
# dev.off()




## Side issue: Is the delta statistic correlated with character rate? ##########
load("deltas.morph")
load("deltas.eco")

# Ecology
plot(log(eco.rates$rate), log(deltas.eco), main = "Ecology") # Appears to be a power function
wh.zero <- which(eco.rates$rate == 0)
lm.eco <- lm(log(deltas.eco[-wh.zero]) ~ log(eco.rates$rate[-wh.zero]))
abline(lm.eco)
summary(lm.eco) # adj r2 = 0.686 with inverse correlation, p < 1e-9

# Morphology
plot(log(morph.rates$rate), log(deltas.morph), main = "Morphogy") # Appears to be a power function
wh.zero <- which(morph.rates$rate == 0)
lm.morph <- lm(log(deltas.morph[-wh.zero]) ~ log(morph.rates$rate[-wh.zero]))
abline(lm.morph)
summary(lm.morph) # adj r2 = 0.079 with weak but significant inverse correlation, p < 1e-5

# The delta statistic of Borges, et al. (2019) is significantly inversely
# correlated as a power function with the character rate metrics of Lloyd, et
# al. (2012). However, they are not identical so can not be converted to one
# another in a trivial manner.




## RATE TRENDS THROUGH TIME IN TWO DATA BLOCKS (& AMONG BRANCHES & CLADES) #####

# Here we calculate the rate of character evolution in each block separately.
# Because we are only interested in the rate and not in modeling rate
# partitions, we simplify the test. We also remove invariant characters.

# Isolate each block
morph.data <- prune_cladistic_matrix(cladistic_matrix = data, blocks2prune = 2,
                                     remove_invariant = TRUE)
eco.data <- prune_cladistic_matrix(cladistic_matrix = data, blocks2prune = 1,
                                   remove_invariant = TRUE)

# Choose partitions
tbins <- seq(from = min(TimeBins), to = max(TimeBins), length.out = 50)
mean(diff(tbins)) # 2.0 Myr bins
time.partitions <- list(c(list(1), list(2:(length(tbins) - 1))))
branch.part <- list(c(list(1)))
branch.part <- lapply(X = as.list(seq(1, length(tree$edge.length))), as.list)
clade.part <- as.list(seq(ape::Ntip(tree) + 1, ape::Ntip(tree) + ape::Nnode(tree)))

morph.rates <- test_rates2(cladistic_matrix = morph.data, time_tree = tree, 
                           time_bins = tbins, branch_partitions = branch.part,
                           time_partitions = time.partitions,
                           clade_partitions = clade.part)
eco.rates <- test_rates2(cladistic_matrix = eco.data, time_tree = tree,
                         time_bins = tbins, branch_partitions = branch.part,
                         time_partitions = time.partitions,
                         clade_partitions = clade.part)
beep()

# Save output
# save(morph.rates, file = "morph.rates")
# save(eco.rates, file = "eco.rates")
# load("morph.rates")
# load("eco.rates")


# Visualize rates through time (ignoring the partitions)
par(op)
plot_rates_time(test_rates_output = morph.rates, model_number = 1)
mtext("morphology", side = 3)
plot_rates_time(test_rates_output = eco.rates, model_number = 1)
mtext("ecology", side = 3)
# Initially high rates (and last low rates) are because of very small bin
# duration and few number of branches in these "edges".


# Visualize clade partitions

# Find lowest AIC models for clade partitions
morph.aics <- sapply(1:length(morph.rates$clade_test_results), function(x) morph.rates$clade_test_results[[x]]$aic[1])
# 10 best models
data.frame(model = order(morph.aics)[1:10], AIC = round(morph.aics[order(morph.aics)[1:10]]))
(morph.clade.model <- which.min(morph.aics))
eco.aics <- sapply(1:length(eco.rates$clade_test_results), function(x) eco.rates$clade_test_results[[x]]$aic[1])
# 10 best models
data.frame(model = order(eco.aics)[1:10], AIC = round(eco.aics[order(eco.aics)[1:10]]))
(eco.clade.model <- which.min(eco.aics))

par(op)
# pdf("clade_rates_morphology.pdf")
plot_rates_tree(test_rates_output = morph.rates, model_type = "clade", 
                model_number = morph.clade.model, cex = 0.75)
legend("topright", inset = c(0.1, 0), legend = "morphology", bty = "n", cex = 1.5)
# dev.off()
# Best model has morphological rate slowing significantly among asterozoans,
# echinoids, holothuroids, and ophiocistoids. 2nd best model throws in unusual
# edrioasteroid Camptostroma. Most other viable models have higher rates at
# origins of blastozoans + crinozoans.

# pdf("clade_rates_ecology.pdf")
plot_rates_tree(test_rates_output = eco.rates, model_type = "clade", 
                model_number = eco.clade.model, cex = 0.75)
legend("topright", inset = c(0.1, 0), legend = "ecology", bty = "n", cex = 1.5)
# dev.off()

# Best model has ecological rate slowing significantly among homostelean
# cinctans. Other top models focus on higher rates among everything other than
# ctenocystids + cinctans, higher rate among small group of asteroids, higher
# rate among edrioasteroids + asteroids (plus ophios, holothuroids, echinoids, +
# cyclocystoids).

# Visualize edge rates
# pdf("branch_rates_morphology.pdf")
edge_rate_values <- log(morph.rates$branch_rates[, 2] + 1)
phytools::plotBranchbyTrait(tree = tree, x = edge_rate_values, mode = "edge", 
                            xlims = c(0, max(edge_rate_values)), 
                            title = "log(Changes per lineage myr + 1)", 
                            leg = max(nodeHeights(tree)), show.tip.label = FALSE)
legend("topright", inset = c(0.1, 0), legend = "morphology", bty = "n", cex = 1.5)
# dev.off()

# pdf("branch_rates_ecology.pdf")
edge_rate_values <- log(eco.rates$branch_rates[, 2] + 1)
phytools::plotBranchbyTrait(tree = tree, x = edge_rate_values, mode = "edge", 
                            xlims = c(0, max(edge_rate_values)), 
                            title = "log(Changes per lineage myr + 1)", 
                            leg = max(nodeHeights(tree)), show.tip.label = FALSE)
legend("topright", inset = c(0.1, 0), legend = "ecology", bty = "n", cex = 1.5)
# dev.off()


# Plot rates on same timescale and axis

# Modification of geoscale::geoscalePlot to allow ICS 2020 timescale and
# multipanel plots
source("~/Manuscripts/CamOrdEchino/geoscalePlot2.R")

# Import ICS 2020 timescale to use in plotting
ICS2020 <- read.csv("~/Manuscripts/CamOrdEchino/timescales2020.csv", 
                    stringsAsFactors =  TRUE)

time_bin_mids <- 
  (morph.rates$time_bins_used[2:length(x = morph.rates$time_bins_used)] 
   + morph.rates$time_bins_used[1:(length(x = morph.rates$time_bins_used) - 1)]) / 2
lim <- c(0, max(c(morph.rates$time_rates[, 2], eco.rates$time_rates[, 2])))

geoscalePlot2(tbins, rep(lim[1], length(tbins)), units = c("Epoch", "Period"), 
              tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
              cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
              ts.col = TRUE, label = "Character changes per lineage million years", 
              timescale = ICS2020, type = "n", abbrev = "Period")
mtext(text = "Character rate", side = 3, cex = 1.25)
cols <- viridisLite::plasma(3)         # Uses RGB-sensitive red and blue
lines(time_bin_mids, eco.rates$time_rates[, 2], lwd = 4, lty = 1, col = cols[1])
lines(time_bin_mids, morph.rates$time_rates[, 2], lwd = 4, lty = 6, col = cols[2])
legend("topright", legend = c("ecology", "morphology"), col = cols, 
       bty = "n", lty = c(1, 6), lwd = 4, inset = 0.02, cex = 1.5)


# Same, but as per-character rates (dividing by number of characters in each block)
# pdf(file = "rate_through_time.pdf")
lim <- c(0, max(c(morph.rates$time_rates[, 2] / 413, eco.rates$time_rates[, 2] / 40)))
geoscalePlot2(tbins, rep(lim[1], length(tbins)), units = c("Epoch", "Period"), 
              tick.scale = "Period", boxes = "Age", cex.age = 0.65, 
              cex.ts = 0.7, cex.pt = 1, age.lim = c(540, 445), data.lim = lim, 
              ts.col = TRUE, label = "Per-character changes per lineage million years", 
              timescale = ICS2020, type = "n", abbrev = "Period")
mtext(text = "Per-character rate", side = 3, cex = 1.25)
lines(time_bin_mids, eco.rates$time_rates[, 2] / 40, lwd = 4, lty = 1, col = cols[1])
lines(time_bin_mids, morph.rates$time_rates[, 2] / 413, lwd = 4, lty = 6, col = cols[2])
legend("topright", legend = c("ecology", "morphology"), col = cols, 
       bty = "n", lty = c(1, 6), lwd = 4, inset = 0.02, cex = 1.5)
# dev.off()


# Summarize rates:
summary(eco.rates$time_rates[, 2] / 40)
summary(morph.rates$time_rates[, 2] / 413)

hist(eco.rates$time_rates[, 2] / 40, 20)
hist(morph.rates$time_rates[, 2] / 413, 20)

# Both on same histogram
breaks <-
  pretty(c(eco.rates$time_rates[, 2] / 40, morph.rates$time_rates[, 2] / 413), 40)
hist(eco.rates$time_rates[, 2] / 40, main = "Per-character rates", 
     xlab = "Rate (per-character changes / lineage Myr)", ylab = "Density", 
     breaks = breaks, col = "transparent", border = "transparent", prob = TRUE)
hist(eco.rates$time_rates[, 2] / 40, add = TRUE, border = "white", col = "darkgray", 
     breaks = breaks, prob = TRUE)
hist(morph.rates$time_rates[, 2] / 413, add = TRUE, border = "black", col = "transparent", 
     breaks = breaks, prob = TRUE)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)

# Same, without first two intervals
breaks <-
  pretty(c(eco.rates$time_rates[, 2][-(1:2)] / 40, morph.rates$time_rates[, 2][-(1:2)] / 413), 40)
hist(eco.rates$time_rates[, 2][-(1:2)] / 40, main = "Per-character rates", 
     xlab = "Rate (per-character changes / lineage Myr)", ylab = "Density", 
     breaks = breaks, col = "transparent", border = "transparent", prob = TRUE)
hist(eco.rates$time_rates[, 2][-(1:2)] / 40, add = TRUE, border = "white", col = "darkgray", 
     breaks = breaks, prob = TRUE)
hist(morph.rates$time_rates[, 2][-(1:2)] / 413, add = TRUE, border = "black", col = "transparent", 
     breaks = breaks, prob = TRUE)
legend("topright", inset = .05, c("ecology", "morphology"), pch = c(22, 22), 
       pt.bg = c("darkgray", "transparent"), col = c("darkgray", "black"), 
       cex = 1, pt.cex = 2)

# How much greater? (2.14 X)
mean(morph.rates$time_rates[, 2] / 413) / mean(eco.rates$time_rates[, 2] / 40)

# How much greater (excluding 1st or first two interval)?
# (2.41 X)
mean(morph.rates$time_rates[, 2][-1] / 413) / mean(eco.rates$time_rates[, 2][-1] / 40)
# (2.22 X)
mean(morph.rates$time_rates[, 2][-(1:2)] / 413) / mean(eco.rates$time_rates[, 2][-(1:2)] / 40)


# Statistical tests:
wilcox.test(eco.rates$time_rates[, 2] / 40, morph.rates$time_rates[, 2] / 413)
# W = 138, p < 2.2e-16: morphology greater


plot(diff(eco.rates$time_rates[, 2] / 40),
     diff(morph.rates$time_rates[, 2] / 413))
cor(diff(eco.rates$time_rates[, 2] / 40),
    diff(morph.rates$time_rates[, 2] / 413)) # r = + 0.423
lm.diff <- lm(diff(eco.rates$time_rates[, 2] / 40) ~ diff(morph.rates$time_rates[, 2] / 413))
summary(lm.diff) # p = 0.0028, strongly correlated

# Remove first two time intervals because outliers
plot(diff(eco.rates$time_rates[, 2] / 40)[-c(1:2)],
     diff(morph.rates$time_rates[, 2] / 413)[-c(1:2)])
cor(diff(eco.rates$time_rates[, 2] / 40)[-c(1:2)],
    diff(morph.rates$time_rates[, 2] / 413)[-c(1:2)]) # r = + 0.204
lm.diff <- lm(diff(eco.rates$time_rates[, 2] / 40)[-c(1:2)] ~ diff(morph.rates$time_rates[, 2] / 413)[-c(1:2)])
summary(lm.diff) # p = 0.173, so only correlated because of the initially high rates
