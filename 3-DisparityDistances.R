## DISPARITY ANALYSES USING CLADDIS ############################################

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' package.

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
setwd("~/Manuscripts/CamOrdEchino/Data files/NA Reformatted")

# Load packages
library(beepr)      # v. 1.3
library(ade4)       # v. 1.7-15
library(Claddis)    # v. 0.4.1 - Check SI for Lloyd 2018 for walk-through on code functions (in R experiments folder)
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub.\n")


## IMPORT FILES ################################################################

# See 0-DataSummaries.R for description of each data treatment.

# load("mode.anc"); anc <- mode.anc; rm("mode.anc")
# load("constant.anc"); anc <- constant.anc; rm("constant.anc")
# load("raw.anc"); anc <- raw.anc; rm("raw.anc")
# load("morph.anc"); anc <- morph.anc; rm("morph.anc")
anc$Topper
head(anc$Matrix_1$Matrix)

## Import character dependencies
# Two-column MATRIX, where first column is named "DependentCharacter" (e.g., the
# secondary, tertiary, etc. character) and second column is named
# "IndependentCharacter" (e.g., the primary character it depends on). Each row
# is one linkage. Make sure the character values (i.e., CHARLABELS in the Nexus
# files) here match those in the morphological and ecological data matrices, and
# that character values do not restart in subsequent data blocks! (But Graeme's
# user-friendly code does not enforce this requirement, so long as the column
# names are correct.)
CharDepends <- as.matrix(read.csv("CharDepends_Eco.csv", header = TRUE))
# CharDepends <- as.matrix(read.csv("CharDepends_Morph.csv", header = TRUE))
head(CharDepends)




## CALCULATE DISPARITY DISTANCE MATRICES #######################################

# Lloyd (2016) does sensitivity tests of four distance metrics (raw Euclidean
# distance [RED], general Euclidean distance [GED], Gower's coefficient [GC],
# and maximum observable re-scaled distance [MORD, equal to Gower when data are
# all binary or unordered]. All are fine, but best depends on proportion of
# missing data, presence of ordered data, and size of data set. Trying all here,
# a philosophy analogous to that used in non-phylogenetic disparity analyses.

# For the Hopkins and St. Johns distances, using a modified version of
# Claddis::MorphDistMatrix that silences line 396 to avoid triggering
# CharDepends error (which in the morphological case is not really an error
# because the dependencies are not always straightforward when coding
# phylum-wide characters. Thanks, Graeme!
source("~/Manuscripts/CamOrdEchino/MorphDistMatrix2.R")

# MORD: Maximum observable rescaled distance
(t.start0 <- Sys.time())
distances.MORD <- MorphDistMatrix(anc, Distance = "MORD", 
                                  TransformDistances = "arcsine_sqrt")
(Sys.time() - t.start0) # ~ 1 min for life habit, ~ 6 min for morph
beep(3)

# RED: Raw Euclidean distance
(t.start0 <- Sys.time())
distances.RED <- MorphDistMatrix(anc, Distance = "RED")
(Sys.time() - t.start0) # ~ 1 min for life habit, ~ 4 min for morph
beep(3)

# GC: Gower coefficient
(t.start0 <- Sys.time())
distances.GC <- MorphDistMatrix(anc, Distance = "GC", 
                                TransformDistances = "arcsine_sqrt")
(Sys.time() - t.start0) # ~ 1 min for life habit, ~ 4 min for morph
beep(3)

# GED: Wills' Generalized Euclidean distance (ignoring dependencies)
(t.start0 <- Sys.time())
distances.GED <- MorphDistMatrix(anc, Distance = "GED", GEDType = "Wills")
(Sys.time() - t.start0) # ~ 2 min for life habit, ~ 5 min for morph
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 0 (H&SJ 2018)
# "When Alpha = 0 the secondary [and tertiary, etc.] characters contribute
# nothing to the [distance] comparisons of the two taxa." (Hopkins and St. John
# 2018: p. 3)
# Note that this is run with the "silenced" character dependencies warning.
(t.start0 <- Sys.time())
distances.GED0 <- MorphDistMatrix2(anc, Distance = "GED", GEDType = "Wills", 
                                   Alpha = 0, InapplicableBehaviour = "HSJ",
                                   CharacterDependencies = CharDepends)
(Sys.time() - t.start0) # ~ 2 min for life habit, ~ 25 min for morph
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 0.5 (H&SJ 2018)
# "When Alpha = 0.5, the primary character contributes weight according to the
# fraction of shared secondary [and tertiary, etc.] characters." (Hopkins and
# St. John 2018: p. 3)
# Note that this is run with the "silenced" character dependencies warning.
(t.start0 <- Sys.time())
distances.GED.5 <- MorphDistMatrix2(anc, Distance = "GED", GEDType = "Wills", 
                                    Alpha = 0.5, InapplicableBehaviour = "HSJ",
                                    CharacterDependencies = CharDepends)
(Sys.time() - t.start0) # ~ 2 min for life habit, ~ 24 min for morph
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 1 (H&SJ 2018)
# When Alpha = 1, the primary character is not counted in the weight
# separately." (Hopkins and St. John 2018: p. 3)
# Note that this is run with the "silenced" character dependencies warning.
(t.start0 <- Sys.time())
distances.GED1 <- MorphDistMatrix2(anc, Distance = "GED", GEDType = "Wills", 
                                   Alpha = 1, InapplicableBehaviour = "HSJ",
                                   CharacterDependencies = CharDepends)
(Sys.time() - t.start0) # ~ 2 min for life habit, ~ 44 min for morph
beep(3)


# Save distance matrices (appending data treatment as prefix to name)
matrices <- c("distances.RED", "distances.MORD", "distances.GC", "distances.GED",
              "distances.GED0", "distances.GED.5", "distances.GED1")
treatment <- "morph"                    # Change this depending on data treatment
matrix.names <- paste0(treatment, ".", matrices)
# for(f in seq(matrices)){
#   assign(matrix.names[f], get(matrices[f]))
#   save(list = matrix.names[f], file = matrix.names[f])
# }


# Reload using same logic
for(f in seq(matrices)){
  tmp.name <- load(matrix.names[f])
  tmp.obj <- get(tmp.name)
  names(tmp.obj) <- matrix.names[f]
}

# For GED.5, load manually (because subsequently resaved differently)
load("mode.distances.GED.5")
load("constant.distances.GED.5")
load("raw.distances.GED.5")
load("morph.distances.GED.5")


# Compare distance matrices using Mantel test
dist.dists <- matrix(NA, 7, 7)
colnames(dist.dists) <- c("RED", "MORD", "GC", "GED", "GED0", "GED.5", "GED1") 
rownames(dist.dists) <- c("RED", "MORD", "GC", "GED", "GED0", "GED.5", "GED1")
# Much more complicated than needed to handle the 'raw' treatment with many NAs
remove.NAs <- FALSE
remove.zeros <- FALSE
set.seed(3124)
for (r in 1:6) {
  for (c in (r+1):7) {
    d1 <- get(matrix.names[r])$DistanceMatrix
    d2 <- get(matrix.names[c])$DistanceMatrix
    if (remove.NAs) {
      diag(d1) <- diag(d2) <- NA
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
  }
}
beep(3)
round(dist.dists, 3)
min(dist.dists, na.rm = TRUE)
# Across-metric averages:
sort(((rowSums(dist.dists, na.rm = TRUE) + colSums(dist.dists, na.rm = TRUE)))
      / 6)
# RESULTS (only removing NAs/zeros for "raw" treatment"):
# Mode:     All strongly correlated (min R = 0.966); RED == GED == GED1. GED,
#                     GED1, RED, and GED.5 are most correlated.
# Constant: All strongly correlated (min R = 0.954); GED == GED1. GED, GED1,
#                     GED.5, and RED are most correlated.
# Raw:      All strongly correlated (min R = 0.932), but had to remove many NAs 
#                    and zeros. GED = GED1. GED, GED1, GED.5, and RED are most
#                    correlated.
# Morph:    All moderately strongly correlated (min R = 0.733), with RED lowest
#                    among all pairs, and GED most strongly correlated with the
#                    others. GED, GED1, GED.5, and GC are most correlated.

# CONCLUSIONS: Given the many hierarchically dependent characters and the fact
# that Wills' GED is very highly correlated across distance measures, we will
# use GED.5 as our standard metric. Using GED0 means ignoring 89% of the
# morphological characters and using GED1 means ignoring the 11% primary
# characters that are generally coded for all taxa. (Relatedly, it is sensible
# that GED is consistently very highly correlated with GED1 because in the
# morphological data set they use 89% identical data and in the ecological data
# set they use 98% identical data.)


