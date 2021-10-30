## DISPARITY ANALYSES USING CLADDIS ############################################

# Prior to running, run 2-InferAncestralStates.R to infer ancestral states using
# 'Claddis' package.

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
library(ade4)       # v. 1.7-17
library(Claddis)    # v. 0.6.3 - Check SI for Lloyd 2018 for walk-through on code functions
if(packageVersion("Claddis") < "0.4.1")
  stop("wrong version of 'Claddis' Get updated version from GitHub or CRAN.\n")
library(snowfall)   # v. 1.84-6.1 
library(parallel)   # v. 4.1.0



### FUNCTIONS ##################################################################

# Prepare snowfall cluster. A simpler way to run the initialization code. Make
# sure to Stop with sfStop().
#    CPUs = no. CPUs to use. Default is all available.
#    oufile = file to store parallel output status
prepare_SF <- function(CPUs = "max", outfile = "initfile") {
  require(snowfall)
  require(parallel)
  if(CPUs == "max")
    cpus <- parallel::detectCores() else cpus <- CPUs
  sfInit(parallel = TRUE, cpus = cpus, slaveOutfile = outfile)
  stopifnot(sfCpus() == cpus)		    # Confirm set up CPUs properly
  stopifnot(sfParallel() == TRUE)		# Confirm now running in parallel
  sfExportAll()				              # Export all libraries, files, & objects
}

# For the Hopkins and St. Johns distances, using a modified version of
# Claddis::calculate_morphological_distances that silences line 416 (in
# annotated GitHub version) to avoid triggering CharDepends error (which in the
# morphological case is not really an error because the dependencies are not
# always straightforward when coding phylum-wide characters. Thanks, Graeme!
source("~/Manuscripts/CamOrdEchinos/calculate_morphological_distances2.R")


## IMPORT FILES ################################################################

# See 0-DataSummaries.R for description of each data treatment. See
# 2-InferAncestralStates.R for how created.

# load("mode.anc"); anc <- mode.anc; rm("mode.anc")
# load("constant.anc"); anc <- constant.anc; rm("constant.anc")
# load("raw.anc"); anc <- raw.anc; rm("raw.anc")
load("morph.anc"); anc <- morph.anc; rm("morph.anc")
length(anc)
anc[[1]]$topper
head(anc[[1]]$matrix_1$matrix)

## Processing of the distance matrices for the morphological data set tends to
## hit memory allocation errors for the GED distance metrics. Breaking into
## tenths to improve memory efficiency.
# anc1 <- anc[1:5] # Used for all GEDs: 1:7 and greater are too large
# anc2 <- anc[6:10]
# anc3 <- anc[11:15]
# anc4 <- anc[16:20]
# anc5 <- anc[21:25]
# anc6 <- anc[26:30]
# anc7 <- anc[31:35]
# anc8 <- anc[36:40]
# anc9 <- anc[41:45]
# anc10 <- anc[46:50]
# save(anc1, file = "anc1"); save(anc2, file = "anc2"); save(anc3, file = "anc3")
# save(anc4, file = "anc4"); save(anc5, file = "anc5"); save(anc6, file = "anc6")
# save(anc7, file = "anc7"); save(anc8, file = "anc8"); save(anc9, file = "anc9")
# save(anc10, file = "anc10")
# rm(list = c("anc", "anc2", "anc3", "anc4", "anc5", "anc6", "anc7", "anc8", "anc9", "anc10"))
# gc()
# load("anc2") # Load as needed to increase functional memory


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

# Unless specified, behavior for polymorphisms (= "min_difference"), uncertain
# (= "min_difference"), and inapplicable (= "missing") characters use the
# default settings for the function.

# Saving each output in case of crash or outage.

# MORD: Maximum observable rescaled distance
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.MORD <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances, 
                   distance_metric = "mord", 
                   distance_transformation =  "arcsine_sqrt")
sfStop()
(Sys.time() - t.start0) # ~ 12-14 min for life habit, ~ 62 min for morph
save(distances.MORD, file = "distances.MORD")
beep(3)

# RED: Raw Euclidean distance
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.RED <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances, 
                   distance_metric = "red")
sfStop()
(Sys.time() - t.start0) # ~ 14-15 min for life habit, ~ 71 min for morph
save(distances.RED, file = "distances.RED")
beep(3)

# GC: Gower coefficient
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.GC <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances, 
                   distance_metric = "gc", 
                   distance_transformation =  "arcsine_sqrt")
sfStop()
(Sys.time() - t.start0) # ~ 14-17 min for life habit, ~ 64 min for morph
save(distances.GC, file = "distances.GC")
beep(3)

# GED: Wills' Generalized Euclidean distance (ignoring dependencies)
# - If using the morphological data set, change 'anc' to 'anc1' or 'anc2', and
#    do same for output distance matrices. Because of memory limits, can only 
#    handle list of up to 5 trees, so using 5 CPUs, 10 batches, and no 
#    load-balancing.
(t.start0 <- Sys.time())
prepare_SF()
# prepare_SF(CPUs = length(anc1))
sfLibrary(Claddis)
distances.GED <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances, 
                   distance_metric = "ged", ged_type = "wills")
# Following only for morphological data set, updating 'anc1' and '....p1' as relevant
# distances.GED.p2 <- sfLapply(x = anc2, fun = calculate_morphological_distances, distance_metric = "ged", ged_type = "wills")
sfStop()
(Sys.time() - t.start0) # ~ 15-24 min for life habit, ~ 16 * 10 = 160 min for morph
save(distances.GED, file = "distances.GED")
# save(distances.GED.p2, file = "distances.GED.p2")
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 0 (H&SJ 2018)
# - "When Alpha = 0 the secondary [and tertiary, etc.] characters contribute
# nothing to the [distance] comparisons of the two taxa." (Hopkins and St. John
# 2018: p. 3)
# - Note that this is run with the "silenced" character dependencies warning.
# - If using the morphological data set, change 'anc' to 'anc1' or 'anc2', and
#    do same for output distance matrices. Because of memory limits, can only 
#    handle list of up to 5 trees, so using 5 CPUs, 10 batches, and no 
#    load-balancing.
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.GED0 <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances2, 
                   distance_metric = "ged", ged_type = "wills",
                   alpha = 0, inapplicable_behaviour = "hsj",
                   character_dependencies = CharDepends)
# Following only for morphological data set, updating 'anc1' and '....p1' as relevant
# distances.GED0.p1 <- sfClusterApplyLB(x = anc1, fun = calculate_morphological_distances2, distance_metric = "ged", ged_type = "wills", alpha = 0, inapplicable_behaviour = "hsj", character_dependencies = CharDepends)
sfStop()
(Sys.time() - t.start0) # ~  20-27 min for life habit, ~ 44 * 10 = 440 min for morph
save(distances.GED0, file = "distances.GED0")
# save(distances.GED0.p1, file = "distances.GED0.p1")
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 0.5 (H&SJ 2018)
# - "When Alpha = 0.5, the primary character contributes weight according to the
#    fraction of shared secondary [and tertiary, etc.] characters." (Hopkins and
#    St. John 2018: p. 3)
# - Note that this is run with the "silenced" character dependencies warning.
# - If using the morphological data set, change 'anc' to 'anc1' or 'anc2', and
#    do same for output distance matrices. Because of memory limits, can only 
#    handle list of up to 5 trees, so using 5 CPUs, 10 batches, and no 
#    load-balancing.
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.GED.5 <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances2, 
                   distance_metric = "ged", ged_type = "wills",
                   alpha = 0.5, inapplicable_behaviour = "hsj",
                   character_dependencies = CharDepends)
# Following only for morphological data set, updating 'anc1' and '....p1' as relevant
# distances.GED.5.p1 <- sfClusterApplyLB(x = anc1, fun = calculate_morphological_distances2, distance_metric = "ged", ged_type = "wills", alpha = 0.5, inapplicable_behaviour = "hsj", character_dependencies = CharDepends)
sfStop()
(Sys.time() - t.start0) # ~ 24-29 min for life habit, ~ 45 min * 10 = 450 min for morph
save(distances.GED.5, file = "distances.GED.5")
# save(distances.GED.5.p1, file = "distances.GED.5.p1")
beep(3)

# GED: Wills' Generalized Euclidean distance, with alpha = 1 (H&SJ 2018)
# - When Alpha = 1, the primary character is not counted in the weight
#    separately." (Hopkins and St. John 2018: p. 3)
# - Note that this is run with the "silenced" character dependencies warning.
# - If using the morphological data set, change 'anc' to 'anc1' or 'anc2', and
#    do same for output distance matrices. Because of memory limits, can only 
#    handle list of up to 5 trees, so using 5 CPUs, 10 batches, and no 
#    load-balancing.
(t.start0 <- Sys.time())
prepare_SF()
sfLibrary(Claddis)
distances.GED1 <-
  sfClusterApplyLB(x = anc, fun = calculate_morphological_distances2, 
                   distance_metric = "ged", ged_type = "wills",
                   alpha = 1, inapplicable_behaviour = "hsj",
                   character_dependencies = CharDepends)
# Following only for morphological data set, updating 'anc1' and '....p1' as relevant
# distances.GED1.p1 <- sfClusterApplyLB(x = anc1, fun = calculate_morphological_distances2, distance_metric = "ged", ged_type = "wills", alpha = 1, inapplicable_behaviour = "hsj", character_dependencies = CharDepends)
sfStop()
(Sys.time() - t.start0) # ~ 26-33 min for life habit, ~ 47 min * 10 = 470 min for morph
save(distances.GED1, file = "distances.GED1")
# save(distances.GED1.p1, file = "distances.GED1.p1")
beep(3)


# Load and combine (ONLY use for morphological data set)
# load("distances.GED.p1"); load("distances.GED0.p1"); load("distances.GED.5.p1")
# load("distances.GED1.p1"); load("distances.GED.p2"); load("distances.GED0.p2")
# load("distances.GED.5.p2"); load("distances.GED1.p2"); load("distances.GED.p3")
# load("distances.GED0.p3"); load("distances.GED.5.p3"); load("distances.GED1.p3")
# load("distances.GED.p4"); load("distances.GED0.p4"); load("distances.GED.5.p4")
# load("distances.GED1.p4"); load("distances.GED.p5"); load("distances.GED0.p5")
# load("distances.GED.5.p5"); load("distances.GED1.p5"); load("distances.GED.p6")
# load("distances.GED0.p6"); load("distances.GED.5.p6"); load("distances.GED1.p6")
# load("distances.GED.p7"); load("distances.GED0.p7"); load("distances.GED.5.p7")
# load("distances.GED1.p7"); load("distances.GED.p8"); load("distances.GED0.p8")
# load("distances.GED.5.p8"); load("distances.GED1.p8"); load("distances.GED.p9")
# load("distances.GED0.p9"); load("distances.GED.5.p9"); load("distances.GED1.p9")
# load("distances.GED.p10"); load("distances.GED0.p10"); load("distances.GED.5.p10")
# load("distances.GED1.p10")

# Combine (ONLY use for morphological data set)
# distances.GED <- c(distances.GED.p1, distances.GED.p2, distances.GED.p3, distances.GED.p4, distances.GED.p5, distances.GED.p6, distances.GED.p7, distances.GED.p8, distances.GED.p9, distances.GED.p10)
# distances.GED0 <- c(distances.GED0.p1, distances.GED0.p2, distances.GED0.p3, distances.GED0.p4, distances.GED0.p5, distances.GED0.p6, distances.GED0.p7, distances.GED0.p8, distances.GED0.p9, distances.GED0.p10)
# distances.GED.5 <- c(distances.GED.5.p1, distances.GED.5.p2, distances.GED.5.p3, distances.GED.5.p4, distances.GED.5.p5, distances.GED.5.p6, distances.GED.5.p7, distances.GED.5.p8, distances.GED.5.p9, distances.GED.5.p10)
# distances.GED1 <- c(distances.GED1.p1, distances.GED1.p2, distances.GED1.p3, distances.GED1.p4, distances.GED1.p5, distances.GED1.p6, distances.GED1.p7, distances.GED1.p8, distances.GED1.p9, distances.GED1.p10)
# save(distances.GED, file = "distances.GED")
# save(distances.GED0, file = "distances.GED0")
# save(distances.GED.5, file = "distances.GED.5")
# save(distances.GED1, file = "distances.GED1")


# Save distance matrices (appending data treatment as prefix to name)
matrices <- c("distances.RED", "distances.MORD", "distances.GC", "distances.GED",
              "distances.GED0", "distances.GED.5", "distances.GED1")
# Change this depending on data treatment: mode, constant, raw, or morph
treatment <- "morph"
(matrix.names <- paste0(treatment, ".", matrices))
# for(f in seq(matrices)){
#   assign(matrix.names[f], get(matrices[f]))
#   save(list = matrix.names[f], file = matrix.names[f])
# }
# beep(3)

# Reload using same logic
for(f in seq(matrices)){
  tmp.name <- load(matrix.names[f])
  tmp.obj <- get(tmp.name)
  names(tmp.obj) <- matrix.names[f]
}
beep()


## COMPARE DISTANCE METRICS USING MANTEL TEST ##################################

# Create objects to store results
dist.dists <- matrix(NA, 7, 7)
colnames(dist.dists) <- c("RED", "MORD", "GC", "GED", "GED0", "GED.5", "GED1") 
rownames(dist.dists) <- c("RED", "MORD", "GC", "GED", "GED0", "GED.5", "GED1")
mean.dists <- dist.dists
sq <- 1:50
list.dists <- lapply(sq, function(sq) dist.dists)
index <- c(1, seq(from = 5, to = 50, by = 5))

# Compare with Mantel test. (Much more complicated than needed, to handle the
# 'raw' treatment with many NAs. Only set 'TRUE' for the 'raw' treatment.)
remove.NAs <- FALSE
remove.zeros <- FALSE
set.seed(3124) # Set because Mantel test is stochastic
(t.start0 <- Sys.time())
for(l in seq(length(list.dists))) {
# for(l in 36:length(list.dists)) {
    if (l %in% index)
    cat("processing tree", l, "\n")
  for (r in 1:6) {
    for (c in (r+1):7) {
      d1 <- get(matrix.names[r])[[l]]$distance_matrix
      d2 <- get(matrix.names[c])[[l]]$distance_matrix
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
        d1 <- replace(d1, d1 == 0, min(d1[d1 > 0], na.rm = TRUE) / 2)
        d2 <- replace(d2, d2 == 0, min(d2[d2 > 0], na.rm = TRUE) / 2)
        if (length(d1) != length(d2))
          stop("inconsisent lengths when NAs removed.\n")
      }
      list.dists[[l]][r, c] <- ade4::mantel.rtest(d1, d2, nrepet = 99)$obs
    }
  }
}
(Sys.time() - t.start0) # Takes ~ 1 hour to run. Not worth re-writing code to run in parallel.
beep(3)

# Distance matrices 31-35 are problematic for morphology (b/c they have many
# more non-diagonal 0s and NAs than the other ancestral-tree-reconstruction
# distance matrices). The issue applies to all distance metrics. Investigation
# of the previously processed ancestral-state reconstructions and time-scaled
# trees show they are not improper. Re-running samples (e.g., re-estimating
# ancestral states) produces the same outcomes as before. The problem appears to
# be caused by idiosyncratic ancestral reconstructions for these trees, with
# polymorphisms and unknown states contributing to 0 'minimum' distances between
# ancestral nodes (see Claddis:calculate_morphological_distances() for details
# of how distances are calculated in the case of polymorphisms, uncertain, or
# inapplicable character states). Given the diagnostic testing, it seems the
# fact that the 5 idiosyncratic trees are consecutive is simply due to
# coincidence. (Note that the same time trees were used in the 3 life habit
# treatments, without any problems.) Therefore, proceeding as normal, but
# ignoring these 5 distance matrices in the Mantel test.

# Convert list to means and summarize
for (r in 1:6) {
  for (c in (r+1):7) {
    mean.dists[r,c] <- mean(sapply(sq, function(sq) list.dists[[sq]][r,c]), 
                            na.rm = TRUE)
  }
}
round(mean.dists, 3)
min(mean.dists, na.rm = TRUE)

# Across-metric averages (divide by 6 so ignores necessary 0 in row/col sums):
sort(((rowSums(mean.dists, na.rm = TRUE) + colSums(mean.dists, na.rm = TRUE)))/ 6)

# RESULTS (only removing NAs/zeros for "raw" treatment"):
# Mode:     All strongly correlated (min R = 0.966; w/o RED = 0.966); RED == GED. GEDs and RED 
#                    most highly correlated with each other; GC and MORD 
#                    slightly less correlated with others.
# Constant: All strongly correlated (min R = 0.954; w/o RED = 0.954); GEDs and RED most highly 
#                    correlated with each other; GC and MORD slightly less 
#                    correlated with others.
# Raw:      All strongly correlated (min R = 0.931; w/o RED = 0.934), but had to remove many NAs 
#                    and zeros. GEDs and RED most highly correlated with each 
#                    other; GC and MORD slightly less correlated with others.
# Morph:    All moderately strongly correlated (min R = 0.645; w/o RED = 0.819), with RED lowest
#                    among all pairs, and GED most strongly correlated with the
#                    others. GED, GED1, GC, and GED.5 are most correlated, with 
#                    MORD, GED0, and RED lowest, on average.

# CONCLUSIONS: Given the many hierarchically dependent characters and the fact
# that Wills' GED is very highly correlated across distance measures, we will
# use GED.5 as our standard metric. Using GED0 means ignoring 89% of the
# morphological characters and using GED1 means ignoring the 11% primary
# characters that are generally coded for all taxa. (Relatedly, it is sensible
# that GED is consistently very highly correlated with GED1 because in the
# morphological data set they use 89% identical data and in the ecological data
# set they use 98% identical data.)