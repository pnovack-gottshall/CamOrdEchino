## DATA SUMMARIES ##############################################################

# There are four data treatments of two data sets used in the analyses, plus an
# informal supertree of Cambrian - Ordovician echinoderms. The supertree has 366
# genera and was assembled by Colin Sumrall, Brad Deline, and associates. The
# morphological data set was compiled by Brad Deline and associates. The
# ecological (= functional or life-habit) data set was compiled by Phil
# Novack-Gottshall and associates (including co-authors Ali Sultan, Jack
# Purcell, and Isa Ranjha).

# The ecological data set has 40 characters: 5 numerics, 1 ordered factors, and
# 34 binaries. Because the ecology of many early echinoderm genera have not been
# formally studied, missing ecological character states were inferred using
# three data treatments spanning a range of conservative to liberal assumptions.

#  * 'Mode' infers missing life-habit states using the most frequent state of
#      (taxonomically) closely related taxa. This data set has 40 characters: 5
#      numerics, 1 ordered factors, and 34 binaries. One of the characters is
#      dependent on another character. 98.6% of possible tip states are coded,
#      with 1.2% missing (unknown or unable to be inferred by the 'mode' 
#      algorithm) and 0.3% inapplicable.

#  * 'Constant' infers missing life-habit states only when all (taxonomically)
#      closely related taxa share the same state. 91.6% of possible tip states 
#      are coded, with 8.1% missing (unknown or unable to be inferred by the 
#      'constant' algorithm) and 0.2% inapplicable.

#  * 'Raw' only uses life-habit states that are known, making no assumptions
#      about missing data. 52.0% of possible tip states are coded, with 47.9%
#      missing (unknown) and 0.1% inapplicable. For analyses that trim out taxa with
#      uncoded states, this data treatment includes only 256 genera.

#  * 'Morph' includes morphological data (independent of cladistic characters).
#      The data set has 413 characters: 30 numerics, 68 ordered factors, 147
#      unordered factors, and 168 binaries. 369 of the characters are dependent 
#      on other characters, with levels spanning secondary through quinary 
#      dependents. 28.6% of possible tip states are coded, with 9.1% missing 
#      (unknown) and 62.3% inapplicable.


## PREPARATIONS ################################################################
rm(list = ls())
op <- par()

# Set working directory (point to the folder containing the input files on your
# own machine):
# setwd("[filepath to folder containing data files on your personal machine]")
setwd("~/Manuscripts/CamOrdEchinos/Data files/NA Reformatted")



## IMPORT RAW DATA #############################################################
# Note these files are the source data from which all other data files used here
# were prepared. These are archival and used here only to describe the data.

# Morphological data set:
morph <- read.csv(file = "DelineMorphData_NAreformatted.csv",
                  header = TRUE, stringsAsFactors = FALSE)
# Unknown states coded as ? and inapplicable states coded as -
head(morph)

# Ecological data set / treatments:
eco <- read.csv(file = "EchinoLHData_Mode_NAreformatted.csv", 
                header = TRUE, stringsAsFactors = FALSE)
# Order the taxonomic scales:
scales <- c("Species", "Subgenus", "Genus", "Subfamily", "Family", "Superfamily", 
            "Suborder", "Order", "Subclass", "Class", "Subphylum", "Phylum", "", NA)
scales <- factor(scales, levels = scales, ordered = TRUE)
eco$EcologyScale <- factor(eco$EcologyScale, levels = scales, ordered = TRUE)
eco$BodySizeScale <- factor(eco$BodySizeScale, levels = scales, ordered = TRUE)
head(eco)



## SUMMARIZE DATA SETS #########################################################

## How resolved are the life habits?
table(eco$EcologyScale)
round(100 * table(eco$EcologyScale) / nrow(eco), 1)
round(cumsum(100 * table(eco$EcologyScale) / nrow(eco)), 1)
# 55% genera have life habits coded at species level, 62% at genus or better,
# 86% at family or better, 95% at order or better

table(eco$BodySizeScale)
round(100 * table(eco$BodySizeScale) / nrow(eco), 1)
round(cumsum(100 * table(eco$BodySizeScale) / nrow(eco)), 1)
# 65% genera have body sizes measured at species level, 71% at genus or better,
# 91% at family or better, 96% at order or better


## Taxonomic coverage
length(unique(eco$Genus)) # 365 genera
# (2 subgenera Anatiferocystis and Guichenocarpos in genus Anatifopsis treated
# as different entries in all data sets.)

# 55 orders, none with more than 33 genera
length(unique(eco$Order))
sort(table(eco$Order))

# 25 classes, only crinoids (136 genera) with more than 36 genera
length(unique(eco$Class))
sort(table(eco$Class))


## How many unique morphotypes and life habits?
# Using the Wills GED distance matrices as they better accommodate differences
# due to missing or inapplicable data. Using time-scaled tree #50 as an example.

# See 3-DisparityDistances.R for code used to build these objects.
load("~/Manuscripts/CamOrdEchinos/Data files/NA reformatted/mode.distances.GED.5")
load("~/Manuscripts/CamOrdEchinos/Data files/NA reformatted/constant.distances.GED.5")
load("~/Manuscripts/CamOrdEchinos/Data files/NA reformatted/raw.distances.GED.5")
load("~/Manuscripts/CamOrdEchinos/Data files/NA reformatted/morph.distances.GED.5")

nrow(unique(mode.distances.GED.5[[50]]$distance_matrix))
# 203 unique life habits in mode treatment

nrow(unique(constant.distances.GED.5[[50]]$distance_matrix))
# 309 unique life habits in constant treatment

nrow(unique(raw.distances.GED.5[[50]]$distance_matrix))
# 395 unique life habits in raw treatment

nrow(unique(morph.distances.GED.5[[50]]$distance_matrix))
# 731 unique morphotypes (all unique!)
