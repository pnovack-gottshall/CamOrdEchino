## FUNCTIONAL DIVERSITY FUNCTIONS ##############################################
# Functions from FD::dbFD() to directly calculate functional diversity metrics
# FRic (plus qual.FRic), FEve, FDis, and FDiv given a pre-calculated PCoA
# eigenvector matrix and associated distance matrix. Also allows usage of the
# Dineen, et al. (2019) standardization of PCoA eigenvectors by eigenvalue
# magnitude. Modified from FD::dbFD written by Etienne Laliberté. 
# Thanks, Etienne!
# 
#               pcoa = output from ape::pcoa, with $correction noting the  
#                      eigenvalue correction used (e.g., Cailliez, etc.), and, 
#                      if so, using vectors.cor for corrected eigenvectors.
#        dist.matrix = distance matrix used to create pcoa.
#                  m = number of PCoA eigenvalues to reduce dimensionality to.
#         stand.pcoa = logical indicating whether to standardize PCoA  
#                      eigenvectors according to magnitude of eigenvalues. 
#                      Default = TRUE.
# calc.FRic.and.FDiv = whether to calculate FRic and FDiv. Default = TRUE.

simple.dbFD <- function(pcoa = NULL, dist.matrix = NULL, m = NA, 
                        stand.pcoa = TRUE, calc.FRic.and.FDiv = TRUE) {
  
  # Functional to apply Dineen, et al. (2019, Biology Lettters) standardization
  # of PCoA eigenvectors according to the magnitude of the first eigenvalue.
  # (Their eq. 1 in Supplementary Information.)
  stand.pcoa.fn <- function(vectors = NULL, eigenvalues = NULL) {
    nc <- ncol(vectors)
    first.e <- eigenvalues[1]
    for (c in 1:nc) {
      vectors[, c] <- vectors[, c] * eigenvalues[c] / first.e
    }
    return(vectors)
  }
  
  # Set default values (if incalculable) as NA
  FRic <- FEve <- FDis <- FDiv <- qual.FRic <- NA
  
  if (!"dist" %in% class(dist.matrix))
    dist.matrix <- as.dist(dist.matrix)
  
  # Calculate FDis first, because does not use PCoA vectors
  
  # But first confirm the distance matrix is complete
  if (all(!is.na(dist.matrix))) {
    
    ## FDis: Functional dispersion
    # No modification needed. Use dbFD::fdisp(d).
    #    d = distance matrix (in class dist)
    FDis <- unname(FD::fdisp(d = dist.matrix)$FDis)
    
  }
  
  # Confirm pcoa exists before calculating rest
  if (!is.null(pcoa)) {
    
    # Get components (with switch for whether corrected or not)
    if (pcoa$correction[1] == "none") {
      vectors <- pcoa$vectors
      eigenvalues <- pcoa$values$Eigenvalues
    } else {
      vectors <- pcoa$vectors.cor
      eigenvalues <- pcoa$values$Corr_eig
    }
      
    # Standardize the eigenvectors according to magnitude of eigenvalues (if desired)
    if(stand.pcoa)
      vectors <- stand.pcoa.fn(vectors = vectors, eigenvalues = eigenvalues)
    if (is.na(m))
      stop("Must supply 'm' to calculate functional diversity metrics.")
    
    # The eigenvectors with 'm' reduced dimensionality:
    m.vectors <- vectors[, 1:m]
    
    ## FRic: Functional richness
    # With earlier processing of reduced dimensionality and standardized PCoA, no
    # novel function needed.
    FRic <- FDiv <- FEve <- NA
    
    if (calc.FRic.and.FDiv) {
      FRic <- geometry::convhulln(m.vectors, "FA")$vol
    
      ## FDiv: Functional divergence
      #    v = eigenvector matrix output from ape::pcoa (with 'm' reduced 
      #        dimensionality and stand.pcoa pre-standardized).
      calc.FDiv <- function(v = NULL) {
        # calculate vertices
        vert0 <- geometry::convhulln(v, "Fx TO 'vert.txt'")
        vert1 <- scan("vert.txt", quiet = TRUE)
        vert2 <- vert1 + 1
        vertices <- vert2[-1]
        # traits values of vertices
        trvertices <- v[vertices, ]
        # coordinates of the center of gravity of the vertices (Gv)
        baryv <- apply(trvertices, 2, mean)
        # Euclidean distances to Gv (dB)
        S <- nrow(v)
        distbaryv <- rep(0, S)
        for (j in 1:S) {
          distbaryv[j] <- (sum((v[j, ] - baryv) ^ 2)) ^ 0.5
        }
        # mean of dB values
        meandB <- mean(distbaryv)
        # deviations to mean of db
        devdB <- distbaryv - meandB
        # computation of FDiv (using equal-weighting, i.e., ignoring rel. abund.)
        out <- (sum(devdB) + meandB) / (sum(abs(devdB)) + meandB)
        return(out)
      }
    
      FDiv <- calc.FDiv(v = m.vectors)
    }
    
    ## FEve: Functional evenness
    #    v = eigenvector matrix output from ape::pcoa (with all 'm' eigenvalues
    #        and with stand.pcoa pre-standardized).
    calc.FEve <- function(v = NULL) {
      # computation of minimum spanning tree and conversion of the 'mst' matrix
      # into 'dist' class
      tr.dist <- dist(v)
      # The next step is very slow (for large data sets)!!!
      linkmst <- ape::mst(tr.dist)
      mstvect <- as.dist(linkmst)
      # computation of the pairwise cumulative relative abundances and conversion
      # into 'dist' class
      S <- nrow(v)
      # Use equal-weighting (ignoring relative abundance)
      abund2 <- matrix(1, nrow = S, ncol = S)
      abund2vect <- as.dist(abund2)
      # computation of EW for the (S - 1) segments to link the S points
      EW <- rep(0, S - 1)
      flag <- 1
      for (m in 1:((S - 1) * S / 2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m] / (abund2vect[m])
          flag <- flag + 1
        }
      }
      # computation of the PEW and comparison with 1 / S - 1, finally computation
      # of FEve
      minPEW <- rep(0, S - 1)
      OdSmO <- 1 / (S - 1)
      for (l in 1:(S - 1)) {
        minPEW[l] <- min((EW[l] / sum(EW)), OdSmO)
      }
      FEve <- ((sum(minPEW)) - OdSmO) / (1 - OdSmO)
      return(FEve)
    }
    
    FEve <- calc.FEve(v = vectors)
    
    ## qual.FRic: Quality of reduced-dimensionality hyperspace used in FRic
    #    v = eigenvector matrix output from ape::pcoa (with 'm' reduced 
    #        dimensionality and stand.pcoa pre-standardized).
    calc.qual.FRic <- function(pcoa = NULL, dist.matrix = NULL, m = NA) {
  
      # What 'correction' was used in ape::pcoa to handle negative eigenvalues?
      corr <- pcoa$correction[1]
      
      eigenvalues <- pcoa$values$Corr_eig
      
      if (corr != "none")
        qual.FRic <- sum(eigenvalues[1:m]) / sum(eigenvalues)
      
      if (corr == "none") {
        # compute the eigenvalues (including negative ones)
        delta <- -0.5 * bicenter.wt(dist.matrix * dist.matrix)
        lambda <- eigen(delta, symmetric = TRUE, only.values = TRUE)$values
        # calculate corrected R^2-like ratio
        sum.m <- sum(lambda[1:m])
        sum.n <- sum(lambda)
        lambda.neg <- c(lambda[lambda < 0])
        max.neg <- abs(min(lambda.neg))
        qual.FRic <- (sum.m + (length(lambda[1:m]) * max.neg)) /
          (sum.n + ((length(lambda) - 1) * max.neg))
      }
      
      return(qual.FRic)
    }
    
    qual.FRic <- calc.qual.FRic(pcoa = pcoa, dist.matrix = dist.matrix, m = m)
  }
      
  # Clean up output
  out <- list(FRic = FRic, FEve = FEve, FDis = FDis, FDiv = FDiv, 
              qual.FRIC = qual.FRic)
  return(out)
}