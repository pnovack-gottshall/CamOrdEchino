## Modification of ecospace::calc_metrics() to (1) allow direct import of
## pre-calculated distance matrix, character matrix (including inferred
## ancestral states), and PCoA output, and passage to stream-lined dbFD2. (2)
## Includes modification for calculation of 'H' (life habit/morphotype richness)
## using distance matrix (instead of character matrix), which offers improved
## handling of missing data. (3) Adds arg to standardize PCoA eigenvectors
## according to magnitude of eigenvalues. (4) Streamlines error-checking to make
## faster. Code by Phil Novack-Gottshall and taken from
## https://github.com/pnovack-gottshall/ecospace/blob/master/R/calc_metrics.R.
## (5) Adds sum of total ranges (using option to use Dineen, et al. (2019)
## standardization of PCoA eigenvectors).

calc_metrics2 <- function(sample = NA,
                          dist.sam = NULL,
                          pcoa = NULL,
                          m = 3,
                          stand.pcoa = TRUE,
                          calc.FRic.and.FDiv = TRUE) {

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
  
  s <- nrow(sample)
  odir <-
    setwd(tempdir())     # Specify the pre-built (and CPU-process unique) temp directory for storage of vert.txt temp files for convex hull calculations
  on.exit(setwd(odir))
  sam.out <- data.frame(S = NA, H = NA, D = NA, M = NA, V = NA, R = NA, FRic = NA, FEve = NA, FDiv = NA, FDis = NA, qual.FRic = NA)
  sam <- sample[1:s,]
  sam.out$S <- s
  H <- nrow(unique(dist.sam))
  sam.out$H <- H
  if (!any(is.nan(dist.sam)) | length(dist.sam) != 0) {
    sam.out$D <- mean(dist.sam, na.rm = TRUE)
    sam.out$M <- max(dist.sam, na.rm = TRUE)
    sam.out$V <- sqrt(sum(apply(sam, 2, var, na.rm = TRUE), na.rm = TRUE))
    
    # Confirm 'pcoa' exists before calculating sum of ranges
    if (!is.null(pcoa)) {
  
      # Get PCoA components (with switch for whether corrected or not)
      if (pcoa$correction[1] == "none") {
        vectors <- pcoa$vectors
        eigenvalues <- pcoa$values$Eigenvalues
      } else {
        vectors <- pcoa$vectors.cor
        eigenvalues <- pcoa$values$Corr.eig
      }
      
      # Standardize the eigenvectors according to magnitude of eigenvalues (if desired)
      if (stand.pcoa)
        vectors <- stand.pcoa.fn(vectors = vectors, eigenvalues = eigenvalues)
    
      maxes <- apply(vectors, 2, max, na.rm = TRUE)
      mins <- apply(vectors, 2, min, na.rm = TRUE)
      sam.out$R <- sum(maxes - mins, na.rm = TRUE)
    }
    
    # Move on to functional diversity metrics
    if (s > m | H > m) {
      FD <- simple.dbFD(pcoa = pcoa, dist.matrix = dist.sam, m = m, 
                        stand.pcoa = stand.pcoa, 
                        calc.FRic.and.FDiv = calc.FRic.and.FDiv)
      sam.out$FRic <- FD$FRic
      sam.out$FEve <- FD$FEve
      sam.out$FDiv <- FD$FDiv
      sam.out$FDis <- FD$FDis
      sam.out$qual.FRic <- FD$qual.FRIC
    }
  }
  return(sam.out)
}
