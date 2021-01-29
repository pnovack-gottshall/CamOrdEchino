## Modification of ecospace::calc_metrics() to (1) allow direct import of
## pre-calculated distance matrix, character matrix (including inferred
## ancestral states), and PCoA output, and passage to stream-lined dbFD2. (2)
## Includes modification for calculation of 'H' (life habit/morphotype richness)
## using distance matrix (instead of character matrix), which offers improved
## handling of missing data. (3) Adds arg to standardize PCoA eigenvectors
## according to magnitude of eigenvalues. (4) Streamlines error-checking to make
## faster. Code by Phil Novack-Gottshall and taken from
## https://github.com/pnovack-gottshall/ecospace/blob/master/R/calc_metrics.R.

calc_metrics2 <- function(sample = NA,
                          dist.sam = NULL,
                          pcoa = NULL,
                          m = 3,
                          stand.pcoa = TRUE,
                          calc.FRic.and.FDiv = TRUE) {
  s <- nrow(sample)
  odir <-
    setwd(tempdir())     # Specify the pre-built (and CPU-process unique) temp directory for storage of vert.txt temp files for convex hull calculations
  on.exit(setwd(odir))
  sam.out <- data.frame(S = integer(1), H = integer(1), D = numeric(1), 
                        M = numeric(1), V = numeric(1), FRic = numeric(1), 
                        FEve = numeric(1), FDiv = numeric(1), FDis = numeric(1),
                        qual.FRic = numeric(1))
  sam <- sample[1:s,]
  sam.out$S <- s
  H <- nrow(unique(dist.sam))
  sam.out$H <- H
  if (any(is.nan(dist.sam)) | length(dist.sam) == 0)
    next
  sam.out$D <- mean(dist.sam, na.rm = TRUE)
  sam.out$M <- max(dist.sam, na.rm = TRUE)
  sam.out$V <- sqrt(sum(apply(sam, 2, var, na.rm = TRUE), na.rm = TRUE))
  if (s <= m | H <= m)
    next
  FD <- simple.dbFD(pcoa = pcoa, dist.matrix = dist.sam, m = m, 
                    stand.pcoa = stand.pcoa, 
                    calc.FRic.and.FDiv = calc.FRic.and.FDiv)
  sam.out$FRic <- FD$FRic
  sam.out$FEve <- FD$FEve
  sam.out$FDiv <- FD$FDiv
  sam.out$FDis <- FD$FDis
  sam.out$qual.FRic <- FD$qual.FRIC
  return(sam.out)
}
