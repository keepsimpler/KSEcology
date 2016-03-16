#' @title the second raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu2, the second raw moment \mu^2
Mu2 <- function(r,n) {
  mu2 =  r^2 * factorial(2/n) * gamma(3/n) /
    (factorial(1/n) * gamma(4/n))
  return(mu2)
}

#' @title the fourth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu4, the fourth raw moment \mu^4
Mu4 <- function(r,n) {
  mu4 = 8 * r^4 * 2 * gamma(2/n) * gamma(5/n) /
    (3 * gamma(1/n) * gamma(6/n))
  return(mu4)
}

#' @title the sixth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu6, the sixth raw moment \mu^6
Mu6 <- function(r,n) {
  mu4 = 8 * r^6 * 2 * gamma(2/n) * gamma(7/n) /
    ( gamma(1/n) * gamma(8/n))
  return(mu4)
}


## Copyed from Allesina's code
PDFSuperEllipseSymmetric <- function(x, r, n){
  LogDenominator <- log(2) + (1/n) * log((r)^n - abs(x)^n)
  LogNumerator <- log(4) + 2 * log(r) + 2 * lgamma(1 + 1 / n) - lgamma(1 + 2/n)
  res <- (exp(LogDenominator - LogNumerator))
  res[is.nan(res)] <- 0
  return(res)
}

FindSECounts <- function(Hbreaks, Hcounts, r, n){
  data.points <- length(Hcounts)
  for (i in 1:data.points){
    Hcounts[i] <- integrate(PDFSuperEllipseSymmetric, Hbreaks[i], Hbreaks[i+1], r, n)$value
  }
  return(Hcounts)
}
