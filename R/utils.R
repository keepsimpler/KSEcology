#' @title Emulate ggplot2 default color palette
#' @param n number of colors
#' @references http://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_palette_default <- function(n) {
  require(graphics)
  palette('default') # reset back to the default
  rep(palette(), ceiling(n / length(palette())))[1:n]
}

#' @title another form of uniform distribution between [mean - sd, mean + sd]
runif2 <- function(n, mean, sd) {
  runif(n) * 2 * sd + (mean - sd)
}

#' @title transfer an incidence matrix to an adjacency matrix
#' @param inc, an incidence matrix
#' @return adj, an adiacency matrix
inc_to_adj <- function(inc){
  p <- dim(inc)[1]  # number of Plants
  a <- dim(inc)[2]  # number of Animals
  s <- p + a  # number of all Species
  adj <- matrix(0, s, s)  # initialize the adjacency matrix as a zero-matrix
  adj[1:p, (p + 1):s] <- inc  # the upper right sub-matrix is the incidence matrix
  adj <- adj + t(adj)  # the lower left sub-matrix is transpose of the incidence matrix
  return(adj)
}

#' @title get stability measurements according to the Jacobian at equilibrium
#' @param Phi, the Jacobian at equilibrium
#' @return a list of stability measurements:
#'   Rinf -- Asymptotic resilience (largest eigenvalue)
#'   R0 -- Initial resilience (reactivity)
#'   Is -- Stochastic Invariability
#'   Id -- Deterministic Invariability
#'   R0 < Is < Id < Rinf
#' @references 'Resilience, reactivity and variability :A mathematical comparison of ecological stability measures'
get_stability_measurements <- function(Phi) {
  s = dim(Phi)[1]
  I = diag(1, s)
  Rinf = - max(Re(eigen(Phi)$values))
  Is = 1 / (2 * norm(- solve(kronecker(I, Phi) + kronecker(Phi, I)), type = '2'))
  Id = 1 / norm(- solve(Phi), type = '2')
  R0 = - max(Re(eigen(Phi + t(Phi))$values)) / 2
  c(Rinf = Rinf, Is = Is, Id = Id, R0 = R0)
}

#' @title get covariance matrix of multivariate OU process
#' @param Phi, community matrix
#' @param Sigma, the covariance matrix of environmental fluctuating
get_covariance_mou <- function(Phi, Sigma) {
  s = dim(Phi)[1]
  I = diag(1, s)
  - matrix(solve(kronecker(I, Phi) + kronecker(Phi, I)) %*% as.vector(Sigma), nrow = s, ncol = s)
}

