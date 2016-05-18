##### parameterize for the models -------------------------------

#' @title assign coefficients used to generate parameters of the models
#' @description all the coefficients have their default values except \code{n1,n2}. All species are splitted into two groups. The intrinsic growth rates of species are assigned by groups. The case of all species in one group is achieved by considering species number of group 2 is zero.
#' \describe{
#'   \item{n1}{species number of group 1}
#'   \item{n2}{species number of group 2}
#'   \item{alpha.row.mu, alpha.row.sd}{the intrinsic growth rates of species in group 1}
#'   \item{alpha.col.mu, alpha.col.sd}{the intrinsic growth rates of species in group 2}
#'   \item{beta0.mu, beta0.sd}{the intra-species self-regulation competition in species}
#'   \item{beta1.mu, beta1.sd}{the interspecies competition between species}
#'   \item{gamma.mu, gamma.sd}{the interspecies mutualism between species}
#'   \item{h.mu, h.sd}{the Handling time, saturating coefficient of mutualistic interactions}
#'   \item{delta}{trade-off between strength and number of mutualistic interactions}
#' }
#' @return a list of coefficients
params_coeffs <- function(n1, n2,
                          alpha.row.mu = 0, alpha.row.sd = 0,
                          alpha.col.mu = 0, alpha.col.sd = 0,
                          beta0.mu = 0, beta0.sd = 0,
                          beta1.mu = 0, beta1.sd = 0,
                          gamma.mu = 0, gamma.sd = 0,
                          h.mu = 0, h.sd = 0, delta = 0) {
  list(n1 = n1, n2 = n2, alpha.row.mu = alpha.row.mu, alpha.row.sd = alpha.row.sd, alpha.col.mu = alpha.col.mu, alpha.col.sd = alpha.col.sd, beta0.mu = beta0.mu, beta0.sd = beta0.sd, beta1.mu = beta1.mu, beta1.sd = beta1.sd, gamma.mu = gamma.mu, gamma.sd = gamma.sd, h.mu = h.mu, h.sd = h.sd, delta = delta)
}

#' @title parameters for model \code{\link{model_lv2_cm}}
#' @description assign parameters for model \code{\link{model_lv2_cm}} according to a structural network and a couple of coefficients
#' @param graph, the adjacency matrix of mutualistic-competitive-mixed network, in which the competition part is marked by Negative values(-1), the cooperation part is marked by Positive values(1)
#' @param coeffs coefficients assigned in \code{\link{params_coeffs}}
#' @return a list of parameters for model \code{\link{model_lv2_cm}}
params_lv2_cm <- function(graph, coeffs) {
  with(coeffs, {
    n = n1 + n2
    C = graph # the competition part
    C[C > 0] = 0 # remove the mutualistic part
    C[C < 0] = 1 # assert the competitive part is only adjacency
    C = runif2(n * n, beta1.mu, beta1.sd) * C
    #C = matrix(0, nrow = n, ncol = n)
    #C[1:n1, 1:n1] = runif2(n1 * n1, beta1.mu, beta1.sd)
    #C[(n1+1):n, (n1+1):n] = runif2(n2 * n2, beta1.mu, beta1.sd)
    diag(C) = runif2(n, beta0.mu, beta0.sd)

    M = graph
    M[M < 0] = 0  # remove the mutualistic part
    edges = sum(M > 0)  # the number of all mutualistic interactions(edges)
    degrees = rowSums(M) # the degrees of all species
    M[M > 0] = runif2(edges, gamma.mu, gamma.sd) # values of interspecies cooperation
    M = M / degrees^delta  # trade-off of mutualistic strength and number

    h = runif2(n, h.mu, h.sd)
    r = c(runif2(n1, alpha.row.mu, alpha.row.sd), runif2(n2, alpha.col.mu, alpha.col.sd))
    list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
  })
}


#' @title parameters for model \code{\link{model_lv2_cm}}
#' @description being replaced by (\code{\link{params_lv2_bipartite_new}})
#' @param graph, the incident matrix of mutualistic networks which are bipartite
# params_lv2_bipartite <- function(graph, coeffs) {
#   s1 = dim(graph)[1]
#   s2 = dim(graph)[2]
#   s = s1 + s2
#   with(coeffs, {
#     C = matrix(0, nrow = s, ncol = s)
#     C[1:s1, 1:s1] = runif2(s1 * s1, beta1.mu, beta1.sd)
#     C[(s1+1):s, (s1+1):s] = runif2(s2 * s2, beta1.mu, beta1.sd)
#     diag(C) = runif2(s, beta0.mu, beta0.sd)
#
#     edges = sum(graph > 0)  # the number of edges
#     M = inc_to_adj(graph)  # transform to adjacency matrix (function of package [bipartite])
#     degrees = rowSums(M)
#     M[M > 0] = runif2(2 * edges, gamma.mu, gamma.sd) # values of interspecies cooperation
#     M = M / degrees^delta  # trade-off of mutualistic strength
#
#     h = runif2(s, h.mu, h.sd)
#     if(is.null(nstar)) {  # assign [r]
#       r = c(runif2(s1, alpha.row.mu, alpha.row.sd), runif2(s2, alpha.col.mu, alpha.col.sd))
#     } else {  # calculate [r] by model equation in equilibrium
#       r = C %*% nstar - (M %*% nstar) / (1 + h * M %*% nstar)
#       r = c(r)
#     }
#     list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
#   })
# }


