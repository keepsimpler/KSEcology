##########################Optimization Solve LV##############################


##########################Competition-Mutualism##############################

#' @title Lotka-Volterra (LV) Equations of Holling type II for a community mixed by Competition and Mutualism interactions
#' @param time, time step of simulation
#' @param init, the initial state of the LV system, a vector
#' @param params, parameters passed to LV model, a list of:
#' \describe{
#'   \item{r}{a vector of the intrinsic growth rates of species}
#'   \item{C}{a matrix of intra-species and inter-species competitions}
#'   \item{M}{a matrix of mutualism interactions among species}
#'   \item{h}{the saturate coefficient, handling time of species feed}
#' }
#' @return the derivation
#' @import deSolve
model_lv2_cm <- function(time, init, params, ...) {
  N = init  # initial state
  with(params, {
    dN <- N * ( r - C %*% N + (M %*% N) / (1 + h * M %*% N) )
    list(c(dN))
  })
}

#' @title parameters for model model_lv2_cm
#' @param graph, the incident matrix of mutualistic networks which are bipartite
#' @param alpha.row.mu, alpha.row.sd, the intrinsic growth rate of Plants (i.e., the ROW part)
#' @param alpha.col.mu, alpha.col.sd, the intrinsic growth rate of Animals (i.e., the COL part)
#' @param nstar, species abundances in equilibrium. When it is given, intrinsic growth rates are calculated rather then assigned.
#' @param beta0.mu, beta0.sd, the intraspecies competition
#' @param beta1.mu, beta1.sd, the interspecies competition
#' @param gamma.mu, gamma.sd, the interspecies cooperation
#' @param h.mu, h.sd, the Handling time, saturate coefficient
#' @param delta, trade-off between pairwise mutualistic interaction strengths
#' @return a list of parameters
params_lv2_bipartite <- function(graph, coeffs) {
  s1 = dim(graph)[1]
  s2 = dim(graph)[2]
  s = s1 + s2
  with(coeffs, {
    C = matrix(0, nrow = s, ncol = s)
    C[1:s1, 1:s1] = runif2(s1 * s1, beta1.mu, beta1.sd)
    C[(s1+1):s, (s1+1):s] = runif2(s2 * s2, beta1.mu, beta1.sd)
    diag(C) = runif2(s, beta0.mu, beta0.sd)

    edges = sum(graph > 0)  # the number of edges
    M = inc_to_adj(graph)  # transform to adjacency matrix (function of package [bipartite])
    degrees = rowSums(M)
    M[M > 0] = runif2(2 * edges, gamma.mu, gamma.sd) # values of interspecies cooperation
    M = M / degrees^delta  # trade-off of mutualistic strength

    h = runif2(s, h.mu, h.sd)
    if(is.null(nstar)) {  # assign [r]
      r = c(runif2(s1, alpha.row.mu, alpha.row.sd), runif2(s2, alpha.col.mu, alpha.col.sd))
    } else {  # calculate [r] by model equation in equilibrium
      r = C %*% nstar - (M %*% nstar) / (1 + h * M %*% nstar)
      r = c(r)
    }
    list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
  })
}

params_lv2_bipartite_new <- function(graph, coeffs) {
  with(coeffs, {
    s = s1 + s2
    C = matrix(0, nrow = s, ncol = s)
    C[1:s1, 1:s1] = runif2(s1 * s1, beta1.mu, beta1.sd)
    C[(s1+1):s, (s1+1):s] = runif2(s2 * s2, beta1.mu, beta1.sd)
    diag(C) = runif2(s, beta0.mu, beta0.sd)

    M = graph
    edges = sum(M > 0)  # the number of edges
    degrees = rowSums(M)
    M[M > 0] = runif2(edges, gamma.mu, gamma.sd) # values of interspecies cooperation
    M = M / degrees^delta  # trade-off of mutualistic strength

    h = runif2(s, h.mu, h.sd)
    r = c(runif2(s1, alpha.row.mu, alpha.row.sd), runif2(s2, alpha.col.mu, alpha.col.sd))
    list(r = r, C = C, M = M, h = h)  # the [parms] of ode model
  })
}


#' @title Simulate ODE dynamics of non-autonomous systems.A example is ecosystems under "press" perturbations. The dynamic is iteration of successive ODE dynamics of automous sytems (\code{\link{sim_ode_auto}}), while at each iterating step, the parameters and/or state values of systems are changed to reflect "press" perturbations.
#' @param perturb a function that change the parameters and state values after each iteration step
#' @param iter_steps possiblely maximum iteration steps
#' @param isout if output the transiting trajectory of each ODE iterate step
#' @param ... any arguments which are transfered to perturbation function
#' @return a list of lists :
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{params}{parameters assigned to the model}
#'   \item{species.survived}{a vector of species that survived}
#' }
sim_ode_press <- function(model, params, init, times, perturb, perturbNum = 500, isout = FALSE, extinct_threshold = 1e-10, ...) {
  ode.outs = list()
  for(i in 1:perturbNum) {
    ode.out = ode(init, times, model, params, method = "ode45", atol = 1e-14, rtol = 1e-12) #  method = "ode45", atol = 1e-13, rtol = 1e-13   # , rootfun = rootfun, method = "lsodar"
    nstar = as.numeric(ode.out[nrow(ode.out), 2:ncol(ode.out)]) # species biomass at equilibrium
    nstar[nstar < extinct_threshold] = 0  # species with biomass less than extinct threshold is considered to be extinct
    species.survived = which(nstar > 0)  # survived species

    flag = 0
    # if all species are extinct, will end the simulation
    if (length(species.survived) == 0) flag = 1
    # if any species' abundance is NaN, that means the ODE dynamic is unstable, the simulation will also be ended
    if (any(is.nan(nstar))) flag = 2

    Phi = jacobian.full(y = nstar, func = model, params = params) # community matrix, Jacobian matrix at equilibrium
    if (isout) {
      ret = list(out = ode.out, nstar = nstar, Phi = Phi, params = params, species.survived = species.survived, flag = flag)
    }
    else {
      ret = list(nstar = nstar, Phi = Phi, params = params, species.survived = species.survived, flag = flag)
    }
    ode.outs[[length(ode.outs) + 1]] = ret
    # if all species are extinct, end the simulation
    if (flag == 1 || flag == 2)
      break;

    # perturbation that returns new parameters and initial values
    perturb.ret = perturb(params, nstar, ...)
    params = perturb.ret$params
    init = perturb.ret$nstar
  }
  return(ode.outs)
}

rootfun <- function(time, init, params) {
  dstate <- unlist(model_lv2_cm(time, init, params))
  return(sum(abs(dstate)) - 1e-13)
}

#' @title calculate Jocobian matrix in equilibrium of model [model_lv2_cm]
#' @param params, parameters of model
#' @param nstar, species densities in equilibrium
#' @return Jacobian matrix
get_jacobian <- function(params, nstar = NULL) {
  with(params, {
    # Competitive part
    PhiC <- - diag(nstar) %*% C
    #PhiM <-  diag(c(nstar / (1 + h * M %*% nstar)^2)) %*% M
    PhiM <-  diag(nstar / ((1  + h * rowSums(M %*% diag(nstar)))^2)) %*% M
    Phi <- PhiC + PhiM
    return(Phi)
  })
}

#' @title analytically calculate [dot] and [ellipse] eigenvalue
#' @param s, species number
#' @param coeff, coefficient about
#' @param km, number of mutualistic neighbors
#' @param nstar, species abundance in equilibrium
#' @return [dot] eigenvalue and [ellipse] eigenvalue
get_jacobian_lev <- function(s, coeff, km, nstar) {
  # mean of the off-diagonal elements
  kc <- (s/2 - 1)  # number of competitive neighbors
  c <- coeff$beta1.mu  # strength of competitive interactions
  beta <- c * nstar
  #km <- 5  # number of mutualistic neighbors
  m <- coeff$gamma.mu  # strength of mutualistic interactions
  h <- coeff$h.mu
  gamma <- m * nstar / (1 + h * km * m * nstar)^2
  E <- (- kc * beta + km * gamma) / (s - 1)

  # variance of the off-diagonal elements
  E2 <- (s * kc * beta^2 + s * km * gamma^2) / (s * (s - 1))
  VAR <- E2 - E^2
  # the sample variance and the population variance have a scale relation
  # * (1 - 1/(s*(s-1))), So the population variance is:
  VAR2 <- VAR / (1 - 1/(s*(s-1)))

  # dot eigenvalue
  dot <- (s - 1) * E - coeff$beta0.mu * nstar
  # ellipse eigenvalue
  ellipse <- - E - coeff$beta0.mu * nstar + sqrt(s * VAR) * (1 + 1)
  c(dot = dot, ellipse = ellipse)
}
