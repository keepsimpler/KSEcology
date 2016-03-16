#' @title perturbations that effect on species by increasing/decreasing the intrinsic growth rates of all species
#' @param params parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param r.delta deviation of intrinsic growth rates at each iterating step
perturb_growthrate <- function(params, nstar, r.delta.mu = 0.01, r.delta.sd = 0.01) {
  #set.seed(1)
  params$r = params$r - runif2(length(params$r), r.delta.mu, r.delta.sd)
  list(params = params, nstar = nstar)
}

#' @title perturbations that effect on species by increasing/decreasing the intrinsic growth rates of a part of species
#' @param params parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param r.delta deviation of intrinsic growth rates at each iterating step
#' @param perturbed_species the index of perturbed species
perturb_growthrate_part <- function(params, nstar, r.delta.mu = 0.01, r.delta.sd = 0.01,  perturbed_species) {

  params$r[perturbed_species] = params$r[perturbed_species] - runif(length(perturbed_species), min = r.delta.mu - r.delta.sd, max = r.delta.mu + r.delta.sd)
  list(params = params, nstar = nstar)
}

#' @title perturbations that effect on mutualistic interactions by increasing/decreasing strengths of them
#' @param params parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param gamma.delta deviation of mutualistic interaction strengths at each iterating step
perturb_mutualistic_strength <- function(params, nstar, gamma.delta.mu = 0.01, gamma.delta.sd = 0.01) {
  edges = length(params$M[params$M > 0])
  params$M[params$M > 0] = params$M[params$M > 0] - runif(edges, min = gamma.delta.mu - gamma.delta.sd, max = gamma.delta.mu + gamma.delta.sd)
  list(params = params, nstar = nstar)
}

#' @title perturbations that remove one species
#' @param params parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param extinct_species the removed species
perturb_primary_extinct <- function(params, nstar, extinct_species) {
  nstar = nstar[- extinct_species]  # primary extinction
  params$r = params$r[- extinct_species]
  params$C = params$C[- extinct_species, - extinct_species]
  params$M = params$M[- extinct_species, - extinct_species]
  params$h = params$h[- extinct_species]
  list(params = params, nstar = nstar)
}

