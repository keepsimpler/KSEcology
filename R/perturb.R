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
perturb_mutualistic_strength <- function(params, nstar, gamma.delta.mu = 0.01, gamma.delta.sd = 0.0, nstar.sd = 0) {
  edges = length(params$M[params$M > 0])
  params$M[params$M > 0] = params$M[params$M > 0] + runif(edges, min = gamma.delta.mu - gamma.delta.sd, max = gamma.delta.mu + gamma.delta.sd)
  # random perturb nstars
  if (nstar.sd > 0)
    nstar <- nstar * runif2(length(nstar), 1, nstar.sd)
  list(params = params, nstar = nstar)
}

#' @title perturbations that effect on competitive interactions by increasing/decreasing strengths of them
#' @param params parameters assigned to the ODE model
#' @param nstar state values at equilibrium
#' @param beta1.delta deviation of competitive interaction strengths at each iterating step
#' @param nstar.immigration, simulating a immigration factor to (re)establish extinct species
#' @param nstar.sd, randomly proportionally perturb species abundances at the end of each step of ode simulation to simulate(explore) the Initial Value Problem and to break the problem of saddle point(if the initial values are Very close to the equilibrium values, then the system may keep in equilibrium even it is not local stable. Because simulation error??)
perturb_competitive_strength <- function(params, nstar, beta1.delta.mu = 0.01, beta1.delta.sd = 0.0, nstar.immigration = 0, nstar.sd = 0) {
  beta0 <- diag(params$C) # backup the self-regulation
  diag(params$C) <- 0 # remove the self-regulation
  edges = length(params$C[params$C > 0])
  params$C[params$C > 0] = params$C[params$C > 0] + runif(edges, min = beta1.delta.mu - beta1.delta.sd, max = beta1.delta.mu + beta1.delta.sd)
  diag(params$C) <- beta0 # restore the self-regulation
  # reestablish extinct species by simulating a immigration factor
  if (nstar.immigration > 0)
    nstar[nstar == 0] <- nstar[nstar == 0] + nstar.immigration
  # random perturb nstars
  if (nstar.sd > 0)
    nstar <- nstar * runif2(length(nstar), 1, nstar.sd)
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

