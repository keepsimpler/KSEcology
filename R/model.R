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


model_lv1 <- function(time, init, params, ...) {
  N = init  # initial state
  with(params, {
    dN <- N * ( r - C %*% N + (M %*% N) )
    list(c(dN))
  })
}

#' @title Simulate ODE dynamics of non-autonomous systems.A example is ecosystems under "press" perturbations. The dynamic is iteration of successive ODE dynamics of automous sytems, while at each iterating step, the parameters and/or state values of systems are changed to reflect "press" perturbations.
#' @param perturb a function that change the parameters and state values after each iteration step
#' @param perturbNum possiblely maximum iteration steps
#' @param method of ODE simulation
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
sim_ode_press <- function(model, params, init, times, perturb, perturbNum = 500, method = 'lsoda', isout = FALSE, extinct_threshold = 1e-10, ...) {
  ode.outs = list()
  for(i in 1:perturbNum) {
    ode.out = ode(init, times, model, params, method = method, atol = 1e-14, rtol = 1e-14) #  method = "ode45", atol = 1e-13, rtol = 1e-13   # , rootfun = rootfun, method = "lsodar"
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

# rootfun <- function(time, init, params) {
#   dstate <- unlist(model_lv2_cm(time, init, params))
#   return(sum(abs(dstate)) - 1e-13)
# }

