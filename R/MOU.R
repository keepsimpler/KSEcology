#####################################################################
# Simulation of Multivariate OU (Ornstein - Uhlenbeck) process
# dX = \Theta X dt + \Sigma dW
# \Theta is the drift vector,
# \Sigma is the diffusion matrix, we assume it's a diagonal matrix
# W is a vector of Wiener process


library(yuima)
library(plyr)

# library(doMC)  #
# registerDoMC(20)  # register Multi Cores
# getDoParWorkers()  # get available Cores

#' @title construct names of solve variables according to the dimension
#' @param m dimension of variables
#' @return a vector of characters
#' @example c('x1', 'x2')
set_solve_variables <- function(m) {
  sapply(1:m, function(i) paste('x', i, sep = ''))
}

#' @title construct the drift vector of MOU process according to the dimension
#' @param m dimension of variables
#' @return drift vector of MOU process
#' @example c('theta11 * x1 + theta12 * x2', 'theta21 * x1 + theta22 * x2')
set_drift_mou <- function(m) {
  drift.ij <- outer(1:m, 1:m, FUN = paste, sep="")
  apply(drift.ij, 1, function(row) {
    drift <- paste('theta', row, '*x', 1:m, sep = '', collapse = '+');
    })
}

#' @title construct the diffusion matrix of MOU process according to the dimension
#' @param m dimension of variables
#' @return diffusion matrix of MOU process
#' @example c('sigma1', '0', '0', 'sigma2')
set_diffusion_mou <- function(m) {
  apply(diag(1:m), c(1, 2), function(ij)
    if (ij != 0)
      paste('sigma', ij, sep = '')
    else '0'
    )
}

#' @title set true values for the coefficients in drift and diffusion of MOU process
#' @param m dimension of variables
#' @param Theta drift coefficient matrix
#' @param Sigma diffusion coefficient matrix(a diagonal matrix)
#' @return a coefficient list, whose names setted by \code{set_drift_mou} and \code{set_diffusion_mou}
#' @example list(theta11 = -1, theta21 = 0.1, theta12 = 0.2, theta22 = -1,
#' sigma1 = 0.1, sigma2 = 0.1)
set_true_parameters_mou <- function(m, Theta, Sigma) {
  # first, transform drift coefficient matrix [Theta] to list
  params.theta <- as.list(Theta)
  drift.ij <- outer(1:m, 1:m, FUN = paste, sep="")
  names(params.theta) <- sapply(drift.ij, function(ij) paste('theta', ij, sep = ''))
  params.sigma <- as.list(diag(Sigma))
  names(params.sigma) <- paste('sigma', 1:m, sep = '')
  c(params.theta, params.sigma)
}

#' @title multivariate OU process simulation
#' @param m, the number(dimension) of variables
#' @param Theta, the drift coefficient matrix
#' @param Sigma, the diffusion coefficient matrix
#' @param Xinit, the initialized values of variables
#' @param simnum, number of simulations
#' @param steps, number of time steps
#' @param stepwise, length of one step
#' @return a 3 dimensional array:
#'         first dimension: simnum, number of simulations
#'         second dimension: steps + 1, the simulated steps of one simulation
#'         third dimension: m, the number of variables
sim_mou <- function(m, Theta, Sigma, Xinit, simnum = 1, steps = 10000, stepwise = 0.01) {
  grid = setSampling(Terminal = steps * stepwise, n = steps)
  drift = set_drift_mou(m)
  diffusion = set_diffusion_mou(m)
  solve.variable = set_solve_variables(m)
  mod = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, xinit = Xinit)
  parameters = set_true_parameters_mou(m, Theta, Sigma)
  Xs = laply(1:simnum, .parallel = FALSE, function(i) {
    print(i)
    X = simulate(mod, true.parameter = parameters, sampling = grid)
    matrix(X@data@original.data, ncol = m)
  })
  Xs
}



# m = 1
# Theta = matrix(-1)
# Sigma = matrix(0.8)
# Xinit = c(0)
# mou.out = sim_mou(simnum = 1, steps = 100000, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)

# m = 2
# Theta = matrix(c(-1, 0.1, 0.8, -1), ncol = m)
# Sigma = diag(0.8, m)
# Xinit = c(0, 0)
# mou.out = sim_mou(simnum = 1, steps = 100000, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# var = diag(Sigma)^2 / 2 * solve(Theta)  # expected variance-covariance matrix for the symmetric [Theta]

# m = 3
# Theta = matrix(c(-1, 0.1, 0.1, 0.8, -1, 0.1, 0.1, 0.1, -1), ncol = m)
# Sigma = diag(0.8, m)
# Xinit = c(0, 0, 0)
# mou.out = sim_mou(simnum = 1, steps = 10000, stepwise = 0.01, m = m, Theta = Theta, Sigma = Sigma, Xinit = Xinit)
# matplot(mou.out[2,,], type = 'l')
# mou.XMeans = aaply(mou.out, .margins = c(2, 3), mean)
# mou.XVars = aaply(mou.out, .margins = c(2), var)
# mou.XVars.self = aaply(mou.XVars, .margins = c(1), function(XVar) {
#   c(diag(XVar), XVar[lower.tri(XVar)])
# })
# matplot(mou.XVars.self, type = 'l')


# grid = setSampling(Terminal = 10, n = 10000)
# simnum = 100
#
# ## One dimensional OU process
# mod = setModel( drift = "mu - theta * x", diffusion = "sigma", state.var = "x", time.var = "t", solve.var = "x")
# X = simulate(mod, xinit = 0, true.parameter = list(theta = 1.5, sigma = 0.5, mu = 1), sampling = grid)

# ## Two dimensional OU process
# sol = c('x1', 'x2')
# drift = c('-theta11 * x1 - theta12 * x2', '-theta21 * x1 - theta22 * x2')
# diffusion = matrix(c('sigma1', '0', '0', 'sigma2'), 2, 2)
# mod2 = setModel(drift = drift, diffusion = diffusion, solve.variable = sol)
# X2 = simulate(mod2, true.parameter = list(theta11 = 1, theta12 = 0.1, theta21 = 0.1, theta22 = 1, sigma1 = 1, sigma2 = 1),
#               sampling = grid)
