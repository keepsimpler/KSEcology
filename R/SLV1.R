#####################################################################
# Simulation of  Stochastic LV1 model process
# dX = X (\Alpha + \Theta X) dt + X \Sigma dW
# X (\Alpha + \Theta X) is the drift vector, which include:
#   \Alpha - vector of intrinsic growth rates
#   \Theta - matrix of (mutualistic/competitive) interactions
# X \Sigma is the diffusion matrix, we assume it's a diagonal matrix
# W is a vector of Wiener process


library(yuima)
library(plyr)


#' @title construct the drift vector of SLV1 process according to the dimension
#' @param m dimension of variables
#' @return drift vector of SLV1 process
#' @example c("alpha1*x1+theta11*x1*x1+theta12*x1*x2",
#'            "alpha2*x2+theta21*x2*x1+theta22*x2*x2")
set_drift_slv1 <- function(m) {
  sapply(1:m, function(row) {
    theta <- paste('theta', row, 1:m, '*x', row, '*x', 1:m, sep = '', collapse = '+')
    paste('alpha', row, '*', 'x', row, '+', theta, sep = '')
  })
}

#' @title construct the diffusion matrix of SLV1 process according to the dimension
#' @param m dimension of variables
#' @return diffusion matrix of SLV1 process
#' @example c('sigma1*x1', '0', '0', 'sigma2*x2')
set_diffusion_slv1_geometric <- function(m) {
  apply(diag(1:m), c(1, 2), function(ij)
    if (ij != 0)
      paste('sigma', ij, '*x', ij, sep = '')
    else '0'
  )
}
set_diffusion_slv1_wiener <- function(m) {
  apply(diag(1:m), c(1, 2), function(ij)
    if (ij != 0)
      paste('sigma', ij, sep = '')
    else '0'
  )
}

#' @title set true values for the coefficients in drift and diffusion of SLV1 process
#' @param m dimension of variables
#' @param Alpha, Theta drift coefficient matrix
#' @param Sigma diffusion coefficient matrix(a diagonal matrix)
#' @return a coefficient list, whose names setted by \code{set_drift_slv1} and \code{set_diffusion_slv1}
#' @example list(alpha1 = 1, alpha2 = 1, theta11 = -1, theta21 = 0.1, theta12 = 0.2, theta22 = -1,
#' sigma1 = 0.1, sigma2 = 0.1)
set_true_parameters_slv1 <- function(m, Alpha, Theta, Sigma) {
  params.alpha <- as.list(Alpha)
  names(params.alpha) <- paste('alpha', 1:m, sep = '')
  # transform drift coefficient matrix [Theta] to list
  params.theta <- as.list(Theta)
  drift.ij <- outer(1:m, 1:m, FUN = paste, sep="")
  names(params.theta) <- sapply(drift.ij, function(ij) paste('theta', ij, sep = ''))
  params.sigma <- as.list(diag(Sigma))
  names(params.sigma) <- paste('sigma', 1:m, sep = '')
  c(params.alpha, params.theta, params.sigma)
}


sim_slv1 <- function(m, Alpha, Theta, Sigma, diffusion.type = 'wiener', Xinit, steps = 10000, stepwise = 0.01) {
  # set model
  drift = set_drift_slv1(m)
  if (diffusion.type == 'wiener') {
    diffusion = set_diffusion_slv1_wiener(m)
  }
  else if (diffusion.type == 'geometric') {
    diffusion = set_diffusion_slv1_geometric(m)
  }
  else
    stop('Please assign a proper diffusion type: wiener process or geometric brownian motion.')
  solve.variable = set_solve_variables(m)
  mod.slv1 = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, state.variable = solve.variable)
  # set sampling
  samples = setSampling(Terminal = steps * stepwise, n = steps)
  # initialize yuima object
  yuima.slv1 <- setYuima(model = mod.slv1, sampling = samples)
  # assign values for parameters
  parameters = set_true_parameters_slv1(m, Alpha, Theta, Sigma)
  # simulate and update yuima object
  yuima.slv1 = simulate(yuima.slv1, true.parameter = parameters, xinit = Xinit)
  yuima.slv1
}


