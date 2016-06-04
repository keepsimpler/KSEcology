#####################################################################
# Simulation of  Stochastic LV2 model process
# Note: we support the bipartite network of mutualistic interactions
# dX = X (\Alpha + \Theta X + \Gamma X / (1 + h \Gamma X)) dt + X \Sigma dW
# X (\Alpha + diag(\Theta) X + \Theta X / (1 + h \Theta X)) is the drift vector,
# which include:
#   \Alpha - vector of intrinsic growth rates
#   \Theta - matrix of (diagonal) competitive interactions and (off-diagonal) mutualistic interactions
#   h - handling time
# Note: we only support intra-species competition, do not support inter-species competition
# X \Sigma is the diffusion matrix, we assume it's a diagonal matrix
# W is a vector of Wiener process


library(yuima)
library(plyr)

#' @title construct the drift vector of SLV2 process according to the dimension
#' @param m dimension of variables
#' @return drift vector of SLV2 process
#' @example c("alpha1*x1+theta11*x1*x1+theta12*x1*x2/(1+h*theta12*x2)",
#'            "alpha2*x2+theta21*x2*x1/(1+h*theta21*x1)+theta22*x2*x2")
set_drift_slv2 <- function(m) {
  rows <- 1:m
  sapply(rows, function(row) {
    theta <- paste('theta', row, rows[-row], '*x', rows[-row], sep = '', collapse = '+')
    paste('alpha', row, '*', 'x', row, '+', 'theta', row, row, '*x', row,
          '*x', row, '+x', row, '*(', theta, ')/(1+h*(', theta,')', sep = '')
  })
}


getParameters <- function(m, Alpha, Theta, Sigma, h) {
  parameters = list()
  for(i in 1:m) {
    for(j in 1:m) {
      paraname = paste('theta', i, j, sep = '')
      parameters[paraname] = Theta[i, j]
    }
  }
  for(i in 1:m) {
    paraname = paste('sigma', i, sep = '')
    parameters[paraname] = Sigma[i, i]
  }
  for(i in 1:m) {
    paraname = paste('alpha', i, sep = '')
    parameters[paraname] = Alpha[i]
  }
  parameters['h'] = h
  parameters
}

# return a 3 dimensional array:
# first dimension: simnum, number of simulations
# second dimension: t+1, the simulating steps of one simulation
# third dimension: m, the dimensions of multivariates
sim.slv2 <- function(simnum = 1000, steps = 1000, stepwise = 0.01, m, Alpha, Theta, Sigma, h = 0.1, Xinit) {
  grid = setSampling(Terminal = steps * stepwise, n = steps)
  drift = getDrift(m)
  diffusion = getDiffusion(m)
  solve.variable = getSolveVariables(m)
  mod3 = setModel(drift = drift, diffusion = diffusion, solve.variable = solve.variable, xinit = Xinit)
  parameters = getParameters(m, Alpha, Theta, Sigma, h)
  #X3 = simulate(mod3, true.parameter = parameters, sampling = grid)
  Xs = laply(1:simnum, .parallel = TRUE, function(i) {
    print(i)
    X = simulate(mod3, true.parameter = parameters, sampling = grid)
    matrix(X@data@original.data, ncol = m)
  })
  Xs
}


# m = 3
# Alpha = c(1, 1, 1)
# Theta = matrix(c(1, -0.1, -0.1, -0.1, 1, -0.1, -0.1, -0.1, 1), ncol = m)
# Sigma = diag(0.05, m)
# Xinit = c(1, 1, 1)
# h = 0.01
# slv2.out = sim.slv2(simnum = 1, steps = 1000, m = m, Alpha = Alpha, Theta = Theta, Sigma = Sigma, h = h, Xinit = Xinit)
# matplot(slv2.out, type = 'l', lwd = 0.7)
