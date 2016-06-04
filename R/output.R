##### process the output of ODE simulation---------------------------

#' @title plot the output of ode simulation
#' @param the output of ode simulation
#' \describe{
#'   \item{out}{output of one ODE simulation, including the trajectory of values of state variables}
#'   \item{nstar}{the values of state variables in equilibrium}
#'   \item{Phi}{the Jacobian matrix in equilibrium}
#'   \item{params}{parameters assigned to the model}
#'   \item{species.survived}{a vector of species that survived}
#' }
plot_ode_output_1 <- function(out) {
  out.nstars = laply(out, function(one) {
    one$nstar
  })
  matplot(out.nstars, type = 'l', lwd = 1.8, xlab = 'Steps of gradual pressure', ylab = 'Abundances of species')
}

#' @title ggplot function to plot equilibrium trajectory under press
#' @param out output of press simulation
ggplot_ode_output_1 <- function(out, title = NULL, reverse = FALSE) {
  out.nstars = laply(out, function(one) {
    one$nstar
  })
  n <- dim(out.nstars)[2]
  out.nstars <- data.frame(out.nstars)
  require(reshape2)
  # reshape to a long format
  out.nstars$c <- sapply(out, function(one) one$params$C[1,2])
  out.nstars.long <- melt(data = out.nstars, id.vars = 'c')
  p <-
  ggplot(out.nstars.long, aes(x = c,y = value, group = variable, colour = variable, lty = variable)) +
    geom_line(size = 1.) + theme_bw() +
#    scale_colour_manual(name = 'ews', values = gg_color_hue(n)) +
    labs(title = title, x = 'Competitive Strength(c)', y = 'Species  Abundances') +
    theme(legend.position = "none")
  if (reverse == TRUE) {
    p <- p + scale_x_reverse()
  }
  p
}

plot_ode_output_2 <- function(out, step) {
  ode.out <- out[[step]]$out
  matplot(ode.out[, -1], type = 'l', lwd = 1.8, xlab = 'Time', ylab = 'Abundances of species')

}

#' @title ggplot function to plot state trajactory under fixed parameters
#' @param out output of press simulation
#' @param step a step which fix the parameters
ggplot_ode_output_2 <- function(out, step, title = NULL, timesteps = NULL) {
  ode.out <- out[[step]]$out
  ode.out <- data.frame(ode.out)
  if (!is.null(timesteps) && timesteps > nrow(ode.out)) {
    stop('Given time steps larger than the maximam value')
  }
  else {
    ode.out <- ode.out[1:timesteps, ]
  }
  require(reshape2)
  # reshape to a long format
  ode.out.long <- melt(data = ode.out, id.vars = 'time')
  ggplot(ode.out.long, aes(x = time,y = value, group = variable, colour = variable, lty = variable)) +
    geom_line(size = 1.) + theme_bw() +
    labs(title = title, x = 'Time  Steps', y = 'Species  Abundances') +
    theme(legend.position = "none")
}

get_nstars <- function(out, steps) {
  out.nstars = laply(out[steps], function(one) {
    one$nstar
  })
  return(out.nstars)
}

get_jacobian <- function(out, step) {
  params <- out[[step]]$params
  nstar <- out[[step]]$nstar
  get_jacobian_inner(params, nstar)
}

#' @title calculate Jocobian matrix in equilibrium of model [model_lv2_cm]
#' @param params, parameters of model
#' @param nstar, species densities in equilibrium
#' @return Jacobian matrix
get_jacobian_inner <- function(params, nstar = NULL) {
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
