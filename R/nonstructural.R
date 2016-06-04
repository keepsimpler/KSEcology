#' @title calculate several measurements of stability
#'
#' @example calc_measurements(rmax = 0, m = 1.5, s = 1, h = 0.1, sigma = 0.1, steps = 1000)
calc_measurements <- function(rmax, m, s, h, sigma = 0.1, steps = 1000) {
  rmin = - (sqrt(m) - sqrt(s))^2 / (h * m)
  xmin = sqrt(s) * (sqrt(m) - sqrt(s)) / (s * h * m)
  #xr0 <- -(s - m) / (s * h * m)
  r = seq(from = rmax, to = rmin, length.out = steps + 1)
  A = s * h * m
  B = s - m - r * h * m
  C = - r
  delta = B^2 - 4 * A * C
  X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
  X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
  J =  (- s + m / (1 + h * m * X1)^2) # Jm + Jc
  # calculate the largest eigenvalue
  lambda1 <- X1 * J
  # calculate sensitivity
  dxdr <- - 1 / J
  dr2dx2 <- (2 * h * m^2) / (1 + h * m * X1)^3
  # calculate variance
  variance.wiener <- - sigma^2 / 2 / lambda1
  variance.geometric <- - (sigma * X1)^2 / 2 / lambda1
  #variance.geometric.2 <- - (sigma * X1^2)^2 / 2 / lambda1
  basin <- sqrt(B^2 - 4 * A * C) / A
  data.frame(r = r, X1 = X1, X2 = X2, J = J, lambda1 = lambda1, dxdr = dxdr, dr2dx2 = dr2dx2, variance.wiener = variance.wiener, variance.geometric = variance.geometric, basin = basin)
}

#' @title Simulation of stochastic process of the nonstructural model
#' @param sim.type simulation type: 'transient', 'stationary'
sim_nonstructural_model_2 <- function(r, s, m, h, sim.type = 'transient', sigma, diffusion.type = 'wiener', xinit = NULL, steps = 10000, stepwise = 0.01) {
  if (sim.type == 'stationary') {
    delta = 0
  }
  else if (sim.type == 'transient') {
    rmax = r
    rmin = - (sqrt(m) - sqrt(s))^2 / (h * m)
    terminal_time = steps * stepwise  # time interval
    delta = (rmax - rmin) / terminal_time # increment of parameter at each step
  }
  else {
    stop('Only support \'transient\' and \'stationary\' simulation')
  }
  sim_nonstructural_model(r = r, s = s, m = m, h = h, delta = delta, sigma = sigma, diffusion.type = diffusion.type, xinit = xinit, steps = steps, stepwise = stepwise)
}

#' @title Simulation of stochastic process of the nonstructural model
#' @description dx = x(r - delta t - s x + m x / (1 + h m x) ) dt + \sigma dw
#' @param r intrinsic growth rate, s self-regulation, m mutualistic strength, h handing time,
#' @param delta (delta = 0, stationary simulation; delta > 0, transient simulation, delta = (r - rmin) / (steps * stepwise)), sigma diffusion
#' @param diffusion.type 'wiener' - \sigma dw, 'geometric' - \sigma x dw
#' @param xinit steps stepwise
#' @import yuima package
sim_nonstructural_model <- function(r, s, m, h, delta, sigma, diffusion.type = 'wiener', xinit = NULL, steps = 10000, stepwise = 0.01) {
  # a flaw: when h = 0 or m = 0, then A = 0, xinit can not be calculated
  if (is.null(xinit)){
    A = s * h * m
    B = s - m - r * h * m
    C = - r
    stopifnot(B^2 - 4 * A * C > 0)
    xinit <- (- B + sqrt(B^2 - 4 * A * C)) / (2 * A)
  }
  # set model
  drift = 'x * (r - delta * t - s * x + m * x / (1 + h * m * x))'
  if (diffusion.type == 'wiener') {
    #drift = 'x * (r - s * x + m * x / (1 + h * m * x)) - delta * t'
    diffusion = 'sigma'
  }
  else if (diffusion.type == 'geometric') {
    diffusion = 'sigma * x'
  }
  else
    stop('Please assign a proper diffusion type: wiener process or geometric brownian motion.')
  mod.nonstructural = setModel(drift = drift, diffusion = diffusion) # , solve.variable = 'x, state.variable = 'x', time.variable = 't'
  # set sampling
  samples = setSampling(Terminal = steps * stepwise, n = steps)
  # initialize yuima object
  yuima.nonstructural <- setYuima(model = mod.nonstructural, sampling = samples)
  # assign values for parameters
  parameters = list(r = r, s = s, m = m, h = h, delta = delta, sigma = sigma)
  # simulate and update yuima object
  yuima.nonstructural = simulate(yuima.nonstructural, true.parameter = parameters, xinit = xinit)
  yuima.nonstructural
}



#####################################################################
# Functions of the non-structural model
# dx = x (r - s x - k_c c x + k_m m x / (1 + h k_m m x))

#' @title plot the feasible(x>0) and infeasible(x<0) equilibriums of the nonstructural model I
#'  control param is r, state variable is x.
#'  r_max = (\sqrt(m) - \sqrt(s))^2 / h * m
#'  X1_max = (\sqrt(2)-1)(\sqrt(s) - \sqrt(m)) / (\sqrt(s) * h * m)
#'  X2_max = -(\sqrt(2)+1)(\sqrt(s) - \sqrt(m)) / (\sqrt(s) * h * m)
#' @param m, total mutualistic strength
#' @param s, a vector of total competitive strengths, from 0.5*m to 2*m
#' @param h, handling time
#' @example plot_equilibrium_rx(m = 1, s = c(0.5, 1, 1.5) * m, h = 0.3)
plot_equilibrium_rx <- function(m, s, h) {
  r_max = (sqrt(m) - sqrt(2 * m))^2 / h * m
  tmp <- ldply(s, function(s) {
    A = s * h * m
    B = s - m - r_max * h * m
    C = - r_max
    delta = B^2 - 4 * A * C
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    X = seq(from = X2, to = X1, length.out = 500)
    r = s * X - m * X / (1 + h * m * X)
    data.frame(s = s, X = X, r = r) # s = s, X1 = X1, X2 = X2,
  } )
  #colnames(tmp) <- c('strong', 'equal', 'weak')
  ggplot(data = subset(tmp, X >= 0), aes(x = X, y = r, group = factor(s), color = factor(s))) +
    geom_line(size = 1.) +
    geom_line(data = subset(tmp, X <= 0), aes(x = X, y = r, group = factor(s), color = factor(s)), size = 1., linetype = 3) +
    geom_vline(aes(xintercept = 0), size = 0.3) +
    geom_hline(aes(yintercept = 0), size = 0.3) +
    scale_color_manual(values = gg_color_hue(3),
                       labels = c("Strong",
                                  "Equal",
                                  "Weak")) +
    guides(color = guide_legend(title = "Type")) +
    coord_flip() +
    theme_bw()
}

#' @example ggplot_fold_rx_h(m = 1.5, s = 1, h = c(0.5, 1, 1.5) * m)
ggplot_fold_rx_h <- function(m, s, h) {
  rmax = 0
  tmp <- ldply(h, function(h) {
    rmin = - (sqrt(m) - sqrt(s))^2 / (h * m)
    xmin = sqrt(s) * (sqrt(m) - sqrt(s)) / (s * h * m)
    xr0 <- -(s - m) / (s * h * m)
    r = seq(from = rmin, to = rmax, length.out = 2000)
    A = s * h * m
    B = s - m - r * h * m
    C = - r
    delta = B^2 - 4 * A * C
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    data.frame(h = h, xmin = xmin, rmin = rmin, xr0 = xr0, X1 = X1, X2 = X2, r = r) # s = s, X1 = X1, X2 = X2,
  } )
  tmp <- subset(tmp, !is.nan(X1))
  require(grid) # define arrow
  ggplot(data = tmp, aes(x = r, y = X1, group = factor(h), color = factor(h))) +
    geom_line(size = 1.) +
    geom_line(data = tmp, aes(x = r, y = X2, group = factor(h), color = factor(h)), size = 1., linetype = 3) +
    #    geom_vline(aes(xintercept = 0), size = 0.3) +
    geom_segment(aes(x = rmin, y = 0, xend = 0, yend = 0), color = 'black', size = 1) +
    geom_segment(aes(x = rmin, y = xmin, xend = rmin, yend = 0), color = 'black', size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = xr0), color = 'black', size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
    scale_color_manual(values = gg_color_hue(3),
                       labels = c("0.5*m", "1.0*m", "1.5*m")) +
    guides(color = guide_legend(title = "Handling\nTime")) +
    labs(y = 'X') +
    theme_bw()
}

#' @example ggplot_fold_rx_m(m = c(1.5, 2, 2.5) * s, s = 1, h = 0.3)
ggplot_fold_rx_m <- function(m, s, h) {
  rmax = 0
  tmp <- ldply(m, function(m) {
    rmin = - (sqrt(m) - sqrt(s))^2 / (h * m)
    xmin = sqrt(s) * (sqrt(m) - sqrt(s)) / (s * h * m)
    xr0 <- -(s - m) / (s * h * m)
    r = seq(from = rmin, to = rmax, length.out = 2000)
    A = s * h * m
    B = s - m - r * h * m
    C = - r
    delta = B^2 - 4 * A * C
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    data.frame(m = m, xmin = xmin, rmin = rmin, xr0 = xr0, X1 = X1, X2 = X2, r = r) # s = s, X1 = X1, X2 = X2,
  } )
  tmp <- subset(tmp, !is.nan(X1))
  require(grid) # define arrow
  ggplot(data = tmp, aes(x = r, y = X1, group = factor(m), color = factor(m))) +
    geom_line(size = 1.) +
    geom_line(data = tmp, aes(x = r, y = X2, group = factor(m), color = factor(m)), size = 1., linetype = 3) +
#    geom_vline(aes(xintercept = 0), size = 0.3) +
    geom_segment(aes(x = rmin, y = 0, xend = 0, yend = 0), color = 'black', size = 1) +
    geom_segment(aes(x = rmin, y = xmin, xend = rmin, yend = 0), color = 'black', size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = xr0), color = 'black', size = 0.3, arrow = arrow(length = unit(0.02, "npc"))) +
    scale_color_manual(values = gg_color_hue(3),
                       labels = c("1.5*s", "2*s", "2.5*s")) +
    guides(color = guide_legend(title = "Mutualism\nStrength")) +
    labs(y = 'X') +
    theme_bw()
}

#' @title plot a fold bifurcation, for the nonstructural model I,
#'  control param is r, state variable is x
#' @param m, total mutualistic strength
#' @param s, total competitive strength
#' @param h, handling time
#' @example plot_fold_rx(m = 1.5, s = 1, h = 0.3)
plot_fold_rx <- function(m, s, h) {
  rmin = round(- (sqrt(m) - sqrt(s))^2 / (h * m), 12)
  xmin = sqrt(s) * (sqrt(m) - sqrt(s)) / (s * h * m)
  xr0 <- -(s - m) / (s * h * m)
  rs = seq(from = rmin, to = 0, length.out = 200)
  require(plyr)
  tmp <- ldply(rs, function(r) {
    A = s * h * m
    B = s - m - r * h * m
    C = - r
    delta = B^2 - 4 * A * C
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    c(r = r, X1 = X1, X2 = X2)
  })
  matplot(tmp$r, tmp[,c('X1', 'X2')], type = 'l', lwd = 2, xlab = 'r', ylab = 'x', xlim = c(round(rmin, 2)- 0.01, 0+0.01), xaxp = c(round(rmin, 2) - 0.01, 0, 2), ylim = c(0, round(xr0, 2)+0.01), yaxp = c(0, round(xr0, 2)+0.01, 2))
  lines(c(rmin*1.05, 0), c(0, 0), lwd = 2)
  r1 <- tmp[nrow(tmp)/3,]$r
  x1r1 <- tmp[nrow(tmp)/3,]$X1
  x2r1 <- tmp[nrow(tmp)/3,]$X2
  r2 <- tmp[nrow(tmp)*2/3,]$r
  x1r2 <- tmp[nrow(tmp)*2/3,]$X1
  x2r2 <- tmp[nrow(tmp)*2/3,]$X2
  arrows(-0.0, xr0*0.2, -0.0, xr0*0.8, length = 0.1, lty = 1)
  arrows(r1, x2r1+(x1r1-x2r1)*0.1, r1, x2r1+(x1r1-x2r1)*0.9, length = 0.1, lty = 1)
  arrows(r2, x2r2+(x1r2-x2r2)*0.1, r2, x2r2+(x1r2-x2r2)*0.9, length = 0.1, lty = 1)
  arrows(rmin, xmin*0.8, rmin, xmin*0.2, length = 0.1, lty = 1)
  arrows(r1, x2r1*0.9, r1, x2r1*0.1, length = 0.1, lty = 1)
  arrows(r2, x2r2*0.9, r2, x2r2*0.1, length = 0.1, lty = 1)
  points(0, 0, pch = 16)
  text(0, 0, '(0,0)', pos = 3)
  points(0, xr0, pch = 16)
  text(0, xr0, expression(paste("(0,",x[0],")")), pos = 1)
  points(rmin, xmin, pch = 16)
  text(rmin, xmin, expression(paste("(", r[min], ',', x[min],")")), pos = 4)
}

#' @title plot a fold bifurcation, control param is m, state variable is x
#' @param r, s, h parameters
#' @example plot_fold_mx(r = -0.1, s = 1, h = 0.5)
plot_fold_mx <- function(r, s, h) {
  mmin = s * (sqrt(-4 * r * h) - r * h + 1) / (r * h + 1)^2
  mmin = mmin + 1e-12
  xmin = sqrt(s) * (sqrt(mmin) - sqrt(s)) / (s * h * mmin)
  #xr0 <- -(s - m) / (s * h * m)
  ms = seq(from = mmin, to = 2 * mmin, length.out = 200)
  require(plyr)
  tmp <- ldply(ms, function(m) {
    A = s * h * m
    B = s - m - r * h * m
    C = - r
    delta = B^2 - 4 * A * C
    X1 = (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
    X2 = (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
    c(m = m, X1 = X1, X2 = X2)
  })
  matplot(tmp$m, tmp[,c('X1', 'X2')], type = 'l', lwd = 2, xlab = 'm', ylab = 'x', xlim = c(round(mmin, 2)- 0.01, round(2 * mmin, 2)+0.01), xaxp = c(round(mmin, 2) - 0.01, round(2 * mmin, 2)+0.01, 2), ylim = c(0, round(max(tmp$X1), 2)+0.01), yaxp = c(0, round(max(tmp$X1), 2)+0.01, 2))
  lines(c(mmin, 0), c(2 * mmin, 0), lwd = 2)
}
