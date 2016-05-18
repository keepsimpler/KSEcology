#' @title estimate the largest eigenvalue of a bipartite competition-cooperation-mixture network
#'
estimate_second_eigenvalue_mixture <- function(graph) {
  # check if a regular graph
  stopifnot(length(unique(rowSums(graph))) == 1)
  # the dot eigenvalue
  lambda.dot <- unique(rowSums(graph))
  # the real largest eigenvalue
  eigenvalues <- eigen(graph)$values
  lambda1 <- eigenvalues[1]
  # if the real largest eigenvalue is the dot eigenvalue
  # then the second eigenvalue is the semicircle eigenvalue
  # if the real largest eigenvalue is larger than the dot
  # then the largest eigenvalue is the semicircle eigenvalue
  # the real largest eigenvalue is less than the dot, Impossible
  if (isTRUE(all.equal(lambda1, lambda.dot))) {
    lambda.semicircle <- eigenvalues[2]
  }
  else if (lambda1 > lambda.dot) {
    lambda.semicircle <- lambda1
  }
  # estimate the semicircle eigenvalue
  if (lambda1 > lambda.dot) {
    trace2 <- sum(diag(graph %*% graph))
    trace4 <- sum(diag(graph %*% graph %*% graph %*% graph))
    lambda.dot.abs <- unique(rowSums(abs(graph)))
    lambda.semicircle.trace2 <- 2 * sqrt((trace2 - lambda.dot.abs^2)/(n - 1))
    lambda.semicircle.trace4 <- 2 * sqrt(sqrt((trace4 - lambda.dot.abs^4)/(2 * (n - 1))))
  }
  else {
    trace2 <- sum(diag(graph %*% graph))
    trace4 <- sum(diag(graph %*% graph %*% graph %*% graph))
    lambda.dot.abs <- unique(rowSums(abs(graph)))
    lambda.semicircle.trace2 <- 2 * sqrt((trace2 - lambda.dot.abs^2 - lambda.dot^2)/(n - 2))
    lambda.semicircle.trace4 <- 2 * sqrt(sqrt((trace4 - lambda.dot.abs^4 - lambda.dot^4)/(2 * (n - 2))))
  }
  #competition_part <- graph
  #competition_part[competition_part > 0] <- 0
  #mutualism_part <- graph
  #mutualism_part[mutualism_part < 0] <- 0
  c(lambda1 = lambda1, lambda.dot = lambda.dot, lambda.semicircle = lambda.semicircle, lambda.semicircle.trace2 = lambda.semicircle.trace2, lambda.semicircle.trace4 = lambda.semicircle.trace4)
}
#' @title estimate the second largest eigenvalue of a regular bipartite graph
#' @param graph the adjacency matrix of a regular bipartite graph
#' @return the estimated second largest eigenvalue
estimate_second_eigenvalue_bipartite <- function(graph) {
  # check if a regular graph
  stopifnot(length(unique(rowSums(graph))) == 1)
  n = dim(graph)[1] # species number
  # for a bipartite regular graph,
  # the largest eigenvalue is always the row sum
  lambda1 <- unique(rowSums(graph))
  # estimate the second largest eigenvalue
  # semicircle plus twins method
  # compute trace(A^2) and trace(A^4)
  trace2 <- sum(diag(graph %*% graph))
  trace4 <- sum(diag(graph %*% graph %*% graph %*% graph))
  lambda2.trace2 <- 2 * sqrt((trace2 - 2 * lambda1^2)/(n - 2))
  lambda2.trace4 <- 2 * sqrt(sqrt((trace4 - 2 * lambda1^4)/(2 * (n - 2))))
  # alon's method
  lambda2.alon <- 2 * sqrt(lambda1 - 1)
  # a empirical method
  lambda2.empirical <- n/2 - lambda1
  # the real second eigenvalue
  lambda2 <- eigen(graph)$values[2]
  return(c(real = lambda2, alon = lambda2.alon,
           trace2 = lambda2.trace2,
           trace4 = lambda2.trace4,
           empirical = lambda2.empirical))
}

#' @title the second raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu2, the second raw moment \mu^2
Mu2 <- function(r,n) {
  mu2 =  r^2 * factorial(2/n) * gamma(3/n) /
    (factorial(1/n) * gamma(4/n))
  return(mu2)
}

#' @title the fourth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu4, the fourth raw moment \mu^4
Mu4 <- function(r,n) {
  mu4 = 8 * r^4 * 2 * gamma(2/n) * gamma(5/n) /
    (3 * gamma(1/n) * gamma(6/n))
  return(mu4)
}

#' @title the sixth raw moment of semisuperellipse distribution (r, n)
#' @param r  n,  the parameters of semi-superellipse distribution
#' @return mu6, the sixth raw moment \mu^6
Mu6 <- function(r,n) {
  mu4 = 8 * r^6 * 2 * gamma(2/n) * gamma(7/n) /
    ( gamma(1/n) * gamma(8/n))
  return(mu4)
}


## Copyed from Allesina's code
PDFSuperEllipseSymmetric <- function(x, r, n){
  LogDenominator <- log(2) + (1/n) * log((r)^n - abs(x)^n)
  LogNumerator <- log(4) + 2 * log(r) + 2 * lgamma(1 + 1 / n) - lgamma(1 + 2/n)
  res <- (exp(LogDenominator - LogNumerator))
  res[is.nan(res)] <- 0
  return(res)
}

FindSECounts <- function(Hbreaks, Hcounts, r, n){
  data.points <- length(Hcounts)
  for (i in 1:data.points){
    Hcounts[i] <- integrate(PDFSuperEllipseSymmetric, Hbreaks[i], Hbreaks[i+1], r, n)$value
  }
  return(Hcounts)
}

###############################################################################
#' @title  Example of the semi-circle plus twins law of bipartite REGULAR random graph
#' @param n, number of nodes
#' @param km, node degree
#' @param title of figure
#' @example
#' p1 <- semicircle.plus.twins.of.bipartite.regular(n = 1000, k = 10, title = 'n = 1000, k = 10')
###############################################################################
semicircle.plus.twins.of.bipartite.regular <- function(n, km, title) {
  #graph <- gen_competition_cooperation_bipartite(s1 = n1, s2 = n2, kc = kc, km = km, typec = typec, typem = typem)
  graph <- gen_regular_competition_cooperation(n = n, kc = 0, km = km, sc = -1, sm = 1)
  eigenvalues = eigen(graph)$values
  lev = eigenvalues[1]
  sev = eigenvalues[2]
  H <- hist(eigenvalues, breaks = 75, plot = FALSE)
  Hbreaks <- H$breaks
  Hmids <- H$mids
  Hcounts <- H$counts
  temp <- data.frame(
    x = Hmids,
    y = Hcounts,
    lev = max(Hmids[Hmids < lev]), # the closest Hmids(hist midpoints) that less than lev, in order to display more concisely.
    sev = min(Hmids[Hmids > sev])
  )
  y = FindSECounts(Hbreaks, Hbreaks[-1], r = sev, n = 2) * 2 * (n - 2)
  temp2 <- data.frame(
    x = Hmids,
    y = y
  )
  require(ggplot2)
  p <- ggplot(data = temp, aes(x = x, ymax = y, ymin = 0)) +  # p1 p2 p3
    geom_linerange(size = 0.75) +
    geom_line(data = temp2, aes(x = x, y = y), colour = "red") +
    geom_vline(aes(xintercept = lev), linetype = 'dashed', size = 0.75, color = 'red')  +
#    geom_vline(aes(xintercept = - lev), linetype = 'dashed', size = 0.75, color = 'red')  +
    geom_vline(aes(xintercept = sev), linetype = 'dashed', size = 0.75, color = 'green')  +
    scale_x_continuous(expression(lambda)) +
    scale_y_continuous('Density(Counts)') +
#    geom_text(label = expression(lambda[1]), x = sev, y = 5, colour = "red") +
    ggtitle(title) +
    theme_bw()
  p
}


semicircle.plus.one.of.mixture.regular <- function(n, km, title, kc = NULL) {
  if (is.null(kc))
    kc = n - 1
  #graph <- gen_competition_cooperation_bipartite(s1 = n1, s2 = n2, kc = kc, km = km, typec = typec, typem = typem)
  graph <- gen_regular_competition_cooperation(n = n, kc = kc, km = km, sc = -1, sm = 1)
  eigenvalues = eigen(graph)$values
  lev = unique(rowSums(graph))
  sev = eigenvalues[1]
  H <- hist(eigenvalues, breaks = 75, plot = FALSE)
  Hbreaks <- H$breaks
  Hmids <- H$mids
  Hcounts <- H$counts
  temp <- data.frame(
    x = Hmids,
    y = Hcounts,
    lev = max(Hmids[Hmids < lev]), # the closest Hmids(hist midpoints) that less than lev, in order to display more concisely.
    sev = max(Hmids[Hmids < sev])
  )
  y = FindSECounts(Hbreaks, Hbreaks[-1], r = sev, n = 2) * 2 * (n - 1)
  temp2 <- data.frame(
    x = Hmids,
    y = y
  )
  require(ggplot2)
  p <- ggplot(data = temp, aes(x = x, ymax = y, ymin = 0)) +  # p1 p2 p3
    geom_linerange(size = 0.75) +
    geom_line(data = temp2, aes(x = x, y = y), colour = "red") +
    geom_vline(aes(xintercept = lev), linetype = 'dashed', size = 0.75, color = 'red')  +
    #    geom_vline(aes(xintercept = - lev), linetype = 'dashed', size = 0.75, color = 'red')  +
    geom_vline(aes(xintercept = sev), linetype = 'dashed', size = 0.75, color = 'green')  +
    scale_x_continuous(expression(lambda)) +
    scale_y_continuous('Density(Counts)') +
    #    geom_text(label = expression(lambda[1]), x = sev, y = 5, colour = "red") +
    ggtitle(title) +
    theme_bw()
  p
}
