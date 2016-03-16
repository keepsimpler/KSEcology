#' @title generate regular competition-cooperation-mixed network
#' @description The adjacency matrix of network should include two parts: the competition part and the cooperation part. The network is constructed through such a schema: two groups of nodes with same node number, nodes in the same group have competitive interactions, nodes between groups have cooperative interactions.
#' @param n number of nodes in each group
#' @param kc competition degree
#' @param km cooperation degree
#' @param sc competition strength
#' @param sm cooperation strength
#' @return the adjacency matrix
gen_regular_competition_cooperation <- function(n, kc, km, sc, sm) {
  stopifnot(sc < 0 & sm > 0 & kc + km <= (n-1))
  G = matrix(0, nrow = n, ncol = n)
  for ( i in 1:n) {
      kc.count <- sum(G[i, ] < 0)
      km.count <- sum(G[i, ] > 0)
      while(kc.count < kc) {
        j <- sample(1:n, 1)
        if (j != i & G[i, j] == 0 & G[j, i] ==0) {
          G[i, j] = sc
          G[j, i] = sc
          kc.count = kc.count + 1
        }
      }
      while(km.count < km) {
        j <- sample(1:n, 1)
        if (j != i & G[i, j] == 0 & G[j, i] ==0) {
          G[i, j] = sm
          G[j, i] = sm
          km.count = km.count + 1
        }
      }
  }
  G
}

#' @title Graph generators
#' @param type graph type: 'bipartite-regular', 'bipartite-gnm', 'bipartite-configuration-model'
#' @param coeffs, coefficients, a list including s1, s2, k
#' @return the adjacency matrix of a graph
gen_graph <- function(type, coeffs) {
  s1 = coeffs$s1
  s2 = coeffs$s2
  k = coeffs$k
  if (type == 'bipartite-regular') {
    if (s1 != s2)
      stop('s1 == s2, if you want a regular bipartite graph')
    G = sample_k_regular(s1, k)
    graph <- as.matrix(as_adjacency_matrix(G))
    graph <- inc_to_adj(graph)  # here a trick
  }
  else if (type == 'bipartite-gnm') {
    G = sample_bipartite(n1 = s1, n2 = s2, type = 'gnm', m = k * (s1 + s2))
    graph <- as.matrix(as_adjacency_matrix(G))
  }
  else if (type =='regular') {
    G = sample_k_regular(no.of.nodes = s1, k = k)
    graph <- as.matrix(as_adjacency_matrix(G))
  }
  graph
}


###############################################################################
#' @title Generate a connected graph using package [igraph]
#'
#' @param s, size of network.
#' if graph type is bipartite, s[1], s[2] represent size of two groups; else s is size of network
#' @param k, average degree for the network.
#' 1 < k < s for unipartite network, 1 < k < s[1]*s[2]/(s[1]+s[2]) for bipartite network.
#' @param gtype, Graph type generated: 'bipartite', 'sf', 'er', 'dag', 'regular'.
#' @param maxtried, the maximum number of tried times.
#' If have tried [maxtried] times, the function will return no matter whether the connected graph is generated.
#' @param expower exponent coefficient of Scale-Free network
#' @param ... the params conform to the implementation functions of [igraph]
#' @return the connected graph
#' @details .
#' @import igraph
gen_connected_graph <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
  #library(igraph)
  if (gtype == 'bipartite' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated.
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  count = 0
  repeat {  # generate a connected graph
    if (gtype == 'bipartite') {
      G = bipartite.random.game(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2])))
    } else if (gtype == 'sf') {
      G = static.power.law.game(s, k * s, exponent.out = expower)
    }
    else if (gtype == 'er') {
      G = sample_gnm(s, k * s)
    }
    else if (gtype == 'regular') {
      G = k.regular.game(s, k)
    }
    else if (gtype == 'complete') {
      G = graph.full(s)
    }
    else if (gtype == 'dag') {
      require('spacejam')  # generate random directed Acyclic graphs
      G = rdag(s, s * k * 2)
    }
    if (igraph::is.connected(G)) break  # until a connected graph is generated
    count = count + 1
    if (count == maxtried) {
      warning(paste('Tried', maxtried, 'times, But connected graph still cannot be generated.'))
      break
    }
  }
  G
}
#plot(G, layout = layout.bipartite)

#' @title another form of uniform distribution between [mean - sd, mean + sd]
runif2 <- function(n, mean, sd) {
  runif(n) * 2 * sd + (mean - sd)
}

#' @title transfer an incidence matrix to an adjacency matrix
#' @param inc, an incidence matrix
#' @return adj, an adiacency matrix
inc_to_adj <- function(inc){
  p <- dim(inc)[1]  # number of Plants
  a <- dim(inc)[2]  # number of Animals
  s <- p + a  # number of all Species
  adj <- matrix(0, s, s)  # initialize the adjacency matrix as a zero-matrix
  adj[1:p, (p + 1):s] <- inc  # the upper right sub-matrix is the incidence matrix
  adj <- adj + t(adj)  # the lower left sub-matrix is transpose of the incidence matrix
  return(adj)
}

#' @title Bipartite Configuration Model
#' @description
