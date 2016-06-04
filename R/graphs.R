##### generate graphs needed by the models --------------------------

#' @title generate the adjacency matrix of a kind of competition-cooperation-mixed network. The competition part and cooperation part are mixed with each other.
#' @description Unlike \code{gen_competition_cooperation_bipartite}, the competition part and cooperation part can not be generated separately. Instead, a random graph is generated first, then randomly assign competition(--) or cooperation(++) to each interaction according to their proportions.
#' @param n species number
#' @param k average degree(number of interactions) of species
#' @param type network types: 'gnm', 'sf', ...
#' @param pc proportion of competition interactions
gen_competition_cooperation_unipartite <- function(n, k, type, pc, ...) {
  stopifnot(pc >= 0., pc <= 1.)
  G = gen_connected_graph(n, k, type, ...)  # generate a connected graph
  graph = as.matrix(get.adjacency(G))  # transform to matrix form
  # split the graph to two sub-graphs - competition and mutualism graphs according to the probability of occurance of two different types of interactions
  # adjacency matrix is symmetric, first erase the lower-tringle
  graph[lower.tri(graph)] <- 0
  ps = runif(sum(graph))
  graph_mark <- graph
  graph_mark[graph_mark > 0] <- ps
  graph[graph_mark > 0 & graph_mark <= pc] <- -1 # competition
  return(graph + t(graph))
}

#' @title generate the adjacency matrix of a kind of competition-cooperation-mixed network. The competition part and cooperation part are separated
#' @description The adjacency matrix should include two parts: the competition part and the cooperation part. The competition part is marked by Negative values(-1), the cooperation part is marked by Positive values(1).
#' The network is constructed through such a schema: two groups of nodes, nodes in the same group have competitive interactions, nodes between groups have cooperative interactions.
#' @param s1, s2 number of nodes of two groups
#' @param kc, km average degree of nodes of two groups
#' @param typec, graph type of the competition part, 'full', 'regular', 'gnm'
#' @param typem, graph type of the cooperation part, 'bipartite-regular', 'bipartite-regular2', 'bipartite-gnm', 'bipartite-nest?'
gen_competition_cooperation_bipartite <- function(s1, s2, kc, km, typec = 'full', typem = 'bipartite-regular') {
  # ensure that the cooperation part is a bipartite graph
  stopifnot(substr(typem, start = 1, stop = 9) == 'bipartite')
  # the mutualistic part
  graphm <- gen_adjacency_matrix(c(s1, s2), km, typem)
  # the competitive part
  graphc = matrix(0, nrow = s1 + s2, ncol = s1 + s2)
  if (typec == 'full') {
    graphc[1:s1, 1:s1] = as.matrix(as_adjacency_matrix(make_full_graph(s1))) * -1
    graphc[(s1+1):(s1+s2), (s1+1):(s1+s2)] = as.matrix(as_adjacency_matrix(make_full_graph(s2))) * -1
  }
  else if (typec == 'regular') {
    graphc[1:s1, 1:s1] = as.matrix(as_adjacency_matrix(sample_k_regular(s1, kc))) * -1
    graphc[(s1+1):(s1+s2), (s1+1):(s1+s2)] = as.matrix(as_adjacency_matrix(sample_k_regular(s2, kc))) * -1
  }
  else if (typec == 'gnm') {
    graphc[1:s1, 1:s1] = as.matrix(as_adjacency_matrix(sample_gnm(s1, kc * s1 / 2))) * -1
    graphc[(s1+1):(s1+s2), (s1+1):(s1+s2)] = as.matrix(as_adjacency_matrix(sample_gnm(s2, kc * s2 / 2))) * -1
  }
  graph = graphc + graphm
  return(graph)
}


#' @title generate regular competition-cooperation-mixed network
#' @description The adjacency matrix of network should include two parts: the competition part and the cooperation part. The network is constructed through such a schema: two groups of nodes with same node number, nodes in the same group have competitive interactions, nodes between groups have cooperative interactions.
#' @param n number of nodes in each group
#' @param kc competition degree
#' @param km cooperation degree
#' @param sc competition strength
#' @param sm cooperation strength
#' @return the adjacency matrix
gen_regular_competition_cooperation <- function(n, kc, km, sc, sm) {
  stopifnot(sc <= 0 & sm >= 0 & kc <= n - 1 & km <= n - 1 )
  # initialize a empty graph with 2 group nodes
  G = matrix(0, nrow = 2 * n, ncol = 2 * n)
  # generate the competition part of regular graph
  G[1:n, 1:n] = as.matrix(as_adjacency_matrix(sample_k_regular(n, kc))) * sc
  G[(n+1):(2*n), (n+1):(2*n)] = as.matrix(as_adjacency_matrix(sample_k_regular(n, kc))) * sc
  # generate the cooperation part of regular graph
  G[1:n, (n+1):(2*n)] = as.matrix(as_adjacency_matrix(sample_k_regular(n, km))) * sm
  G[(n+1):(2*n), 1:n] = t(G[1:n, (n+1):(2*n)])
  G
}


# #' @title Graph generators
# #' @param type graph type: 'bipartite-regular', 'bipartite-gnm', 'bipartite-configuration-model'
# #' @param coeffs, coefficients, a list including s1, s2, k
# #' @return the adjacency matrix of a graph
# gen_graph <- function(type, coeffs) {
#   s1 = coeffs$s1
#   s2 = coeffs$s2
#   k = coeffs$k
#   if (type == 'bipartite-regular') {
#     if (s1 != s2)
#       stop('s1 == s2, if you want a regular bipartite graph')
#     G = sample_k_regular(s1, k)
#     graph <- as.matrix(as_adjacency_matrix(G))
#     graph <- inc_to_adj(graph)  # here a trick
#   }
#   else if (type == 'bipartite-gnm') {
#     G = sample_bipartite(n1 = s1, n2 = s2, type = 'gnm', m = k * (s1 + s2))
#     graph <- as.matrix(as_adjacency_matrix(G))
#   }
#   else if (type =='regular') {
#     G = sample_k_regular(no.of.nodes = s1, k = k)
#     graph <- as.matrix(as_adjacency_matrix(G))
#   }
#   graph
# }


#' @title generate adjacency matrix of different types of graphs
#' @param s, node number. For bipartite graphs, [s] is a length 2 vector.
#' @param k, average node degree.
#' @param type, graph types:
#' Bipartite graphs: [s] is a length 2 vector (s1,s2)
#'    - bipartite-regular: bipartite regular graph that each node has the same degree, which impliciply require two group are equal, s1 = s2
#'    - bipartite-regular2: bipartite regular graph that nodes in each separated group have their same degrees. In this situation, [k] is the node degree in group 1, the node degree in group 2 is calculated by s1 * k / s2.
#'    - bipartite-gnm: bipartite random graph
#' Unipartite graphs: [s] is a integer.
#'    -
gen_adjacency_matrix <- function(s, k, type, ...) {
  stopifnot((type == 'bipartite-regular' && length(s) == 2 && s[1] == s[2]) ||
              (type == 'bipartite-regular2' && length(s) == 2) ||
              (type == 'bipartite-gnm' && length(s) == 2) ||
              (type != 'bipartite-gnm' && type != 'bipartite-regular' && type != 'bipartite-regular2'))
  if (type == 'bipartite-regular') {
    G = sample_k_regular(s[1], k)
    graph <- as.matrix(as_adjacency_matrix(G))
    graph <- inc_to_adj(graph)  # here a trick
  }
  else if (type == 'bipartite-regular2') {
    graph <- gen_bipartite_regular(s, k)
  }
  else {
    G = gen_connected_graph(s, k, type, ...)
    graph <- as.matrix(as_adjacency_matrix(G))
  }
  return(graph)
}

#' @title generate a bipartite regular graph
#' @description Bipartite regular graph can be distinguished to two types. If species number in two grups are same, all species have the same number of interactions. If species number in two groups are different, species in group 1 have the same number of interactions \code{k}, species in group 2 have the same number of interactions \code{s1 * k / s2}.
#' The generation of the second type bipartite regular graph need to use \code{Bipartite configuration model} implemented in \code{NetworkX} package written by \code{Python}, thus we need \code{rPython} package to call \code{Bipartite configuration model} in \code{NetworkX}
#' @import rPython
gen_bipartite_regular <- function(s, k, maxtried = 100) {
  require(rPython)
  stopifnot(length(s) == 2 && s[1] <= s[2])
  if (s[1] == s[2]) {
    G = sample_k_regular(s[1], k)
    graph <- as.matrix(as_adjacency_matrix(G))
    graph <- inc_to_adj(graph)  # here a trick
  }
  else {
    # ensure that s[1] * k is a integer multiplication of s[2]
    # that is required to keep regular degree of s[2]
    stopifnot((s[1] * k) %% s[2] == 0 # remainder
              && (s[1] * k) %/% s[2] >= 1) # mod, 1 is the least degree of group 2
    k1 = k  # regular degree of group 1
    k2 = (s[1] * k) %/% s[2] #regular degree of group 2
    aseq = rep(k1, s[1])
    bseq = rep(k2, s[2])

    python.exec('import networkx as nx')
    python.exec('import networkx.algorithms.bipartite as bipartite')
    python.assign( "aseq",  aseq )
    python.assign( "bseq",  bseq )

    count = 0
    repeat{
      print(count)
      python.exec('G = bipartite.configuration_model(aseq, bseq, create_using=nx.Graph())')
      python.exec('graph = nx.adj_matrix(G)')
      python.exec('graph = graph.toarray().tolist()')
      graph <- python.get('graph')
      graph <- do.call(cbind, graph)
      flag = all(rowSums(graph)[1:s[1]] == aseq) && all(rowSums(graph)[(s[1]+1):(s[1]+s[2])] == bseq)
      if (flag) break  # until a regular bipartite graph is generated
      count = count + 1
      if (count == maxtried) {
        warning(paste('Tried', maxtried, 'times, But regular bipartite graph still cannot be generated.'))
        break
      }
    }
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
#' @param gtype, Graph type generated: 'bipartite-gnm', 'sf', 'er', 'dag', 'regular'.
#' @param maxtried, the maximum number of tried times.
#' If have tried [maxtried] times, the function will return no matter whether the connected graph is generated.
#' @param expower exponent coefficient of Scale-Free network
#' @param ... the params conform to the implementation functions of [igraph]
#' @return the connected graph
#' @details .
#' @import igraph
gen_connected_graph <- function(s, k, gtype, maxtried = 100, expower = 2.5, ...) {
  #library(igraph)
  if (gtype == 'bipartite-gnm' && is.na(s[2])) {  # the bipartite graph need size of two groups of nodes
    warning('sizes of TWO groups of nodes should be designated.
            we have assumed the size of second group equal to the size of first group.')
    s[2] = s[1]  # if missed second size, we assume it equal to the first size.
  }
  count = 0
  repeat {  # generate a connected graph
    if (gtype == 'bipartite-gnm') {
      G = sample_bipartite(s[1], s[2], type = 'gnm', m = ceiling(k * (s[1] + s[2])  / 2))
    } else if (gtype == 'sf') {
      G = sample_fitness_pl(s, k * s / 2, exponent.out = expower)
    }
    else if (gtype == 'er') {
      G = sample_gnm(s, k * s / 2)
    }
    else if (gtype == 'regular') {
      G = sample_k_regular(s, k)
    }
    else if (gtype == 'complete') {
      G = make_full_graph(s)
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
