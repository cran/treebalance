#' Calculation of the rooted quartet index for rooted trees
#'
#' This function calculates the rooted quartet index \eqn{rQI(T)} for a given
#' rooted tree \eqn{T}. The tree must not necessarily be binary.\cr\cr
#' Let \eqn{T} be a rooted tree, whose leaves are \eqn{1,...,n}. Let
#' \eqn{P_4} denote the set of all subsets of \eqn{\{1,...,n\}} that have
#' cardinality 4. Let \eqn{T(Q)} denote the rooted quartet on \eqn{Q\in P_4}
#' that is obtained by taking the subgraph of \eqn{T} that is induced by
#' \eqn{Q} and supressing its outdegree-1 vertices. \eqn{T(Q)} can have one of
#' the five following shapes:\cr\cr
#'     - \eqn{Q_0^*}: This is the caterpillar tree shape on 4 leaves, i.e.
#'     \code{"(,(,(,)));"} in Newick format. It has 2 automorphisms.\cr
#'     - \eqn{Q_1^*}: This is the tree shape on 4 leaves that has three pending
#'     subtrees rooted at the children of the root of \eqn{T}, one of them being
#'     a cherry and the other two
#'     being single vertices, i.e. \code{"((,),,);"} in Newick format. It has 4
#'     automorphisms.\cr
#'     - \eqn{Q_2^*}: This is the tree shape on 4 leaves that has two pending
#'     subtrees rooted at the children of the root of \eqn{T}, one of them being
#'     a star tree shape on 3 leaves
#'     and the other one being a single vertex, i.e. \code{"((,,),);"} in Newick
#'     format. It has 6 automorphisms.\cr
#'     - \eqn{Q_3^*}: This is the fully balanced binary tree shape on 4 leaves,
#'     i.e. \code{"((,),(,));"} in Newick format. Its has 8 automorphisms.\cr
#'     - \eqn{Q_4^*}: This is the star tree shape on 4 leaves, i.e.
#'     \code{"(,,,);"} in Newick format. It has 24 automorphisms.\cr\cr
#' \eqn{T(Q)} is assigned an rQI-value based on its shape, i.e. \eqn{rQI(T(Q))=q_i}
#' if \eqn{T(Q)} has the shape \eqn{Q_i^*}. The values \eqn{q_0,...,q_4} are
#' chosen in such a way that they increase with the symmetry of the shape as
#' measured by means of its number of automorphisms. Coronado et al. (2019)
#' suggested the values \eqn{q_0=0} and \eqn{q_i=i} or \eqn{q_i=2^i} for \eqn{i=1,...,4}.\cr
#' The rooted quartet index \eqn{rQI(T)} of the tree \eqn{T} is then defined as
#' the sum of the rQI-values of its rooted quartets:
#' \deqn{rQI(T)=\sum_{Q\in P_4} rQI(T(Q))}{rQI(T)=\sum rQI(T(Q)) over Q in P_4}
#' The rooted quartet index is a balance index.
#'
#' @param tree A rooted tree in phylo format.
#' @param shapeVal A vector of length 5 containing the shape values \eqn{q_0,...,q_4}.
#' Default is \eqn{(q_0,q_1,q_2,q_3,q_4)=(0,1,2,3,4)}.
#'
#' @return \code{rQuartetI} returns the rooted quartet index of the given tree based on the chosen shape values (see description for details).
#'
#' @author Sophie Kersting
#'
#' @references T. M. Coronado, A. Mir, F. Rosselló, and G. Valiente.  A balance index for phylogenetic trees based on rooted quartets. Journal of Mathematical Biology, 79(3):1105-1148, 2019. doi: 10.1007/s00285-019-01377-w. URL https://doi.org/10.1007/s00285-019-01377-w.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' rQuartetI(tree)
#'
#'@export
rQuartetI <- function(tree, shapeVal=c(0,1,2,3,4)){
  #Check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")

  n <- length(tree$tip.label)
  if(n<4) return(0)

  q0 <- shapeVal[1]
  q1 <- shapeVal[2]
  q2 <- shapeVal[3]
  q3 <- shapeVal[4]
  q4 <- shapeVal[5]

  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs, root=n+1, n=n)
  nodes_botup <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  Ancs <- getAncVec(tree)

  n_v <- rep(NA,n+tree$Nnode)
  ypsilon_v <- rep(NA,n+tree$Nnode)
  rQI_v <- rep(NA,n+tree$Nnode)

  for(v in nodes_botup){
    if(v<=n){ #if v is leaf
      n_v[v] <- 1
      ypsilon_v[v] <- 0
      rQI_v[v] <- 0
    }else{ #v is inner node
      n_v[v] <- sum(n_v[stats::na.omit(Descs[v,])])

      E_3_v <- auxE_l_X(3,n_v[stats::na.omit(Descs[v,])])
      E_4_v <- auxE_l_X(4,n_v[stats::na.omit(Descs[v,])])

      ypsilon_v[v] <- sum(ypsilon_v[stats::na.omit(Descs[v,])])+E_3_v

      rQI_v[v] <- sum(rQI_v[stats::na.omit(Descs[v,])])+
        q4 * E_4_v +
        q3 * auxE_l_X(2,choose(n_v[stats::na.omit(Descs[v,])],2))+
        q2 * ( n_v[v]*(ypsilon_v[v]-E_3_v)-
                 sum(n_v[stats::na.omit(Descs[v,])]*
                       ypsilon_v[stats::na.omit(Descs[v,])]) ) +
        q1 * ( 1/2* E_3_v * auxE_l_X(1,n_v[stats::na.omit(Descs[v,])]) -2*E_4_v-
                 3/2 * E_3_v)
    }
  }

  return(rQI_v[n+1])
}
