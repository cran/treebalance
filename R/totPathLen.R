#' Calculation of the total path length for rooted trees
#'
#' This function calculates the total path length \eqn{TPL(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{TPL(T)} is defined as
#' \deqn{TPL(T)=\sum_{x\in V(T)} \delta(x)}{depth(x)} in
#' which \eqn{V(T)} denotes the set of vertices of \eqn{T}, and \eqn{\delta(x)}{depth(x)}
#' denotes the depth of the vertex \eqn{x}. The total path length is an
#' imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{TPL(T)=0} and a warning. \cr\cr
#' For details on the total path length, see 
#' also Chapter 23 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_23).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{totPathLen} returns the total path length of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references see e.g. R. P. Dobrow, J. A. Fill. Total path length for random recursive trees. Combinatorics, Probability and Computing, 8(4):317–333, 1999. doi: 10.1017/S0963548399003855.
#' @references see e.g. L. Takacs. On the total heights of random rooted trees. Journal of Applied Probability, 29(3):543–556, 1992. doi: 10.2307/3214892.
#' @references see e.g. L. Takacs. On the total heights of random rooted binary trees. Journal of Combinatorial Theory, Series B, 61(2):155–166, 1994. ISSN 0095-8956. doi: 10.1006/jctb.1994.1041.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,,),),(,)));")
#' totPathLen(tree)
#'
#'@export
totPathLen <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if (n==1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  } else{ # get the depth of each vertex in the tree
    allnodeDepths <- getNodesOfDepth(mat=getDescMatrix(tree), root=n+1, n=n)
    # summarize the depths of all vertices
    return(sum(allnodeDepths$nodeDepths, na.rm=TRUE))
  }
}
