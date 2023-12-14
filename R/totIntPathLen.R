#' Calculation of the total internal path length for rooted trees
#'
#' This function calculates the total internal path length \eqn{TIP(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{TIP(T)} is defined as
#' \deqn{TIP(T)=\sum_{x\in V_{in}(T)} \delta(x)}{TIP(T)=\sum_{x in V_in(T)} depth(x)} in
#' which \eqn{V_{in}(T)}{V_in(T)} denotes the set of inner vertices of \eqn{T}, and \eqn{\delta(x)}{depth(x)}
#' denotes the depth of the vertex \eqn{x}. The total internal path length is an
#' imbalance index. \cr\cr
#' For details on the total internal path length, see 
#' also Chapter 23 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_23).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{totIntPathLen} returns the total internal path length of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references D. E. Knuth. The art of computer programming: fundamental algorithms, volume 1. Addison-Wesley, Reading, Mass, 3rd edition, 1997. ISBN 9780201896831.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,,),),(,)));")
#' totIntPathLen(tree)
#'
#'@export
totIntPathLen <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if (n==1) {
    return(0)
  } else{ # get the depth of each vertex in the tree
    allnodeDepths <- getNodesOfDepth(mat=getDescMatrix(tree), root=n+1, n=n)
    # summarize the depths of the inner vertices
    return(sum(allnodeDepths$nodeDepths[(n+1):length(allnodeDepths$nodeDepths)], na.rm=TRUE))
  }
}
