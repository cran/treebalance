#' Calculation of the average vertex depth for rooted trees
#'
#' This function calculates the average vertex depth \eqn{AVD(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{AVD(T)} is defined as
#' \deqn{AVD(T)=\frac{1}{|V(T)|}\cdot\sum_{x\in V(T)} \delta(x)}{AVD(T)=1/|V(T)|*\sum_{x\in V(T)} \delta(x)} in
#' which \eqn{V(T)} denotes the set of vertices of \eqn{T}, and \eqn{\delta(x)}{depth(x)}
#' denotes the depth of the vertex \eqn{x}. The average vertex depth is an
#' imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{AVD(T)=0} and a warning.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{avgVertDep} returns the average vertex depth of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references A. Herrada et al. Scaling properties of protein family phylogenies. BMC Evolutionary Biology, 11(1), June 2011. doi: 10.1186/1471-2148-11-155.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,,),),(,)));")
#' avgVertDep(tree)
#'
#'@export
avgVertDep <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if (n==1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  } else{ # get the depth of each vertex in the tree
    allnodeDepths <- getNodesOfDepth(mat=getDescMatrix(tree), root=n+1, n=n)
    # summarize the depths of all vertices
    return(sum(allnodeDepths$nodeDepths, na.rm=TRUE)/(n+tree$Nnode))
  }
}
