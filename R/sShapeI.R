#' Calculation of the s-shape statistic for rooted trees
#'
#' This function calculates the s-shape statistic \eqn{sShape(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary, however \eqn{sShape} only fulfils
#' the definition of an imbalance index on the space of binary trees. \eqn{sShape(T)} is defined as
#' \deqn{sShape(T)=\sum_{u\in V_{in}(T)} log(n_u-1)}{sShape(T)=\sum_{u in V_in(T)} log(n_u-1)} in
#' which \eqn{V_{in}(T)}{V_in(T)} denotes the set of inner vertices of \eqn{T}
#' and \eqn{n_u} denotes the number of leaves
#' in the pending subtree that is rooted at \eqn{u}. An arbitrary logarithm base can be used
#' (for binary trees it is common to use base 2).\cr\cr
#' For \eqn{n=1} the function returns \eqn{sShape(T)=0} and a warning.
#'
#' @param tree A rooted tree in phylo format.
#' @param logbase The logarithm base that shall be used.
#'
#' @return \code{sShapeI} returns the s-shape statistic of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references M.G. Blum and O. Francois. Which random processes describe the tree of life? a large-scale study of phylogenetic tree imbalance. Systematic Biology, 2006.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' sShapeI(tree)
#'
#' @export
sShapeI <- function(tree, logbase=2){
  #check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n == 1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }
  numdescleaves <- get.subtreesize(tree)[(n+1):(n+tree$Nnode)]
  s_val <- sum(sapply(numdescleaves, function(x) log(x-1, base=logbase)))
  return(s_val)
}
