#' Calculation of the Sackin index for rooted trees
#'
#' This function calculates the Sackin index \eqn{S(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{S(T)} is defined as
#' \deqn{S(T)=\sum_{x\in V_L(T)} \delta(x)=\sum_{u\in V_{in}(T)} n_u}{S(T)=\sum_{x in V_L(T)} depth(x)=\sum_{u in V_in(T)} (n_u)} in
#' which \eqn{V_L(T)} denotes the leaf set of \eqn{T}, \eqn{\delta(x)}{depth(x)}
#' denotes the depth of the leaf \eqn{x}, \eqn{V_{in}(T)}{V_in(T)} denotes the set of
#' inner vertices in \eqn{T} and \eqn{n_u} denotes the number of leaves
#' in the pending subtree that is rooted at \eqn{u}. The Sackin index is an
#' imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{S(T)=0} and a warning.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{sackinI} returns the Sackin index of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references M.J. Sackin. "Good" and "Bad" Phenograms. Systematic Biology, 21(2):225-226, 1972. doi: 10.1093/sysbio/21.2.225.
#' @references K.-T. Shao and R.R. Sokal. Tree Balance. Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' sackinI(tree)
#'
#'@export
sackinI <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n == 1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }
  return(sum(get.subtreesize(tree)) - n)
}
