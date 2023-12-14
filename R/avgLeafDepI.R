#' Calculation of the average leaf depth index for rooted trees
#'
#' This function calculates the average leaf depth \eqn{N(T)} for a given
#' rooted tree \eqn{T}. The tree must not necessarily be binary. \eqn{N(T)} is
#' defined as \deqn{N(T)=\frac{1}{n}\cdot\sum_{u\in V_{in}(T)} n_u}{
#' N(T)=(1/n) * \sum_{u in V_in(T)} n_u} in which \eqn{n} denotes the number of leaves in \eqn{T},
#' \eqn{V_{in}(T)}{V_in(T)} denotes the set of inner nodes of \eqn{T} and
#' \eqn{n_u} denotes the number of leaves in the pending subtree that is rooted
#' at the inner node \eqn{u}. Note that \eqn{N(T)} can also be
#' computed from the Sackin index \eqn{S(T)} as \eqn{N(T)=\frac{1}{n}\cdot S(T)}{N(T)=1/n*S(T)}.
#' The average leaf depth is an imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{N(T)=0} and a warning. \cr\cr
#' For details on the average leaf depth, see also Chapter 6 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_6).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{avgLeafDepI} returns the average leaf depth of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references M. J. Sackin. "Good" and "Bad" Phenograms. Systematic Biology, 21(2):225-226, 1972. doi: 10.1093/sysbio/21.2.225.
#' @references K.-T. Shao and R. R. Sokal. Tree Balance. Systematic Zoology, 39(3):266, 1990. \cr doi: 10.2307/2992186.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' avgLeafDepI(tree)
#'
#'@export
avgLeafDepI <- function(tree)
{
  #check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")

  n <- length(tree$tip.label)

  if(n == 1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  return (sackinI(tree)/n)
}
