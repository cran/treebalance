#' Calculation of the area per pair index for rooted trees
#'
#' This function calculates the area per pair index \eqn{APP(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{APP(T)} is defined as
#' \deqn{APP(T)=\frac{2}{n\cdot(n-1)}\cdot\sum_{1\leq i<j\leq n} d_T(i,j)}{
#' APP(T)=2/(n(n-1))*\sum_{1<=i<j<=n} d_T(i,j)} in which \eqn{n} denotes the
#' number of leaves in \eqn{T}, and
#' \eqn{d_T(i,j)} denotes the number of edges on the path between the two
#' leaves \eqn{i} and \eqn{j}. Note that \eqn{APP(T)} can also be
#' computed from the Sackin index \eqn{S(T)} and the total cophenetic
#' index \eqn{TCI(T)} of \eqn{T} as
#' \eqn{APP(T)=\frac{2}{n}\cdot S(T)-\frac{4}{n(n-1)}\cdot TCI(T)}{APP(T)=2/n*S(T)-4/(n(n-1))*TCI(T)}
#' enabling efficient computation.\cr\cr
#' The area per pair index does not fulfill the definition of an (im)balance
#' index given in "Tree balance indices: a comprehensive survey" (Fischer et al., 2023). \cr\cr
#' For details on the area per pair index, see 
#' also Chapter 24 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_24).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{areaPerPairI} returns the area per pair index of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references T. Araujo Lima, F. M. D. Marquitti, and M. A. M. de Aguiar. Measuring Tree Balance with Normalized Tree Area. arXiv e-prints, art. arXiv:2008.12867, 2020.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' areaPerPairI(tree)
#'
#' @export
areaPerPairI <- function(tree)
{
  #Check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")

  n <- length(tree$tip.label)
  if(n == 1) return(0)

  APP <- 2/n*sackinI(tree) - 4/(n*(n-1))*totCophI(tree)
  return(APP)
}
