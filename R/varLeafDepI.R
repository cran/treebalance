#' Calculation of the variance of leaf depths index for rooted trees
#'
#' This function calculates the variance of leaf depths index \eqn{VLD(T)}
#' for a given rooted tree \eqn{T}. The tree must not necessarily be binary.
#' \eqn{VLD(T)} is defined as
#' \deqn{VLD(T)=\frac{1}{n}\cdot\sum_{x\in V_L(T)} (\delta(x)-N(T))^2}{
#' VLD(T)=1/n * \sum_{x in V_L(T)} (\delta(x)-N(T))^2}
#' in which \eqn{n} denotes the number of leaves of \eqn{T}, \eqn{V_L(T)}{V_L(T)}
#' denotes the set of leaves of \eqn{T}, \eqn{\delta(x)} denotes the depth of
#' the leaf \eqn{x} and \eqn{N(T)} denotes the average leaf depth of \eqn{T}.\cr\cr
#' For \eqn{n=1} the function returns \eqn{VLD(T)=0} and a warning. \cr\cr
#' For details on the variance of leaf depths, see 
#' also Chapter 7 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_7).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{varLeafDepI} returns the variance of leaf depths index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references T. M. Coronado, A. Mir, F. Rossello, and L. Rotger. On Sackin's original proposal: the variance of the leaves' depths as a phylogenetic balance index. BMC Bioinformatics, 21(1), 2020. doi: 10.1186/s12859-020-3405-1. URL https://doi.org/10.1186/s12859-020-3405-1.
#' @references M. J. Sackin. "Good" and "Bad" Phenograms. Systematic Biology, 21(2):225-226, 1972. doi: 10.1093/sysbio/21.2.225.
#' @references K.-T. Shao and R. R. Sokal. Tree Balance. Systematic Zoology, 39(3):266, 1990. doi: 10.2307/2992186.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' varLeafDepI(tree)
#'
#' @export
varLeafDepI <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")

  n <- length(tree$tip.label)
  if(n == 1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  # calculate variance of leaf depths
  Descs <- getDescMatrix(tree)
  nodedepths <- getNodesOfDepth(mat=Descs,root=n+1,n=n)$nodeDepths
  avgLD <- avgLeafDepI(tree)
  return(1/n*sum((nodedepths[1:n]-avgLD)^2))
}
