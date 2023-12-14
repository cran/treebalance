#' Calculation of the maximum width over maximum depth of the tree
#'
#' This function calculates the maximum width over maximum depth \eqn{mWovermD(T)} for a
#' given rooted tree \eqn{T}. The tree must not necessarily be binary. For \eqn{n>1},
#' \eqn{mWovermD(T)} is defined as \deqn{mWovermD(T)=maxWidth(T) / h(T)}
#' in which \eqn{h(T)} denotes the height of the tree \eqn{T}, which is the same as the 
#' maximum depth of any leaf in the tree, and \eqn{maxWidth(T)} denotes
#' the maximum width of the tree \eqn{T}. The maximum width over maximum depth
#' is a balance index. \cr\cr
#' For details on the maximum width over maximum depth, see 
#' also Chapter 23 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_23).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{mWovermD} returns the maximum width over maximum depth of a tree.
#'
#' @author Luise Kuehn
#' 
#' @references C. Colĳn and J. Gardy. Phylogenetic tree shapes resolve disease transmission patterns. Evolution, Medicine, and Public Health, 2014(1):96–108, 2014. doi: 10.1093/emph/eou018.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' mWovermD(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' mWovermD(tree)
#'
#' @export
mWovermD <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (length(tree$tip.label) == 1) {
    stop("mWovermD cannot be computed for n=1.")
  } else {
    return(maxWidth(tree) / maxDepth(tree))
  }
}
