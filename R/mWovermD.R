#' Calculation of the maximum width over maximum depth of the tree
#'
#' This function calculates the maximum width over maximum depth \eqn{mWovermD(T)} for a
#' given rooted tree \eqn{T}. The tree must not necessarily be binary. For \eqn{n>1},
#' \eqn{mWovermD(T)} is defined as \deqn{mWovermD(T)=maxWidth(T) / h(T)}
#' in which \eqn{h(T)} denotes the height of the tree \eqn{T}, which is the same as the 
#' maximum depth of any leaf in the tree, and \eqn{maxWidth(T)} denotes
#' the maximum width of the tree \eqn{T}. The maximum width over maximum depth
#' is a balance index.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{mWovermD} returns the maximum width over maximum depth of a tree.
#'
#' @author Luise Kuehn
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' mWovermD(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' mWovermD(tree)
#'
#'@export
mWovermD <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (maxDepth(tree) == 0) {
    stop("mWovermD cannot be computed for n=1.")
  } else {
    return(maxWidth(tree) / maxDepth(tree))
  }
}
