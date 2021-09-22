#' Calculation of the cherry index for rooted trees
#'
#' This function calculates the cherry index \eqn{ChI(T)} for a
#' given rooted tree \eqn{T}. The tree must not necessarily be binary.
#' \eqn{ChI(T)} is defined as the number of cherries in the tree. A cherry
#' is a pair of leaves that have the same direct ancestor. Note, if a vertex \eqn{u}
#' has \eqn{x} leaves as direct descendants, the number of cherries induced by \eqn{u} is
#' \eqn{binom(x,2)}{binom(x,2)}. \cr\cr
#' The cherry index does not fulfill the definition
#' of an (im)balance index given in "Tree balance indices: a comprehensive survey"
#' (Fischer et al., 2021).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{cherryI} returns the cherry index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references A. McKenzie and M. Steel. Distributions of cherries for two models of trees. Mathematical Biosciences, 164(1):81-92, 2000. doi: 10.1016/s0025-5564(99)00060-7.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' cherryI(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' cherryI(tree)
#' tree <- ape::read.tree(text="((,,,),(,,));")
#' cherryI(tree)
#'
#'@export
cherryI <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n==1){
    return(0)
  }
  Descs <- getDescMatrix(tree)
  numb_cherries <- 0
  for(row in (n+1):(n+tree$Nnode)){
    numbOfDescLeaves <- sum(is.na(Descs[stats::na.omit(Descs[row,]),1]))
    if(numbOfDescLeaves>=2){
      numb_cherries <- numb_cherries + choose(numbOfDescLeaves,2)
    }
  }
  return(numb_cherries)
}

