#' Calculation of the Furnas rank for rooted binary trees
#'
#' This function calculates the Furnas rank \eqn{F(T)} for a given rooted
#' binary tree \eqn{T}. \eqn{F(T)} is the unique rank of the tree \eqn{T}
#' among all rooted binary trees with \eqn{n} leaves in the left-light rooted
#' ordering. For details on the left-light rooted ordering as well as details
#' on how the Furnas rank is computed, see "The generation
#' of random, binary unordered trees" by G.W. Furnas (1984) or "Tree balance
#' indices: a comprehensive survey" by Fischer et al. (2023). The Furnas rank
#' is a balance index.\cr\cr
#' The concept of assigning each rooted binary tree a unique tuple \eqn{(rank, n)}
#' allows to store many trees with minimal storage use. \cr\cr
#' For details on the Furnas rank, see 
#' also Chapter 22 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_22).
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{furnasI} returns the unique Furnas rank of the given tree, i.e.
#' the rank of the tree among all rooted binary trees with \eqn{n} leaves in the
#' left-light rooted ordering. Since the values can get quite large, the function
#' returns them in \eqn{big.z} format (package \eqn{gmp}).
#'
#' @author Luise Kuehn, Lina Herbst
#'
#' @references G. W. Furnas. The generation of random, binary unordered trees. Journal of Classification, 1984. doi: 10.1007/bf01890123. URL https://doi.org/10.1007/bf01890123.
#' @references M. Kirkpatrick and M. Slatkin. Searching for evolutionary patterns in the shape of a phylogenetic tree. Evolution, 1993. doi: 10.1111/j.1558-5646.1993.tb02144.x.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' furnasI(tree)
#' 
#' @export
furnasI <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (!is_binary(tree))        stop("The input tree is not binary.")
  n <- length(tree$tip.label)
  
  # initial conditions
  if(n == 1) return(gmp::as.bigz(1))
  if(n == 2) return(gmp::as.bigz(1))
  
  # return the Furnas rank for the input tree
  allranks <- getfurranks(tree)
  return(allranks[n+1])
}
