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
#' allows to store many trees with minimal storage use.\cr\cr
#' The Furnas rank can be computed for trees with at most 2545 leaves. With
#' \code{type="double"} the function stops for \eqn{n\geq 49}{n>=49}, since the ranks can
#' then no longer be represented exactly as double.
#'
#' @param tree A rooted binary tree in phylo format.
#' @param type A character string specifying whether the rank is returned exactly as
#' big integer ("bigz", default, package \code{gmp}) or as "double" (only for \eqn{n\leq 48}{n<=48}).
#'
#' @return \code{furnasI} returns the unique Furnas rank of the given tree, i.e.
#' the rank of the tree among all rooted binary trees with \eqn{n} leaves in the
#' left-light rooted ordering. Since the values can get quite large, the function
#' returns them by default in \eqn{big.z} format (package \eqn{gmp}).
#'
#' @author Luise Kuehn, Lina Herbst
#'
#' @references G. W. Furnas. The generation of random, binary unordered trees. Journal of Classification, 1984. doi: 10.1007/bf01890123. URL https://doi.org/10.1007/bf01890123.
#' @references M. Kirkpatrick and M. Slatkin. Searching for evolutionary patterns in the shape of a phylogenetic tree. Evolution, 1993. doi: 10.1111/j.1558-5646.1993.tb02144.x.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' furnasI(tree)
#' furnasI(tree, type="double")
#' @export
furnasI <- function(tree, type="bigz"){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (!(type %in% c("bigz", "double"))) stop("The type must be either 'bigz' or 'double'.")
  if (!is_binary(tree))        stop("The input tree is not binary.")
  n <- length(tree$tip.label)
  if (type == "double" && n >= 49)
    stop(paste("For n >= 49 the Furnas rank cannot be represented exactly as double.",
               "Please use type=\"bigz\"."))
  
  # initial conditions
  if(n == 1 || n == 2) {
    furrank <- gmp::as.bigz(1)
  } else {
    # get the Furnas rank for the input tree
    furrank <- getfurranks(tree)[n+1]
  }
  if(type == "double") return(as.numeric(furrank))
  return(furrank)
}
