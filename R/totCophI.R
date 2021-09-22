#' Calculation of the total cophenetic index for rooted trees
#'
#' This function calculates the total cophenetic index \eqn{TCI(T)} of a given
#' rooted tree \eqn{T}. The tree must not necessarily be binary. \eqn{TCI(T)}
#' is defined as \deqn{TCI(T)=\sum_{1\leq i<j\leq n} \delta(lca(i,j))=\sum_{u\in V_{in}(T)\setminus\{\rho\}} binom(n_u,2)}{TCI(T)=\sum_{1<=i<j<=n} depth(lca(i,j))
#' =\sum_{u in V'_in(T)} binom(n_u,2)} in which \eqn{\delta(lca(i,j))}{depth(lca(i,j))} denotes the depth of the last
#' common ancestor of the two leaves \eqn{i} and \eqn{j} and \eqn{V_{in}(T)\setminus\{\rho\}}{V'_in(T)}
#' denotes the set of all inner vertices exept the root and \eqn{n_u} denotes the
#' number of descendant leaves of \eqn{u}. The second formula is useful for efficient
#' computation of \eqn{TCI(T)}. The total cophenetic index is an imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{TCI(T)=0}.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{totCophI} returns the total cophenetic index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references A. Mir, F. Rosselló, and L. Rotger. A new balance index for phylogenetic trees. Mathematical Bio-sciences, 241(1):125-136, 2013. doi: 10.1016/j.mbs.2012.10.005.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' totCophI(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' totCophI(tree)
#' tree <- ape::read.tree(text="((,,,),(,,));")
#' totCophI(tree)
#'
#'@export
totCophI <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n == 1 || n==2) {return(0)}

  nv_vec <- get.subtreesize(tree)[(n+2):(n+tree$Nnode)]
  tci_val <- sapply(nv_vec, function(x) choose(x,2))
  return(sum(tci_val))
}


