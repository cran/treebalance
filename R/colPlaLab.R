#' Calculation of the Colijn-Plazzotta rank for rooted binary trees
#'
#' This function calculates the Colijn-Plazzotta rank \eqn{CP(T)} for a
#' given rooted binary tree \eqn{T}. \eqn{CP(T)} is recursively defined as
#' \eqn{CP(T)=1} if \eqn{T} has only one leaf and otherwise
#' \deqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1} with
#' \eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the ranks of the two pending
#' subtrees rooted at the children of the root of \eqn{T}. The rank
#' of \eqn{T} corresponds to its position in the
#' lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),...
#' The Colijn-Plazzotta rank is an imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{CP(T)=1} and a warning.\cr\cr
#' Note that problems can sometimes arise for trees with 20 leaves or more, due
#' to the limited range of computable values (ranks can reach INF quickly).
#'
#' @param tree A rooted binary tree in phylo format.
#' @param v A vertex of the input tree. The Colijn-Plazzotta rank is calculated
#' for the subtree rooted at \eqn{v}. Default assumes that \eqn{v} is the root.
#'
#' @return \code{colPlaLab} returns the Colijn-Plazotta rank of the given (pending sub)tree.
#'
#' @author Sophie Kersting
#'
#' @references C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, doi: 10.1093/sysbio/syx046.
#' @references N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 2021. doi: 10.1016/j.dam.2020.11.021.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' colPlaLab(tree)
#'
#'@export
colPlaLab <- function(tree, v=length(tree$tip.label)+1){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (!is_binary(tree))        stop("The input tree must be binary.")
  # calculates label of the pending subtree rooted in v
  if(v%%1!=0 || v<1 || v>(2*tree$Nnode+1)){
    stop("The input node is non-sensible. Note: The tree has to be binary.")
  }
  if(length(tree$tip.label) == 1) {
    warning("The function might not deliver accurate results for n=1.")
    return(1)
  }
  Descs <- getDescMatrix(tree)
  if(is.na(Descs[v,1])){
    return(1) #leaves have label 1
  }
  l1 <- colPlaLab(tree,Descs[v,1])
  l2 <- colPlaLab(tree,Descs[v,2])
  k <- max(l1,l2)
  j <- min(l1,l2)
  return (k*(k-1)/2+j+1)
}
