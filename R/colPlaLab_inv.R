#' Generation of the rooted binary tree corresponding to a given Colijn-Plazzotta rank
#'
#' This function generates the unique rooted binary tree \eqn{T} (in phylo
#' format) that corresponds to the given Colijn-Plazzotta rank \eqn{CP(T)}. It
#' is the inverse function of colPlaLab() with the method 'binary'. \cr\cr
#' \code{colPlaLab()}:
#' For a given rooted binary tree \eqn{T}, \eqn{CP(T)} is recursively defined as
#' \eqn{CP(T)=1} if \eqn{T} consists of only one vertex and otherwise
#' \eqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1} with
#' \eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the
#' ranks of the two pending subtrees rooted at the children of \eqn{T}.
#' The rank \eqn{CP(T)} of \eqn{T} corresponds to its position in the
#' lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),... \cr\cr
#' \code{colPlaLab_inv()}:
#' For a given rank \eqn{CP} the corresponding tree \eqn{T} can be reconstructed
#' by starting from one vertex \eqn{\rho} (labelled \eqn{CP}) and recursively
#' splitting vertices whose labels \eqn{h} are greater than 1 into two children with the labels:
#' \deqn{i=\left\lceil\frac{1+\sqrt{8\cdot h-7}}{2}\right\rceil-1}{i=ceil((1+sqrt(8*h-7))/2)-1} and
#' \deqn{j=h-\frac{i\cdot(i-1)}{2}-1}{j=h-(i*(i-1))/2-1}
#' until there are no more vertices to split. \cr
#' For \eqn{CP=1} the function returns the smallest possible tree in the
#' phylo format: the tree consisting of a single edge.\cr\cr
#' Note that problems can arise for extremely high input values (>10e+18).
#'
#' @param rank An integer denoting the Colijn-Plazzotta rank of the sought tree.
#'
#' @return \code{colPlaLab_inv} returns the unique rooted binary tree for the given rank.
#'
#' @author Sophie Kersting
#'
#' @references C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, 67(1):113-126,2018. doi: 10.1093/sysbio/syx046.
#' @references N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 291:88-98, 2021. doi: 10.1016/j.dam.2020.11.021.
#'
#' @examples
#' colPlaLab_inv(22)
#'
#'@export
colPlaLab_inv <- function(rank){
  if(rank == 1) { # phylo format has no tree with a single node (only single edge)
    return(ape::read.tree(text="();"))
  }else{
    return(ape::read.tree(text=paste(cPL_inv(rank),";",sep = "")))
  }
}
