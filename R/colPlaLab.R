#' Calculation of the Colijn-Plazzotta rank for rooted binary trees
#'
#' This function calculates the Colijn-Plazzotta rank \eqn{CP(T)} for a
#' given rooted binary tree \eqn{T}.\cr\cr
#' For a binary tree \eqn{T}, the Colijn-Plazzotta rank \eqn{CP(T)} is
#' recursively defined as \eqn{CP(T)=1} if \eqn{T} consists of only
#' one leaf and otherwise
#' \deqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1}
#' with \eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the ranks of the two pending
#' subtrees rooted at the children of the root of \eqn{T}. This rank
#' of \eqn{T} corresponds to its position in the
#' lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),...
#' The Colijn-Plazzotta rank of binary trees has been shown to be an imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{CP(T)=1} and a warning.\cr\cr
#' Note that problems can sometimes arise even for trees with small leaf numbers due
#' to the limited range of computable values (ranks can reach INF quickly). \cr\cr
#' For details on the Colijn-Plazzotta rank, see 
#' also Chapter 21 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_21).
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{colPlaLab} returns the Colijn-Plazotta rank of the given tree. Since the values can get quite large, the
#' function returns them in big.z format (package gmp).
#'
#' @author Sophie Kersting, Luise Kuehn
#'
#' @references C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, doi: 10.1093/sysbio/syx046.
#' @references N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 2021. doi: 10.1016/j.dam.2020.11.021.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' colPlaLab(tree)
#'
#'@export
colPlaLab <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  if (!is_binary(tree))         stop("The tree has to be binary.")
  n <- length(tree$tip.label)
  
  if (n == 1) {
    warning("The function might not deliver accurate results for n=1.")
    return(gmp::as.bigz(1))
  }
  
  Descs <- getDescMatrix(tree)
  numbOfDescs <- sapply(1:(n+tree$Nnode),function(x) length(stats::na.omit(Descs[x,])))
  depthResults <- getNodesOfDepth(mat = Descs, root = n + 1, n = n)
  nv <- rep(NA, n + tree$Nnode)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  col_pla_labs <- rep(gmp::as.bigz(NA),n+tree$Nnode)
  for (v in nodeorder) {
    if (numbOfDescs[v]==0) {
      col_pla_labs[v] <- gmp::as.bigz(1)
    }
    else {
      descval1 <- col_pla_labs[Descs[v,1]]
      descval2 <- col_pla_labs[Descs[v,2]]
      if(descval1 > descval2) {
        desc_cpl <- c(descval1, descval2)
      } else {
        desc_cpl <- c(descval2, descval1)
      }
      col_pla_labs[v] <- gmp::add.bigz(gmp::div.bigz(gmp::mul.bigz(desc_cpl[1], gmp::sub.bigz(desc_cpl[1], 1)), 2), gmp::add.bigz(desc_cpl[2],1))
      #col_pla_labs[v] <- desc_cpl[1] * (desc_cpl[1] - 1)/2 + desc_cpl[2] + 1
    }
  }
  return(col_pla_labs[n+1])
}
