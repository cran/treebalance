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
#' Note that the ranks grow very quickly with the number of leaves. Thus, they are computed
#' exactly and returned by default in \eqn{big.z} format (package \eqn{gmp}). With
#' \code{type="double"} the function stops for \eqn{n\geq 10}{n>=10}, since the ranks can
#' then no longer be represented exactly as double. The function also stops if the rank
#' would have more than about 5 million digits (e.g. caterpillar trees with more than 27 leaves).
#'
#' @param tree A rooted binary tree in phylo format.
#' @param type A character string specifying whether the rank is returned exactly as
#' big integer ("bigz", default, package \code{gmp}) or as "double" (only for \eqn{n\leq 9}{n<=9}).
#'
#' @return \code{colPlaLab} returns the Colijn-Plazzotta rank of the given tree.
#'
#' @author Sophie Kersting, Luise Kuehn
#'
#' @references C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, doi: 10.1093/sysbio/syx046.
#' @references N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 2021. doi: 10.1016/j.dam.2020.11.021.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' colPlaLab(tree)
#' colPlaLab(ape::read.tree(text="((,),(,(,)));"), type="double")
#'
#'@export
colPlaLab <- function(tree, type="bigz"){
  n <- length(tree$tip.label)
  if (!inherits(tree, "phylo"))
    stop("The input tree must be in phylo-format.")
  if (!(type %in% c("bigz", "double")))
    stop("The type must be either 'bigz' or 'double'.")
  if (type == "double" && n >= 10)
    stop(paste("For n >= 10 the Colijn-Plazzotta rank cannot be represented exactly",
               "as double. Please use type=\"bigz\"."))
  if (n == 1) {
    warning("The function might not deliver accurate results for n=1.")
    if (type == "double") return(1)
    return(gmp::as.bigz(1))
  }
  Descs <- getDescMatrix(tree)
  numbOfDescs <- sapply(1:(n+tree$Nnode),function(x) length(stats::na.omit(Descs[x,])))
  depthResults <- getNodesOfDepth(mat = Descs, root = n + 1, n = n)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  col_pla_labs <- gmp::as.bigz(rep(NA,n+tree$Nnode))
  if(is_binary(tree)){
    for (v in nodeorder) {
      if (numbOfDescs[v]==0) {
        col_pla_labs[v] <- 1
      }
      else {
        desc_cpl <- col_pla_labs[Descs[v,1:2]]
        if (desc_cpl[1] < desc_cpl[2]) desc_cpl <- rev(desc_cpl)
        # the number of digits roughly doubles in each step; stop before gmp runs out of memory
        if (gmp::sizeinbase(desc_cpl[1], 2) > 2^23)
          stop(paste("The Colijn-Plazzotta rank of this tree is too large to be computed",
                     "(it would have more than 5 million digits)."))
        col_pla_labs[v] <- gmp::divq.bigz(desc_cpl[1] * (desc_cpl[1] - 1), 2) + desc_cpl[2] + 1
      }
    }
  } else {
    stop("The tree has to be binary.")
  }
  if (type == "double") return(as.numeric(col_pla_labs[n+1]))
  return(col_pla_labs[n+1])
}
