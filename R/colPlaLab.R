#' Calculation of the Colijn-Plazzotta rank for rooted trees
#'
#' This function calculates the Colijn-Plazzotta rank \eqn{CP(T)} for a
#' given rooted tree \eqn{T}.\cr\cr
#' For a binary tree \eqn{T}, the Colijn-Plazzotta rank \eqn{CP(T)} is
#' recursively defined as \eqn{CP(T)=1} if \eqn{T} consists of only
#' one leaf and otherwise
#' \deqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1}
#' with \eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the ranks of the two pending
#' subtrees rooted at the children of the root of \eqn{T}. This rank
#' of \eqn{T} corresponds to its position in the
#' lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),...
#' The Colijn-Plazzotta rank of binary trees has been shown to be an imbalance index.\cr\cr
#' For an arbitrary tree \eqn{T} whose maximal number of children of any vertex is \eqn{l},
#' the Colijn-Plazzotta rank \eqn{CP(T)} is recursively defined as \eqn{CP(T)=0}
#' if \eqn{T} is the empty tree (with no vertices), \eqn{CP(T)=1} if \eqn{T} consists
#' of only one leaf and otherwise \eqn{CP(T)=\sum_{i=1}^l binom(CP(T_i)+i-1,i)}{CP(T)=\sum_{i=1,...,l} binom{CP(T_i)+i-1,i}}
#' (with \eqn{CP(T_1)\leq ...\leq CP(T_l)}{CP(T_1)<=...<=CP(T_l)}). If there are only \eqn{k<l} pending subtrees
#' rooted at the children of the root of \eqn{T}, then \eqn{T_1,...,T_{l-k}} are empty trees,
#' i.e. \eqn{CP(T_1)=...=CP(T_{l-k})=0}, and \eqn{CP(T_{l-k+1}),...,CP(T_l)} are the
#' increasingly ordered \eqn{CP}-ranks of the \eqn{k} pending subtrees rooted at the
#' children of the root of \eqn{T}. Note that if \eqn{k=l} there are no empty trees.\cr\cr
#' For \eqn{n=1} the function returns \eqn{CP(T)=1} and a warning.\cr\cr
#' Note that problems can sometimes arise even for trees with small leaf numbers due
#' to the limited range of computable values (ranks can reach INF quickly).
#'
#' @param tree A rooted tree in phylo format.
#' @param method The method must be one of: "binary" or "arbitrary". Note that
#' (only) in the arbitrary case vertices of out-degree 1 are allowed.
#'
#' @return \code{colPlaLab} returns the Colijn-Plazotta rank of the given tree
#' according to the chosen method.
#'
#' @author Sophie Kersting, Luise Kuehn
#'
#' @references C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, doi: 10.1093/sysbio/syx046.
#' @references N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 2021. doi: 10.1016/j.dam.2020.11.021.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' colPlaLab(tree, method="binary")
#' tree <- ape::read.tree(text="(((,),(,)),(,),(,));")
#' colPlaLab(tree, method="arbitrary")
#'
#'@export
colPlaLab <- function(tree, method){
  n <- length(tree$tip.label)
  if (!inherits(tree, "phylo"))
    stop("The input tree must be in phylo-format.")
  if (!(method %in% c("binary", "arbitrary")))
    stop("The method must be one of: binary or arbitrary.")
  if (n == 1) {
    warning("The function might not deliver accurate results for n=1.")
    return(1)
  }
  Descs <- getDescMatrix(tree)
  numbOfDescs <- sapply(1:(n+tree$Nnode),function(x) length(stats::na.omit(Descs[x,])))
  depthResults <- getNodesOfDepth(mat = Descs, root = n + 1, n = n)
  nv <- rep(NA, n + tree$Nnode)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  col_pla_labs <- rep(NA,n+tree$Nnode)
  if(method=="binary"){ # for binary trees only
    if(is_binary(tree)){
      for (v in nodeorder) {
        if (numbOfDescs[v]==0) {
          col_pla_labs[v] <- 1
        }
        else {
          desc_cpl <- sort(col_pla_labs[Descs[v,1:2]], decreasing = TRUE)
          col_pla_labs[v] <- desc_cpl[1] * (desc_cpl[1] - 1)/2 + desc_cpl[2] + 1
        }
      }
    } else {
      stop("The tree has to be binary if method 'binary' is selected.")
    }
  } else { # method for arbitrary trees
    maxdescs <- max(numbOfDescs)
    for (v in nodeorder) {
      if (numbOfDescs[v]==0) {
        col_pla_labs[v] <- 1
      }
      else {
        desc_cpl <- rep(0,maxdescs)
        desc_cpl[(maxdescs-numbOfDescs[v]+1):maxdescs] <- sort(
          col_pla_labs[Descs[v,1:numbOfDescs[v]]], decreasing = FALSE)
        col_pla_labs[v] <- sum(choose(desc_cpl+0:(maxdescs-1),1:maxdescs))
      }
    }
  }
  return(col_pla_labs[n+1])
}
