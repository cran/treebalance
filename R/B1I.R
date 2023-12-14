#' Calculation of the B1 index for rooted trees
#'
#' This function calculates the \eqn{B1} index \eqn{B1(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{B1(T)} is defined as
#' \deqn{B1(T)=\sum_{u\in V_{in}(T)\setminus\{\rho\}} h(T_u)^{-1}}{
#' B1(T)=\sum_{u in V'_in(T)} h(T_u)^(-1)} in which
#' \eqn{V_{in}(T)\setminus\{\rho\}}{V'_in(T)} denotes the set of inner vertices of \eqn{T} without the root, and
#' \eqn{h(T_u)} denotes the height of the pending subtree rooted at \eqn{u}.
#' When restricted to binary trees, the \eqn{B1} index is a balance index. For
#' arbitrary trees it does not fulfill the definition of an (im)balance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{B1(T)=0} and a warning. \cr\cr
#' For details on the B1 index, see 
#' also Chapter 10 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_10).
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{B1I} returns the B1 index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references K.-T. Shao and R. R. Sokal. Tree Balance. Systematic Zoology, 39(3):266, 1990. \cr doi: 10.2307/2992186.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' B1I(tree)
#'
#'@export
B1I <- function(tree){
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")

  n <- length(tree$tip.label)
  if(n == 1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }
  if(n == 2) return(0)

  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  maxDepthSubtree <- rep(NA,length(n+tree$Nnode))
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      maxDepthSubtree[v] <- 0 #if leaf
    }else{
      maxDepthSubtree[v] <- max(maxDepthSubtree[stats::na.omit(Descs[v,])])+1
    }
  }
  B1_val <- sum(1/maxDepthSubtree[(n+2):(n+tree$Nnode)])
  return(B1_val)
}
