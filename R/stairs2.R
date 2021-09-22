#' Calculation of the stairs2 value for rooted binary trees
#'
#' This function calculates the stairs2 value \eqn{st2(T)} for a given rooted
#' binary tree \eqn{T}. It is defined as the mean ratio between the leaf
#' numbers of the smaller and larger pending subtree over all inner vertices, more precisely
#' \deqn{st2(T)=\frac{1}{n-1}\cdot\sum_{u \in V_{in}(T)} \frac{n_{u_a}}{n_{u_b}}}{st2(T)=1/(n-1)*\sum_{u in V_in(T)} (n_ua/n_ub)}
#' in which \eqn{V_{in}(T)}{V_in(T)} denotes the set of all inner vertices
#' of \eqn{T}, and in which \eqn{n_{u_a}\geq n_{u_b}}{n_ua >= n_ub} denote the number of leaves
#' in the two pending subtrees that are
#' rooted at the direct descendants of \eqn{u}. The stairs2 value is an imbalance index. \cr\cr
#' Special cases: For \eqn{n=1}, the function returns \eqn{st2(T)=0} and a warning.
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{stairs2} returns the stairs2 value of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references C. Colijn and J. Gardy. Phylogenetic tree shapes resolve disease transmission patterns. Evolution, Medicine, and Public Health, 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' stairs2(tree)
#'
#'@export
stairs2 <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (!is_binary(tree))        stop("The input tree must be binary.")

  n <- length(tree$tip.label)
  if(n==1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  Descs <- getDescMatrix(tree)
  Ancs <- getAncVec(tree)
  depthResults <- getNodesOfDepth(mat=Descs, root=n+1, n=n)
  nv <- rep(NA,length(n+tree$Nnode))
  sum_of_ratios <- 0
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      nv[v] <- 1 #if leaf, then nv=1
    }else{
      desc_nvs <- nv[Descs[v,]]
      nv[v] <- sum(desc_nvs)
      sum_of_ratios <- sum_of_ratios+min(desc_nvs)/max(desc_nvs)
    }
  }
  return(sum_of_ratios/(n-1))
}
