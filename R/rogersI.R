#' Calculation of the Rogers J index for rooted binary trees
#'
#' This function calculates the Rogers J index \eqn{J(T)} for a given rooted
#' binary tree \eqn{T}. It is defined as the number of inner vertices whose balance
#' value is unequal to zero, more precisely
#' \deqn{J(T)=\sum_{u \in V_{in}(T)} (1-I(n_{u_a}=n_{u_b}))}{J(T)=\sum (1-I(n_{u_a}=n_{u_b})) over all u in V_in(T)}
#' in which \eqn{V_{in}(T)}{V_in(T)} denotes the set of all inner vertices
#' of \eqn{T}, and in which \eqn{n_{u_a}}{n_ua}
#' and \eqn{n_{u_b}}{n_ub} denote the number of leaves in the two pending subtrees that are
#' rooted at the direct descendants of \eqn{u}. \cr
#' Special cases: For \eqn{n=1}, the function returns \eqn{J(T)=0} and a warning. \cr\cr
#' For details on the Rogers J index, see 
#' also Chapter 19 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_19).
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{rogersI} returns the Rogers J index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references J. S. Rogers. Central Moments and Probability Distributions of Three Measures of Phylogenetic Tree Imbalance. Systematic Biology, 45(1):99-110, 1996. doi: 10.1093/sysbio/45.1.99.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' rogersI(tree)
#'
#' @export
rogersI <- function(tree){
  if(!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if(!is_binary(tree))        stop("The input tree must be binary.")

  n <- length(tree$tip.label)
  if(n==1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs, root=n+1, n=n)
  nv <- rep(NA, n+tree$Nnode)
  unbalanced_count <- 0
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      nv[v] <- 1 #if leaf, then nv=1
    }else{
      desc_nvs <- nv[Descs[v,1:2]]
      nv[v] <- sum(desc_nvs)
      if(desc_nvs[1]!=desc_nvs[2]){
        unbalanced_count <- unbalanced_count + 1
      }
    }
  }
  return(unbalanced_count)
}

