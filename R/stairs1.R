#' Calculation of the stairs1 value for rooted binary trees
#'
#' This function calculates the stairs1 value \eqn{st1(T)} for a given rooted
#' binary tree \eqn{T}. It is a modified version of the Rogers J index and is
#' defined as the fraction of inner vertices whose balance
#' value is unequal to zero, more precisely
#' \deqn{st1(T)=\frac{1}{n-1}\cdot\sum_{u \in V_{in}(T)} (1-I(n_{u_a}=n_{u_b}))}{st1(T)=1/(n-1)*\sum (1-I(n_ua=n_ub)) over all u in V_in(T)}
#' in which \eqn{V_{in}(T)}{V_in(T)} denotes the set of all inner vertices
#' of \eqn{T}, and in which \eqn{n_{u_a}}{n_ua}
#' and \eqn{n_{u_b}}{n_ub} denote the number of leaves in the two pending subtrees that are
#' rooted at the direct descendants of \eqn{u}. The stairs1 value is an imbalance index. \cr\cr
#' Special cases: For \eqn{n=1}, the function returns \eqn{st1(T)=0} and a warning. \cr\cr
#' For details on the stairs1 value, see 
#' also Chapter 23 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_23).
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{stairs1} returns the stairs1 value of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references M. M. Norstrom, M. C. Prosperi, R. R. Gray, A. C. Karlsson, and M. Salemi. PhyloTempo: A Set of R Scripts for Assessing and Visualizing Temporal Clustering in Genealogies Inferred from Serially Sampled Viral Sequences. Evolutionary Bioinformatics, 8:EBO.S9738, 2012. ISSN 1176-9343, 1176-9343. doi:10.4137/EBO.S9738.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' stairs1(tree)
#'
#'@export
stairs1 <- function(tree){
  if(!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if(!is_binary(tree))        stop("The input tree must be binary.")

  n <- length(tree$tip.label)
  if(n==1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  return(rogersI(tree)/(length(tree$tip.label)-1))
}
