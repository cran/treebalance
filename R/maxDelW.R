#' Calculation of the (modified) maximum difference in widths for a rooted tree
#'
#' This function calculates the maximum difference in widths \eqn{delW(T)} and
#' the modified maximum difference in width \eqn{mdelW(T)} for a
#' given rooted tree \eqn{T}. The tree must not necessarily be binary.
#' \eqn{delW(T)} is defined as \deqn{delW(T)=\max_{i=0,...,h(T)-1} |w(i+1)-w(i)|}{delW(T)=max_{i=0,...,h(T)-1} |w(i+1)-w(i)|}
#' and \eqn{mdelW(T)} is defined as \deqn{mdelW(T)=\max_{i=0,...,h(T)-1} w(i+1)-w(i)}{mdelW(T)=max_{i=0,...,h(T)-1} w(i+1)-w(i)}
#' in which \eqn{h(T)} denotes the height of the tree \eqn{T} and \eqn{w(i)} denotes
#' the number of vertices in \eqn{T} that have depth \eqn{i}. The modified maximum difference
#' in widths is a balance index, while the maximum difference in widths is neither a balance nor imbalance index. \cr\cr
#' Note that there was a spelling error in the previous manual of this function - we wrote "maximum difference in widths" 
#' while the given definition and the R code corresponded to the "modified maximum difference in width". \cr\cr
#' For details on the maximum difference in widths and the modified maximum difference in widths, see 
#' also Chapters 24 and 23 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_24, https://doi.org/10.1007/978-3-031-39800-1_23).
#'
#' @param tree A rooted tree in phylo format.
#' @param method A character string specifying whether the original maximum difference in widths or the 
#' modified maximum difference in widths shall be computed. Can be any of "original" or "modified" (default is modified).
#'
#' @return \code{maxDelW} returns the maximum difference in widths of a tree (if \code{method} is set to \code{original}) 
#' or the modified maximum difference in widths (if \code{method} is set to \code{modified}).
#'
#' @author Sophie Kersting, Luise Kuehn
#'
#' @references C. Colijn, J. Gardy. Phylogenetic tree shapes resolve disease transmission patterns. Evolution, Medicine, and Public Health, 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' maxDelW(tree, method="original")
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' maxDelW(tree, method="modified")
#'
#' @export
maxDelW <- function(tree, method="modified"){
  
  message("Note that there was a spelling error in the previous manual of this function - we wrote 'maximum difference in widths' while the given definition and the R code corresponded to the 'modified maximum difference in width'.")
  
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if (!(method %in% c("original", "modified"))) stop("The method must bei either 'original' or 'modified'.")
  n <- length(tree$tip.label)
  
  if(n==1){
    return(0)
  }
  
  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  widths <- sapply(1:(depthResults$maxdepth+1),
                   function(x) length(stats::na.omit(depthResults$nodesOfDepth[x,])))
  if (method == "original") {
    return(max(abs(widths[2:(depthResults$maxdepth+1)]-widths[1:depthResults$maxdepth])))
  }
  if (method == "modified") {
    return(max(widths[2:(depthResults$maxdepth+1)]-widths[1:depthResults$maxdepth]))
  }
  
}
