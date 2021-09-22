#' Calculation of the maximal width of the tree
#'
#' This function calculates the maximal width \eqn{maxWidth(T)} for a
#' given rooted tree \eqn{T}. The tree must not necessarily be binary.
#' \eqn{maxWidth(T)} is defined as \deqn{maxWidth(T)=\max_{i=0,...,h(T)} w(i)}{maxWidth(T)=max_{i=0,...,h(T)} w(i)}
#' in which \eqn{h(T)} denotes the height of the tree \eqn{T} and \eqn{w(i)} denotes
#' the number of vertices in \eqn{T} that have depth \eqn{i}. The maximal width
#' is a balance index.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{maxWidth} returns the maximal width of a tree.
#'
#' @author Sophie Kersting
#'
#' @references C. Colijn and J. Gardy.  Phylogenetic tree shapes resolve disease transmission patterns. Evolution, Medicine, and Public Health, 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' maxWidth(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' maxWidth(tree)
#'
#'@export
maxWidth <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n==1){
    return(1)
  }
  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  widths <- sapply(1:(depthResults$maxdepth+1),
                   function(x) length(stats::na.omit(depthResults$nodesOfDepth[x,])))
  return(max(widths))
}
