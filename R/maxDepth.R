#' Calculation of the maximal depth of the tree
#'
#' This function calculates the maximal depth of any vertex in a rooted tree \eqn{T}, which
#' is at the same time its height \eqn{h(T)}. The tree must not necessarily be binary. Formally,
#' \eqn{h(T)} is defined as \deqn{h(T)=\max_{v\in V(T)} \delta(v)}{h(T)=max_{v in V(T)} depth(v)}
#' with \eqn{\delta(v)}{depth(v)} being the depth of the vertex \eqn{v}.
#' The maximal depth is an imbalance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{h(T)=0} and a warning.
#'
#' @param tree A rooted tree in phylo format.
#'
#' @return \code{maxDepth} returns the maximal depth, i.e. height, of a tree.
#'
#' @author Luise Kuehn, Sophie Kersting
#'
#' @references C. Colijn and J. Gardy.  Phylogenetic tree shapes resolve disease transmission patterns. Evolution, Medicine, and Public Health, 2014(1):96-108, 2014. ISSN 2050-6201. doi: 10.1093/emph/eou018.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' maxDepth(tree)
#' tree <- ape::read.tree(text="((,),((((,),),),(,)));")
#' maxDepth(tree)
#'
#'@export
maxDepth <- function(tree){
  if (!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n==1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }
  Descs <- getDescMatrix(tree)
  lastNodes <- n+1 #start with root
  current_depth <- 0
  while(length(lastNodes)>0){
    lastNodes <- stats::na.omit(as.vector(Descs[lastNodes,]))
    current_depth <- current_depth +1
  }
  return(current_depth-1)
}

