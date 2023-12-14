#' Calculation of the B2 index for rooted trees
#'
#' This function calculates the B2 index \eqn{B2(T)} for a given rooted
#' tree \eqn{T}. The tree must not necessarily be binary. \eqn{B2(T)} is defined as
#' \deqn{B2(T)=-\sum_{x\in V_L(T)} p_x\cdot log(p_x)}{
#' B2(T)=-\sum_{x in V_L(T)} p_x * log(p_x)} in which \eqn{V_L(T)} denotes the leaf
#' set of \eqn{T}, and in which \deqn{p_x=\prod_{v\in anc(x)} \frac{1}{|child(v)|}}{
#' p_x= \prod_{v in anc(x)} 1/|child(v)|} denotes
#' the probability of reaching leaf \eqn{x} when starting at the root and assuming
#' equiprobable branching at each vertex \eqn{v\in anc(x)}{v in anc(x)} with \eqn{anc(x)}
#' denoting the set of ancestors of \eqn{x} excluding
#' \eqn{x}. \eqn{child(v)} denotes the set of children of the inner vertex \eqn{v}.\cr
#' The \eqn{B2} index is a balance index.\cr\cr
#' For \eqn{n=1} the function returns \eqn{B2(T)=0} and a warning. \cr\cr
#' For details on the B2 index, see 
#' also Chapter 11 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_11).
#'
#' @param tree A rooted tree in phylo format.
#' @param logbase The base that shall be used for the logarithm. For binary
#' trees it is common to use base 2.
#'
#' @return \code{B2I} returns the B2 index of the given tree.
#'
#' @author Sophie Kersting, Luise Kuehn
#'
#' @references K.-T. Shao and R. R. Sokal. Tree Balance. Systematic Zoology, 39(3):266, 1990. \cr doi: 10.2307/2992186.
#' @references P.-M. Agapow and A. Purvis. Power of Eight Tree Shape Statistics to Detect Nonrandom Diversification: A Comparison by Simulation of Two Models of Cladogenesis. Systematic Biology, 51(6):866-872, 2002.doi: 10.1080/10635150290102564. \cr URL https://doi.org/10.1080/10635150290102564.
#' @references M. Hayati, B. Shadgar, and L. Chindelevitch. A new resolution function to evaluate tree shape statistics. PLOS ONE, 14(11):e0224197, 2019. doi: 10.1371/journal.pone.0224197.\cr URL https://doi.org/10.1371/journal.pone.0224197.
#' @references M. Kirkpatrick and M. Slatkin. Searching for evolutionary patterns in the shape of a phylogenetic tree. Evolution, 47(4):1171-1181, 1993. doi: 10.1111/j.1558-5646.1993.tb02144.x.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' B2I(tree)
#'
#'@export
B2I <- function(tree, logbase=2){
  #Check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  if (logbase <= 0)             stop("The logarithm base must be a positive number.")

  n <- length(tree$tip.label)
  if(n == 1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  numchild <- rep(0,n+tree$Nnode)
  for(i in 1:nrow(tree$edge)){
    k <- tree$edge[i,1]
    numchild[k] <- numchild[k] + 1
  }
  if(1 %in% numchild) stop("Vertices with out-degree 1 are not allowed.")

  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs, root=n+1, n=n)
  nodes_topdown <- stats::na.omit(as.vector(t(depthResults$nodesOfDepth)))
  Ancs <- getAncVec(tree)

  p_x <- rep(NA,n+tree$Nnode)
  p_x[n+1] <- 1 #p_x of root is one
  for(x in nodes_topdown[-1]){
    p_x[x] <- p_x[Ancs[x]]/numchild[Ancs[x]]
  }

  B2_val <- - sum(p_x[1:n] * log(p_x[1:n],base=logbase))
  return (B2_val)
}
