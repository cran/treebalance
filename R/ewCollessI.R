#' Calculation of the equal weights Colless index for rooted binary trees
#'
#' This function calculates the equal weights Colless index \eqn{I_2(T)} for a
#' given rooted binary tree \eqn{T}. \eqn{I_2(T)} is defined as
#' \deqn{I_2(T)=\frac{1}{n-2}\cdot\sum_{u\in V_{in}(T), n_u>2} \frac{|n_{u_a}-n_{u_b}|}{n_u-2}}{
#' I_2(T)=1/(n-2)* \sum |n_ua-n_ub|/(n_u-2) over all u in V_in(T) with n_u>2}
#' in which \eqn{V_{in}(T)}{V_in(T)} denotes the set of all inner vertices of \eqn{T},
#' and in which \eqn{n_u}, \eqn{n_{u_a}}{n_ua} and \eqn{n_{u_b}}{n_ub} denote the number of
#' leaves in the pending subtrees that are rooted at \eqn{u} and the two direct
#' descendants of \eqn{u}. The equal weights Colless index is an imbalance index.\cr\cr
#' For \eqn{n=1} and \eqn{n=2} the function returns \eqn{I_2(T)=0} and a warning.
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{ewCollessI} returns the equal weights Colless index of the given tree.
#'
#' @author Luise Kuehn
#'
#' @references A. O. Mooers and S. B. Heard. Inferring Evolutionary Process from Phylogenetic Tree Shape. The Quarterly Review of Biology, 72(1), 1997. doi: 10.1086/419657.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' ewCollessI(tree)
#'
#'@export
ewCollessI <- function(tree){
  #Check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  if(n==1 || n==2){
    warning("The function might not deliver accurate results for n=1 and n=2.")
    return(0)
  }
  Descs <- getDescMatrix(tree)
  # summarize weighted Colless values of all inner vertices with more than two descendant leaves
  numdescleaves <- get.subtreesize(tree=tree)
  ewcol_val <- 0
  for(i in (n+1):(n+tree$Nnode)){
    cur.descs <- stats::na.omit(Descs[i,])
    if(length(cur.descs) != 2) stop("The input tree must be binary.")
    if(numdescleaves[i] <= 2) next
    ewcol_val <- ewcol_val + abs(numdescleaves[cur.descs[1]] - numdescleaves[cur.descs[2]]) / (numdescleaves[i]-2)
  }

  # return the correct value
  return(ewcol_val / (n-2))
}
