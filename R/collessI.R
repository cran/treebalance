#' Calculation of the Colless index for rooted binary trees
#'
#' This function calculates variants of the Colless index for a given rooted
#' binary tree \eqn{T}. All of them are imbalance indices.\cr\cr
#' The original Colless index \eqn{C(T)} is defined as
#' \deqn{C(T)=\sum_{u \in V_{in}(T)} |n_{u_a}-n_{u_b}|}{
#' C(T)=\sum_{u in V_in(T)} |n_{u_a}-n_{u_b}|}
#' in which \eqn{V_{in}(T)}{V_in(T)} denotes the set of all inner vertices
#' of \eqn{T}, and in which \eqn{n_{u_a}}
#' and \eqn{n_{u_b}} denote the number of leaves in the two pending subtrees that are
#' rooted at the direct descendants of \eqn{u}. \cr\cr
#' The corrected Colless index \eqn{I_C(T)} of \eqn{T} is defined as \eqn{I_C(T)=0} for
#' \eqn{n=1} and \eqn{n=2} and for \eqn{n>2} as
#' \deqn{I_C(T)=\frac{2\cdot C(T)}{(n-1)\cdot(n-2)}}{
#' I_C(T)=2*
#' C(T)/((n-1)*(n-2))}
#' in which \eqn{n} denotes the
#' total number of leaves in \eqn{T}. \cr\cr
#' The quadratic Colless index \eqn{QC(T)} of \eqn{T} is defined as
#' \deqn{QC(T)=\sum_{u\in V_{in}(T)} |n_{u_a}-n_{u_b}|^2}{
#' QC(T)=\sum_{u in V_in(T)} |n_{u_a}-n_{u_b}|^2} \cr\cr
#' Special cases: For \eqn{n=1} the function returns \eqn{C(T)=I_C(T)=QC(T)=0} and a warning.
#'
#' @param tree A rooted binary tree in phylo format.
#' @param method A character string specifying the version that shall be computed.
#' It can be one of the following: "original", "corrected", "quadratic"
#'
#' @return \code{collessI} returns the Colless index of the given tree according to the chosen method.
#'
#' @author Luise Kuehn and Sophie Kersting
#'
#' @references D. Colless. Review of Phylogenetics: the theory and practice of phylogenetic systematics. Systematic Zoology, 1982. ISSN 00397989.
#' @references T. M. Coronado, M. Fischer, L. Herbst, F. Rosselló, and K. Wicke. On the minimum value of the Colless index and the bifurcating trees that achieve it. Journal of Mathematical Biology, 2020.doi: 10.1007/s00285-020-01488-9.
#' @references S. B. Heard. Patterns in tree balance among cladistic, phenetic, and randomly generated phylogenetic trees. Evolution, 1992. doi: 10.1111/j.1558-5646.1992.tb01171.x.
#' @references K. Bartoszek, T. M. Coronado, A. Mir, and F. Rosselló. Squaring within the Colless index yields a better balance index. Mathematical Biosciences, 331:108503, 2021. doi: 10.1016/j.mbs.2020.108503.
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' collessI(tree, method="original")
#' collessI(tree, method="corrected")
#' collessI(tree, method="quadratic")
#'
#'@export
collessI <- function(tree, method="original"){
  #Check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  if(!(method %in% c("original", "corrected", "quadratic"))){
    stop("The method must be one of 'original', 'corrected' or 'quadratic'.")
  }

  n <- length(tree$tip.label)
  if(n==1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  Descs <- getDescMatrix(tree)
  if (ncol(Descs)>2){
    stop("The input tree must be binary.")
  }

  # summarize (quadratic) Colless values of all inner vertices
  numdescleaves <- get.subtreesize(tree=tree)
  col_val <- 0
  collessvector <- numdescleaves[Descs[(n+1):(n+tree$Nnode),1]]-
                      numdescleaves[Descs[(n+1):(n+tree$Nnode),2]]
  if(method == "quadratic") {
    col_val <- sum(collessvector * collessvector)
  } else {
    col_val <- sum(abs(collessvector))
  }
  # return the correct value according to the chosen method
  if(method == "corrected") {
    return(col_val / ((n-1)*(n-2)/2))
  }
  return(col_val)
}

