#' Calculation of the Colless-like indices for rooted trees
#'
#' This function calculates the Colless-like index for a given rooted
#' tree \eqn{T} according to the chosen weight function \eqn{f} and dissimilarity \eqn{D}.
#' The Colless-like index \eqn{CL(T)}
#' relative to \eqn{D} and \eqn{f} is the sum of the \eqn{(D,f)}-balance values
#' over all inner vertices of the tree. More precisely,
#' \deqn{CL(T)=\sum_{v\in V_{in}(T)} bal_{D,f}(v)}{CL(T)=\sum_{v in V_in(T)} bal_{D,f}(v)}
#' where \eqn{V_{in}(T)}{V_in(T)} is the
#' set of inner vertices of \eqn{T}. The \eqn{(D,f)}-balance value
#' of \eqn{v} with children \eqn{v_1,...,v_k}{v1,...,vk} is computed as
#' \deqn{bal_{D,f}(v)=D(fs(T_{v_1}),...,fs(T_{v_k}))}{bal_{D,f}(v)=D(fs(T_v1),...,fs(T_vk)))}
#' with \eqn{D} denoting the dissimilarity and \eqn{fs} denoting the f.size.\cr
#' The f.size \eqn{fs(T)} of a tree \eqn{T} uses the function \eqn{f}, which maps any
#' integer to a non-negative real number, to build a weighted sum of
#' the out-degrees of all vertices in \eqn{T}. More precisely,
#' \deqn{fs(T)=\sum_{v\in V(T)} f(deg+(v))}{fs(T)=\sum_{v in V(T)} f(deg+(v))}
#' where \eqn{V(T)} is the set of all
#' vertices of \eqn{T} and \eqn{deg+(v)} denotes the out-degree (i.e. the number of
#' children) of the vertex \eqn{v}. The \eqn{f}-functions that are already
#' implemented are \eqn{f(x)=e^x} and \eqn{f(x)=ln(x+e)}.\cr
#' The dissimilarity \eqn{D(x_1,...,x_k)} of a vector \eqn{x_1,...,x_k} assigns
#' a non-negative value to the vector, is independent of the order of the vector
#' entries and equals zero if and only if \eqn{x_1=...=x_k}. In this
#' implementation the following dissimilarity functions are already built-in:
#' mean deviation from the median (\eqn{mdm}),
#' the sample variance (\eqn{var}) and the sample standard deviation (\eqn{sd}).\cr
#' \code{collesslikeI} also allows the use of other functions for the weight function \eqn{f}
#' and the dissimilarity \eqn{D}.\cr\cr
#' Special cases: For \eqn{n=1} the function returns \eqn{CL(T)=0} and a warning.
#'
#' @param tree A rooted binary tree in phylo format.
#' @param f.size A character string specifying the function \eqn{f} that shall be used to compute the f.size.
#' It can be one of the following: "exp", "ln" or the name of a function as a string.
#' @param dissim A character string specifying the dissimilarity that shall be
#' used. It can be one of the following: "mdm", "var", "sd" or the name of a function as a string.
#'
#' @return \code{collesslikeI} returns the Colless-like index of the given tree according to the chosen
#' function and dissimilarity.
#'
#' @author Luise Kuehn, Sophie Kersting
#'
#' @references A. Mir, L. Rotger, and F. Rosselló. Sound Colless-like balance indices for multifurcating trees. PLOSONE, 13(9):e0203401, 2018. doi: 10.1371/journal.pone.0203401
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' collesslikeI(tree, f.size="exp", dissim="mdm")
#' collesslikeI(tree, f.size="exp", dissim="var")
#' collesslikeI(tree, f.size="ln", dissim="sd")
#' myfsize <- function(x) return(x+1)
#' mydissim <- function(x) return (var(x))
#' collesslikeI(tree, f.size="myfsize",dissim = "mydissim")
#'
#' @export
collesslikeI <- function(tree, f.size, dissim) {
  #Check for errors in input
  if(!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  if(!(f.size %in% c("exp", "ln")) && !exists(f.size, mode = "function")){
    stop("The f.size must be one of the strings 'exp' or 'ln' or the name of an existing function as a string.")}
  if(!(dissim %in% c("mdm", "var", "sd")) && !exists(dissim, mode = "function")){
    stop("The dissimilarity must be one of 'mdm', 'var' or 'sd' or the name of an existing function as a string.")}

  #possible functions for computation of f.size are f(x)=e^x and f(x)=ln(x+e)
  if(f.size=="exp") {
    f.size <- function(x) return(exp(x))
  }else if(f.size=="ln") {
    f.size <- function(x) return(log(x+exp(1)))
  }else {
    f.size <- get(f.size)
  }

  #possible dissimilarities are mean deviation from median, sample variance and sample standard deviation
  if(dissim=="mdm") {
    dissim <- function(x) return(sum(abs(x-stats::median(x))))/length(x)
  } else if(dissim=="var") {
    dissim <- function(x) return(sum((x-mean(x))^2)/(length(x)-1))
  } else if(dissim=="sd") {
    dissim <- function(x) return(sqrt(sum((x-mean(x))^2)/(length(x)-1)))
  } else {
    dissim <- get(dissim)
  }

  n <- length(tree$tip.label)
  if(n==1) {
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  #get outdegree of each vertex
  desc_mat <- getDescMatrix(tree)
  outdegs <- sapply(1:nrow(desc_mat), function(x) sum(!is.na(desc_mat[x,])))

  #get f.size of each subtree
  depthResults <- getNodesOfDepth(mat=desc_mat, root=n+1, n=n)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  fsizes <- rep(NA,n+tree$Nnode)
  for(i in nodeorder) {
    #if i is a leaf
    if(outdegs[i]==0) fsizes[i] <- f.size(0)
    #if i is an inner vertex
    if(outdegs[i]>0)  fsizes[i] <- sum(fsizes[stats::na.omit(desc_mat[i,])]) +
        f.size(outdegs[i])
  }

  #get balance value of each inner vertex
  balvalues <- rep(NA,n+tree$Nnode)
  for(i in nodeorder) {
    #if i is a leaf
    if(outdegs[i]==0) next
    #if i is an inner vertex
    if(outdegs[i]>0)  balvalues[i] <- dissim(fsizes[stats::na.omit(desc_mat[i,])])
  }

  #get Colless-like index
  return(sum(stats::na.omit(balvalues)))
}
