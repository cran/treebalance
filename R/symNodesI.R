#' Calculation of the symmetry nodes index for rooted binary trees
#'
#' This function calculates the symmetry nodes index \eqn{SNI(T)} for a given rooted
#' binary tree \eqn{T}. \eqn{SNI(T)} is defined as the number of inner vertices \eqn{v} that are not
#' symmetry nodes, i.e. the two pending subtrees rooted at the children of \eqn{v} do not
#' have the same tree shape.\cr\cr
#' For \eqn{n=1} the function returns \eqn{SNI(T)=0} and a warning.
#'
#' @param tree A rooted binary tree in phylo format.
#'
#' @return \code{symNodesI} returns the symmetry nodes index of the given tree.
#'
#' @author Sophie Kersting
#'
#' @references S. J. Kersting and M. Fischer. Measuring tree balance using symmetry nodes — A new balance index and its extremal properties. Mathematical Biosciences, page 108690, 2021.  ISSN 0025-5564.  doi:https://doi.org/10.1016/j.mbs.2021.108690
#'
#' @examples
#' tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
#' symNodesI(tree)
#'
#'@export
symNodesI <- function(tree){
  if(!inherits(tree,"phylo")) stop("The input tree must be in phylo-format.")
  if(!is_binary(tree))        stop("The input tree must be binary.")

  n <- length(tree$tip.label)
  if(n==1){
    warning("The function might not deliver accurate results for n=1.")
    return(0)
  }

  Descs <- getDescMatrix(tree)
  Ancs <- getAncVec(tree)
  depthResults <- getNodesOfDepth(mat=Descs, root=n+1, n=n)
  worklab <- matrix(rep(NA,2*(n+tree$Nnode)), ncol = 2)
  inum <- rep(NA, n+tree$Nnode)

  #calcualte working-labels and i-numbers
  for(d in depthResults$maxdepth:0){
    current_Nodes <- stats::na.omit(as.vector(depthResults$nodesOfDepth[d+1,]))
    for(v in current_Nodes){
      if(is.na(Descs[v,1])){ #if leaf
        worklab[v,] <- c(0,0)
      }else{
        worklab[v,] <-sort(c(inum[Descs[v,]]),decreasing = FALSE)
      }
    }
    inum[current_Nodes] <- symBucketLexicoSort(worklab[current_Nodes,])
  }

  #Counting symmetry nodes
  numb_symNodes <- 0
  for(i in (n+1):(n+tree$Nnode)){
    if(worklab[i,1]!=0){ #if not a leaf
      #check if worklab has two times the same entry
      if(worklab[i,1]==worklab[i,2]){
        numb_symNodes <- numb_symNodes + 1
      }
    }
  }
  return(n-1-numb_symNodes) # return number of nodes that are not sym. nodes
}
