#' Calculation of the I-based indices for rooted trees
#'
#' This function calculates \eqn{I}-based indices \eqn{I(T)} for a given rooted
#' tree \eqn{T}. Note that the leaves of the tree may represent single species or
#' groups of more than one species. Thus, a vector is required that contains for
#' each leaf the number of species that it represents.
#' The tree may contain few polytomies, which are not allowed to concentrate in
#' a particular region of the tree (see p. 238 in Fusco(1995)).\cr\cr
#' Let \eqn{v} be a vertex of \eqn{T} that fulfills the following criteria: a)
#' The number of descendant (terminal) species of \eqn{v} is \eqn{k_v>3}
#' (note that if each leaf represents only one species \eqn{k_v} is simply the
#' number of leaves in the pending subtree rooted at \eqn{v}), and
#' b) \eqn{v} has exactly two children.\cr\cr
#' Then, we can calculate the \eqn{I_v} value as follows:
#' \deqn{I_v=\frac{k_{v_a}-\left\lceil\frac{k_v}{2}\right\rceil}{k_v-1-\left\lceil\frac{k_v}{2}\right\rceil}}{I_v=(k_va-ceiling(k_v/2))/(k_v-1-ceiling(k_v/2))}
#' in which \eqn{k_{v_a}}{k_va} denotes the number of descendant (terminal) species
#' in the bigger one of the two pending subtrees rooted at \eqn{v}.\cr\cr
#' As the expected value of \eqn{I_v} under the Yule model depends on \eqn{k_v},
#' Purvis et al. (2002) suggested to take the corrected values \eqn{I'_v} or \eqn{I_v^w} instead.\cr
#' The \eqn{I'_v} value of \eqn{v} is defined as follows: \eqn{I'_v=I_v} if \eqn{k_v} is odd and \eqn{I'_v=\frac{k_v-1}{k_v}\cdot I_v}{(k_v-1)/k_v*I_v}
#' if \eqn{k_v} is even.\cr
#' The \eqn{I_v^w} value of \eqn{v} is defined as follows: \deqn{I_v^w=\frac{w(I_v)\cdot I_v}{mean_{V'(T)} w(I_v)}}{I_v^w=\frac{w(I_v)\cdot I_v}{mean_{V'(T)} w(I_v)}}
#' where \eqn{V'(T)} is the set of inner vertices of \eqn{T} that have precisely
#' two children and \eqn{k_v\geq 4}{k_v>=4}, and \eqn{w(I_v)} is a weight function with
#' \eqn{w(I_v)=1} if \eqn{k_v} is odd and \eqn{w(I_v)=\frac{k_v-1}{k_v}} if \eqn{k_v}
#' is even and \eqn{I_v>0}, and \eqn{w(I_v)=\frac{2\cdot(k_v-1)}{k_v}}{w(I_v)=2*(k_v-1)/k_v}
#' if \eqn{k_v} is even and \eqn{I_v=0}. \cr\cr
#' The \eqn{I}-based index of \eqn{T} can now be calculated using different methods.
#' Here, we only state the version for the \eqn{I'} correction method, but the non-corrected
#' version or the \eqn{I_v^w} corrected version works analoguously.
#' 1) root: The \eqn{I'} index of \eqn{T} equals the \eqn{I'_v} value of the root of
#' \eqn{T}, i.e. \eqn{I'(T)=I'_{\rho}}{I'(T)=I'_\rho}, provided that the root fulfills the two
#' criteria. Note that this method does not fulfil the definition of an (im)balance index.
#' 2) median: The \eqn{I'} index of \eqn{T} equals the median \eqn{I'_v} value of all
#' vertices \eqn{v} that fulfill the two criteria.
#' 3) total: The \eqn{I'} index of \eqn{T} equals the summarised \eqn{I'_v} values of all
#' vertices \eqn{v} that fulfill the two criteria.
#' 4) mean: The \eqn{I'} index of \eqn{T} equals the mean \eqn{I'_v} value of all
#' vertices \eqn{v} that fulfill the two criteria.
#' 5) quartile deviation: The \eqn{I'} index of \eqn{T} equals the quartile
#' deviation (half the difference between third and first quartile) of the \eqn{I'_v} values of all
#' vertices \eqn{v} that fulfill the two criteria.
#'
#' @param tree A rooted tree in phylo format (with possibly few polytomies).
#' @param specnum A vector whose \eqn{i}-th entry is the number of species that
#' the \eqn{i}-th leaf represents. (default is 1,...,1)
#' @param method A character string specifying the method that shall be used to
#' calculate \eqn{I(T)}. It can be one of the following: "root", "median",
#' "total", "mean", "quartdev"
#' @param correction A character string specifying the correction method that
#' shall be applied to the I values. It can be one of the following:
#' "none", "prime", "w"
#' @param logs Boolean value, (default true), determines if the number of suitable
#' nodes (i.e. nodes that fulfill the criteria) and polytomies in the tree should be printed
#'
#' @return \code{IbasedI} returns an \eqn{I}-based balance index of the given tree according to the chosen (correction and) method.
#'
#' @author Luise Kuehn and Sophie Kersting
#'
#' @references G. Fusco and Q. C. Cronk. A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology, 1995. doi: 10.1006/jtbi.1995.0136. URL https://doi.org/10.1006/jtbi.1995.0136.
#' @references A. Purvis, A. Katzourakis, and P.-M. Agapow. Evaluating Phylogenetic Tree Shape: Two Modifications to Fusco & Cronks Method. Journal of Theoretical Biology, 2002. doi: 10.1006/jtbi.2001.2443. URL https://doi.org/10.1006/jtbi.2001.2443.
#'
#' @examples
#' tree <- ape::read.tree(text="(((((,),),),),);")
#' IbasedI(tree, method="mean")
#' IbasedI(tree, method="mean", correction="prime", specnum=c(1,1,2,1,1,1))
#'
#'@export
IbasedI <- function(tree, specnum=rep(1,length(tree$tip.label)), method="mean",
                    correction="none", logs=TRUE){
  # check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  n <- length(tree$tip.label)
  m <- tree$Nnode
  if (length(specnum) != n) stop("The vector specnum must have as many entries as the tree has leaves.")
  if(!(method%in%c("root","median","total","mean","quartdev"))) stop("The method must be one of: root, median, mean, total, quartdev.")
  if(!(correction%in%c("none","prime","w"))) stop("The correction must be one of: none, prime, w.")
  if (sum(specnum) < 4) stop("The tree must comprise at least 4 species.")

  # check for polytomies
  Descs <- getDescMatrix(tree)
  numchild <- rowSums(!is.na(Descs))
  if(sum(numchild[(n+1):(n+m)]!=2)>0 & logs) {
    warning(paste("The tree has",sum(numchild[(n+1):(n+m)]==1),"inner vertices with 1 direct descendant",
                  "and",sum(numchild[(n+1):(n+m)]>2),"polytomies."))
  }

  # get for each vertex the number of its descendant species
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  sub.sizes <- rep(NA,n+tree$Nnode)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      sub.sizes[v] <- specnum[v] #if leaf
    }else{
      sub.sizes[v] <- sum(sub.sizes[stats::na.omit(Descs[v,])])#if inner node
    }
  }
  suitable_nodes <- which(numchild==2 & sub.sizes>=4) #binary nodes with >3 desc species
  suits <- length(suitable_nodes)
  if (suits<1) stop("No binary vertices with >3 descendant species found.")
  if(logs) message(paste("Number of binary vertices with >3 descendant species:",suits))

  # calculate the I, I' and I^w value for each inner vertex v that fulfills the criteria
  I_val <- rep(NA,n+m)
  I_weigths <- rep(NA,n+m)
  for(v in suitable_nodes){
    n1 <- sub.sizes[Descs[v,1]]
    n2 <- sub.sizes[Descs[v,2]]
    I_val[v] <- (max(n1,n2)-ceiling((n1+n2)/2))/(n1+n2-1-ceiling((n1+n2)/2))
    if(correction == "prime"){ #weigh nodes with even sub.size
      if((n1+n2)%%2 == 0) {I_val[v] <- (n1+n2-1)/(n1+n2)*I_val[v]}
    }
    if(correction == "w"){ #get weights
      if((n1+n2)%%2 != 0) { #sub.size odd
        I_weigths[v] <- 1
      } else if(I_val[v]>0){ #sub.size even and I_v>0
        I_weigths[v] <- (n1+n2-1)/(n1+n2)
      } else {#sub.size even and I_v=0
        I_weigths[v] <- 2*(n1+n2-1)/(n1+n2)
      }
    }
  }
  if(correction == "w"){ #apply weights to I values
    I_val <- (I_weigths*I_val)/mean(I_weigths,na.rm = TRUE)
  }
  # return the correct value for each method
  if(method == "root") {
    if (is.na(I_val[n+1])){
      stop("For this method, the root must have exactly two direct descendants.")
    }
    return(I_val[n+1])
  }
  if(method == "median") {
    return(stats::median(I_val, na.rm=TRUE))
  }
  if(method == "mean") {
    return(mean(I_val, na.rm=TRUE))
  }
  if(method == "total") {
    return(sum(I_val, na.rm=TRUE))
  }
  if(method == "quartdev") {
    return(1/2*(stats::quantile(I_val,probs=0.75,na.rm=TRUE,names=FALSE)
                -stats::quantile(I_val,probs=0.25,na.rm=TRUE,names=FALSE)))
  }
}
