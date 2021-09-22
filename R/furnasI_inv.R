#' Calculation of rooted binary tree for tuple (rank, leaf number)
#'
#' This function calculates the unique tree \eqn{T} (in phylo format) for two
#' given integer values \eqn{r} and \eqn{n}, with \eqn{n} denoting the number
#' of leaves of \eqn{T} and \eqn{r} denoting the rank of \eqn{T} in the
#' left-light rooted ordering of all rooted binary trees with \eqn{n} leaves.
#' It is the inverse function of \code{furnasI()}. For details on how to calculate
#' \eqn{T} (including algorithm) see "The generation of random, binary
#' unordered trees" by G.W. Furnas (1984) or "Tree balance indices: a comprehensive
#' survey" by Fischer et al. (2021).\cr\cr
#' \code{furnasI_inv} can be used e.g. to generate random rooted binary trees with a
#' certain number of leaves. Also, the concept of assigning each rooted binary
#' tree a unique tuple \eqn{(rank, n)} allows to store many trees with minimal
#' storage use.
#'
#' @param rank An integer denoting the rank of the sought tree among all rooted
#' binary trees with \eqn{n} leaves.
#' @param n An integer denoting the number of leaves of the sought tree.
#'
#' @return \code{furnasI_inv} returns the unique tree (in phylo format) for
#' the given leaf number and rank.
#'
#' @author Luise Kuehn
#'
#' @references G. W. Furnas. The generation of random, binary unordered trees. Journal of Classification, 1984. doi: 10.1007/bf01890123. URL https://doi.org/10.1007/bf01890123.
#'
#' @examples
#' furnasI_inv(rank=6,n=8)
#'
#'@export
furnasI_inv <- memoise::memoise(function(rank, n){
  if(rank < 1 | rank%%1!=0) stop("Tree cannot be calculated, because rank is not valid.")
  if(n%%1!=0)               stop("Tree cannot be calculated, because number of leaves is no integer.")
  if(rank > we_eth(n))      stop("Tree cannot be calculated, because rank is larger than the available number of trees.")

  # initial condition
  if(n == 1) {return(ape::read.tree(text="();"))}

  # calculate the sizes nL and nR of the two maximal strict pending subtrees
  nL <- 1
  we_sum <- we_eth(1)*we_eth(n-1)
  while(we_sum < rank){
    nL <- nL + 1
    we_sum <- we_sum + we_eth(nL)*we_eth(n-nL)
  }
  nR <- n - nL

  # define help variable
  h <- we_sum - we_eth(nL)*we_eth(n-nL)

  # if nL < nR: calculate the ranks of the two subtrees as follows
  if(nL < nR){
    rL <- floor((rank-1-h)/we_eth(nR) + 1)
    rR <- (rank-1-h)%%we_eth(nR) + 1
  }

  # if nL = nR: calculate the ranks of the two subtrees as follows
  if(nL == nR){
    rL <- 1 + max(floor(1/2*(2*we_eth(nL)-1-sqrt((2*we_eth(nL)-1)^2-8*
                                                   (rank-h-we_eth(nL))))),
                  floor(1/2*(2*we_eth(nL)+1-sqrt((2*we_eth(nL)+1)^2-8*(rank-h-1)))))
    rR <- rank - h - 1 - we_eth(nL)*(we_eth(nL)+1)/2 + (we_eth(nL)-rL+1)*(we_eth(nL)-rL+2)/2 + rL
  }

  # calculate the two subtrees recursively from their ranks and leaf numbers
  tL <- furnasI_inv(rank=rL, n=nL)
  tR <- furnasI_inv(rank=rR, n=nR)

  # merge the two subtrees to form the final tree
  return(tree_merge(tL, tR))
})
