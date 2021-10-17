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
#' @author Sophie Kersting
#'
#' @references G. W. Furnas. The generation of random, binary unordered trees. Journal of Classification, 1984. doi: 10.1007/bf01890123. URL https://doi.org/10.1007/bf01890123.
#'
#' @examples
#' furnasI_inv(rank=6,n=8)
#'
#'@export
furnasI_inv <- memoise::memoise(function(rank, n){
  if (rank < 1 || rank%%1 != 0)
    stop("Tree cannot be calculated, because rank is not valid.")
  if (n<1 || n%%1 != 0)
    stop("Tree cannot be calculated, because number of leaves is no positive integer.")
  if (rank > wedEth[n])
    stop(paste("Tree cannot be calculated, because rank",rank,
               "is larger than the available number of trees",wedEth[n],
               "for n =",n,"."))
  if (n == 1) {
    return(ape::read.tree(text = "();"))
  }
  we_nums_mult <- wedEth[1:ceiling(n/2)]*rev(wedEth[floor(n/2):(n-1)])
  rsums <- cumsum(we_nums_mult)
  alpha <- min(which(rsums>=rank))
  rsums_alpha1 <- ifelse(alpha>1,rsums[alpha-1],0)
  beta <- n-alpha
  if(alpha<beta){
    b_temp <- (rank-rsums_alpha1)%%wedEth[beta]
    a_temp <- (rank-rsums_alpha1-b_temp)/wedEth[beta]
    if(b_temp>0){
      r_alpha <- a_temp + 1
      r_beta <- b_temp
    } else if(b_temp==0) {
      r_alpha <- a_temp
      r_beta <- wedEth[beta]
    }
  } else if(alpha==beta) {
    temp_val <- 2*(wedEth[beta]-rank+rsums_alpha1)+ (1-2*wedEth[beta])^2/4
    m <- max(c(0,ceiling(wedEth[beta]-1/2-sqrt(temp_val))))
    r_alpha <- m+1
    r_beta <- rank-rsums_alpha1-(r_alpha-1)*wedEth[beta]+
      (r_alpha-2)*(r_alpha-1)/2 +r_alpha -1
  }
  tL <- furnasI_inv(rank = r_alpha, n = alpha)
  tR <- furnasI_inv(rank = r_beta, n = beta)
  return(tree_merge(tL, tR))
})
