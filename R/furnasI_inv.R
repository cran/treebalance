#' Calculation of rooted binary tree for tuple (rank, leaf number)
#'
#' This function calculates the unique tree \eqn{T} (in phylo format) for two
#' given integer values \eqn{r} and \eqn{n}, with \eqn{n} denoting the number
#' of leaves of \eqn{T} and \eqn{r} denoting the rank of \eqn{T} in the
#' left-light rooted ordering of all rooted binary trees with \eqn{n} leaves.
#' It is the inverse function of \code{furnasI()}. For details on how to calculate
#' \eqn{T} (including algorithm) see "The generation of random, binary
#' unordered trees" by G.W. Furnas (1984) or "Tree balance indices: a comprehensive
#' survey" by Fischer et al. (2023).\cr\cr
#' \code{furnasI_inv} can be used e.g. to generate random rooted binary trees with a
#' certain number of leaves. Also, the concept of assigning each rooted binary
#' tree a unique tuple \eqn{(rank, n)} allows to store many trees with minimal
#' storage use.
#'
#' @param rank An integer denoting the rank of the sought tree among all rooted
#' binary trees with \eqn{n} leaves. It can be given as \code{bigz} (package \code{gmp}),
#' as character or as double; the latter only up to \eqn{2^{53}}{2^53}, since larger integers
#' cannot be represented exactly as double.
#' @param n An integer denoting the number of leaves of the sought tree (at most 2545).
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
#' furnasI_inv(rank="100000000000000000000",n=60)
#'
#' @importFrom memoise memoise
#' @importFrom gmp as.bigz
#'@export
furnasI_inv <- memoise::memoise(function(rank, n){
  if (n<1 || n%%1 != 0)
    stop("Tree cannot be calculated, because number of leaves is no positive integer.")
  if (!gmp::is.bigz(rank) && !is.character(rank)) {
    if (!is.finite(rank) || rank < 1 || (rank <= 2^53 && rank%%1 != 0))
      stop("Tree cannot be calculated, because rank is not valid.")
    if (rank > 2^53)
      stop(paste("Ranks larger than 2^53 cannot be represented exactly as double.",
                 "Please enter the rank as bigz (package gmp) or as character."))
  }
  rank <- gmp::as.bigz(rank)
  if (is.na(rank) || rank < 1)
    stop("Tree cannot be calculated, because rank is not valid.")
  if (rank > we_eth(n))
    stop(paste("Tree cannot be calculated, because rank",as.character(rank),
               "is larger than the available number of trees",as.character(we_eth(n)),
               "for n =",n,"."))
  if (n == 1) {
    return(ape::read.tree(text = "();"))
  }
  we_nums_mult <- gmp::as.bigz(wedEth_chr()[1:ceiling(n/2)])*
    rev(gmp::as.bigz(wedEth_chr()[floor(n/2):(n-1)]))
  rsums <- cumsum(we_nums_mult)
  alpha <- min(which(rsums>=rank))
  if(alpha>1) {
    rsums_alpha1 <- rsums[alpha-1]
  } else {
    rsums_alpha1 <- gmp::as.bigz(0)
  }
  beta <- n-alpha
  if(alpha<beta){
    b_temp <- gmp::mod.bigz(rank-rsums_alpha1, we_eth(beta))
    a_temp <- gmp::divq.bigz(rank-rsums_alpha1-b_temp, we_eth(beta))
    if(b_temp>0){
      r_alpha <- a_temp + 1
      r_beta <- b_temp
    } else if(b_temp==0) {
      r_alpha <- a_temp
      r_beta <- we_eth(beta)
    }
  } else if(alpha==beta) {
    # r_alpha is the largest a in 1,...,we(beta) with
    # (a-1)*we(beta)-(a-1)*(a-2)/2 < rank-rsums_alpha1,
    # found by binary search (sqrt is not available for bigz)
    lower <- gmp::as.bigz(1)
    upper <- we_eth(beta)
    while(lower < upper) {
      mid <- gmp::divq.bigz(lower+upper+1, 2)
      if((mid-1)*we_eth(beta)-gmp::divq.bigz((mid-1)*(mid-2), 2) < rank-rsums_alpha1) {
        lower <- mid
      } else {
        upper <- mid-1
      }
    }
    r_alpha <- lower
    r_beta <- rank-rsums_alpha1-(r_alpha-1)*we_eth(beta)+
      gmp::divq.bigz((r_alpha-2)*(r_alpha-1), 2) +r_alpha -1
  }
  tL <- furnasI_inv(rank = r_alpha, n = alpha)
  tR <- furnasI_inv(rank = r_beta, n = beta)
  return(tree_merge(tL, tR))
})
