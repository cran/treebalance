#' Wedderburn Etherington numbers (from OEIS)
#'
#' Contains a vector of Wedderburn Etherington numbers for \eqn{n=1} to \eqn{n=2545}.
#' Since these numbers grow very quickly, they are stored exactly as big integers
#' (\code{bigz} format, package \code{gmp}). Single values can also be accessed with
#' \code{we_eth(n)}, which optionally returns them as double (only for \eqn{n\leq 48}{n<=48}).
#'
#' @docType data
#'
#' @format \code{bigz} vector (package \code{gmp}) of length 2545
#'
#' @usage data(wedEth)
#'
#' @keywords datasets
#'
#' @source OEIS Sequence A001190 available at https://oeis.org/A001190
#'
#' @examples
#' data(wedEth)
#' wedEth[5]
#' wedEth[50]
"wedEth"
