% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/furnasI_inv.R
\name{furnasI_inv}
\alias{furnasI_inv}
\title{Calculation of rooted binary tree for tuple (rank, leaf number)}
\usage{
furnasI_inv(rank, n)
}
\arguments{
\item{rank}{An integer denoting the rank of the sought tree among all rooted
binary trees with \eqn{n} leaves.}

\item{n}{An integer denoting the number of leaves of the sought tree.}
}
\value{
\code{furnasI_inv} returns the unique tree (in phylo format) for
the given leaf number and rank.
}
\description{
This function calculates the unique tree \eqn{T} (in phylo format) for two
given integer values \eqn{r} and \eqn{n}, with \eqn{n} denoting the number
of leaves of \eqn{T} and \eqn{r} denoting the rank of \eqn{T} in the
left-light rooted ordering of all rooted binary trees with \eqn{n} leaves.
It is the inverse function of \code{furnasI()}. For details on how to calculate
\eqn{T} (including algorithm) see "The generation of random, binary
unordered trees" by G.W. Furnas (1984) or "Tree balance indices: a comprehensive
survey" by Fischer et al. (2023).\cr\cr
\code{furnasI_inv} can be used e.g. to generate random rooted binary trees with a
certain number of leaves. Also, the concept of assigning each rooted binary
tree a unique tuple \eqn{(rank, n)} allows to store many trees with minimal
storage use. \cr\cr
For details on the Furnas rank, see 
also Chapter 22 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_22).
}
\examples{
furnasI_inv(rank=6,n=8)

}
\references{
G. W. Furnas. The generation of random, binary unordered trees. Journal of Classification, 1984. doi: 10.1007/bf01890123. URL https://doi.org/10.1007/bf01890123.
}
\author{
Sophie Kersting
}
