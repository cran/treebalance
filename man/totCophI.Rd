% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/totCophI.R
\name{totCophI}
\alias{totCophI}
\title{Calculation of the total cophenetic index for rooted trees}
\usage{
totCophI(tree)
}
\arguments{
\item{tree}{A rooted tree in phylo format.}
}
\value{
\code{totCophI} returns the total cophenetic index of the given tree.
}
\description{
This function calculates the total cophenetic index \eqn{TCI(T)} of a given
rooted tree \eqn{T}. The tree must not necessarily be binary. \eqn{TCI(T)}
is defined as \deqn{TCI(T)=\sum_{1\leq i<j\leq n} \delta(lca(i,j))=\sum_{u\in V_{in}(T)\setminus\{\rho\}} binom(n_u,2)}{TCI(T)=\sum_{1<=i<j<=n} depth(lca(i,j))
=\sum_{u in V'_in(T)} binom(n_u,2)} in which \eqn{\delta(lca(i,j))}{depth(lca(i,j))} denotes the depth of the lowest
common ancestor of the two leaves \eqn{i} and \eqn{j} and \eqn{V_{in}(T)\setminus\{\rho\}}{V'_in(T)}
denotes the set of all inner vertices exept the root and \eqn{n_u} denotes the
number of descendant leaves of \eqn{u}. The second formula is useful for efficient
computation of \eqn{TCI(T)}. The total cophenetic index is an imbalance index.\cr\cr
For \eqn{n=1} the function returns \eqn{TCI(T)=0}. \cr\cr
For details on the total cophenetic index, see 
also Chapter 8 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_8).
}
\examples{
tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
totCophI(tree)
tree <- ape::read.tree(text="((,),((((,),),),(,)));")
totCophI(tree)
tree <- ape::read.tree(text="((,,,),(,,));")
totCophI(tree)

}
\references{
A. Mir, F. Rossello, and L. Rotger. A new balance index for phylogenetic trees. Mathematical Bio-sciences, 241(1):125-136, 2013. doi: 10.1016/j.mbs.2012.10.005.
}
\author{
Sophie Kersting
}
