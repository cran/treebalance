% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colPlaLab.R
\name{colPlaLab}
\alias{colPlaLab}
\title{Calculation of the Colijn-Plazzotta rank for rooted binary trees}
\usage{
colPlaLab(tree, v = length(tree$tip.label) + 1)
}
\arguments{
\item{tree}{A rooted binary tree in phylo format.}

\item{v}{A vertex of the input tree. The Colijn-Plazzotta rank is calculated
for the subtree rooted at \eqn{v}. Default assumes that \eqn{v} is the root.}
}
\value{
\code{colPlaLab} returns the Colijn-Plazotta rank of the given (pending sub)tree.
}
\description{
This function calculates the Colijn-Plazzotta rank \eqn{CP(T)} for a
given rooted binary tree \eqn{T}. \eqn{CP(T)} is recursively defined as
\eqn{CP(T)=1} if \eqn{T} has only one leaf and otherwise
\deqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1} with
\eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the ranks of the two pending
subtrees rooted at the children of the root of \eqn{T}. The rank
of \eqn{T} corresponds to its position in the
lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),...
The Colijn-Plazzotta rank is an imbalance index.\cr\cr
For \eqn{n=1} the function returns \eqn{CP(T)=1} and a warning.\cr\cr
Note that problems can sometimes arise for trees with 20 leaves or more, due
to the limited range of computable values (ranks can reach INF quickly).
}
\examples{
tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
colPlaLab(tree)

}
\references{
C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, doi: 10.1093/sysbio/syx046.

N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 2021. doi: 10.1016/j.dam.2020.11.021.
}
\author{
Sophie Kersting
}
