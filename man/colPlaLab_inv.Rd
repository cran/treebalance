% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colPlaLab_inv.R
\name{colPlaLab_inv}
\alias{colPlaLab_inv}
\title{Generation of the rooted binary tree corresponding to a given Colijn-Plazzotta rank}
\usage{
colPlaLab_inv(rank)
}
\arguments{
\item{rank}{An integer denoting the Colijn-Plazzotta rank of the sought tree.}
}
\value{
\code{colPlaLab_inv} returns the unique rooted binary tree for the given rank.
}
\description{
This function generates the unique rooted binary tree \eqn{T} (in phylo
format) that corresponds to the given Colijn-Plazzotta rank \eqn{CP(T)}. It
is the inverse function of colPlaLab(). \cr\cr
\code{colPlaLab()}:
For a given rooted binary tree \eqn{T}, \eqn{CP(T)} is recursively defined as
\eqn{CP(T)=1} if \eqn{T} consists of only one vertex and otherwise
\eqn{CP(T)=\frac{1}{2}\cdot CP(T_1)\cdot(CP(T_1)-1)+CP(T_2)+1}{CP(T)=1/2*CP(T1)(CP(T1)-1)+CP(T2)+1} with
\eqn{CP(T_1) \geq CP(T_2)}{CP(T1)>=CP(T2)} being the
ranks of the two pending subtrees rooted at the children of the root of \eqn{T}.
The rank \eqn{CP(T)} of \eqn{T} corresponds to its position in the
lexicographically sorted list of (\eqn{i,j}): (1),(1,1),(2,1),(2,2),(3,1),... \cr\cr
\code{colPlaLab_inv()}:
For a given rank \eqn{CP} the corresponding tree \eqn{T} can be reconstructed
by starting from one vertex \eqn{\rho} (labelled \eqn{CP}) and recursively
splitting vertices whose labels \eqn{h} are greater than 1 into two children with the labels:
\deqn{i=\left\lceil\frac{1+\sqrt{8\cdot h-7}}{2}\right\rceil-1}{i=ceil((1+sqrt(8*h-7))/2)-1} and
\deqn{j=h-\frac{i\cdot(i-1)}{2}-1}{j=h-(i*(i-1))/2-1}
until there are no more vertices to split. \cr
For \eqn{CP=1} the function returns the smallest possible tree in the
phylo format: the tree consisting of a single edge.\cr\cr
Note that problems can arise for extremely high input values (>10e+18). \cr\cr
For details on the Colijn-Plazzotta rank, see 
also Chapter 21 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_21).
}
\examples{
colPlaLab_inv(22)

}
\references{
C. Colijn and G. Plazzotta. A Metric on Phylogenetic Tree Shapes. Systematic Biology, 67(1):113-126,2018. doi: 10.1093/sysbio/syx046.

N. A. Rosenberg. On the Colijn-Plazzotta numbering scheme for unlabeled binary rooted trees. Discrete Applied Mathematics, 291:88-98, 2021. doi: 10.1016/j.dam.2020.11.021.
}
\author{
Sophie Kersting
}
