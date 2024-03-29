% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighL1dist.R
\name{weighL1dist}
\alias{weighL1dist}
\title{Calculation of weighted l1 distance index for rooted binary trees}
\usage{
weighL1dist(tree)
}
\arguments{
\item{tree}{A rooted binary tree in phylo format.}
}
\value{
\code{weighL1distI} returns the weighted l1 distance index of the given tree.
}
\description{
This function calculates the weighted l1 distance index \eqn{D_{l1}(T)}{Dl1(T)} for a
given rooted binary tree \eqn{T}. \eqn{D_{l1}(T)}{Dl1(T)} is defined as
\deqn{D_{l1}(T)=\sum_{z=2}^n z \cdot |f_n(z)-p_n(z)|}{Dl1(T)= \sum z*|f_n(z)-p_n(z)|
over all possible sizes 2<=z<=n} in which \eqn{n} denotes the
number of leaves of \eqn{T}, \eqn{f_n(z)} denotes the frequency of pending subtrees
of size \eqn{z} in \eqn{T} and \eqn{p_n(z)} is the expected number of
pending subtrees of size \eqn{z} under the Yule model, i.e. \eqn{p_n(z)=\frac{1}{n-1}}{p_n(z)=1/(n-1)}
if \eqn{z=n} and otherwise \eqn{\frac{n}{n-1}\cdot\frac{2}{z\cdot(z+1)}}{n/(n-1)*2/(z*(z+1))}.\cr\cr
For \eqn{n=1} the function returns \eqn{D_{l1}(T)=0}{Dl1(T)=0}. \cr\cr
For details on the weighted l1 distance index, see 
also Chapter 24 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_24).
}
\examples{
tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
weighL1dist(tree)

}
\references{
M. G. Blum and O. Francois. On statistical tests of phylogenetic tree imbalance: The Sackin and other indices revisited. Mathematical Biosciences, 195(2):141-153, 2005. doi: 10.1016/j.mbs.2005.03.003.
}
\author{
Sophie Kersting
}
