% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IbasedI.R
\name{IbasedI}
\alias{IbasedI}
\title{Calculation of the I-based indices for rooted trees}
\usage{
IbasedI(
  tree,
  specnum = rep(1, length(tree$tip.label)),
  method = "mean",
  correction = "none",
  logs = TRUE
)
}
\arguments{
\item{tree}{A rooted tree in phylo format (with possibly few polytomies).}

\item{specnum}{A vector whose \eqn{i}-th entry is the number of species that
the \eqn{i}-th leaf represents. (default is 1,...,1)}

\item{method}{A character string specifying the method that shall be used to
calculate \eqn{I(T)}. It can be one of the following: "root", "median",
"total", "mean", "quartdev".}

\item{correction}{A character string specifying the correction method that
shall be applied to the I values. It can be one of the following:
"none", "prime", "w".}

\item{logs}{Boolean value, (default true), determines if the number of suitable
nodes (i.e. nodes that fulfill the criteria) and polytomies in the tree should be printed.}
}
\value{
\code{IbasedI} returns an \eqn{I}-based balance index of the given tree according to the chosen (correction and) method.
}
\description{
This function calculates \eqn{I}-based indices \eqn{I(T)} for a given rooted
tree \eqn{T}. Note that the leaves of the tree may represent single species or
groups of more than one species. Thus, a vector is required that contains for
each leaf the number of species that it represents.
The tree may contain few polytomies, which are not allowed to concentrate in
a particular region of the tree (see p. 238 in Fusco and Cronk (1995)).\cr\cr
Let \eqn{v} be a vertex of \eqn{T} that fulfills the following criteria: a)
The number of descendant (terminal) species of \eqn{v} is \eqn{k_v>3}
(note that if each leaf represents only one species \eqn{k_v} is simply the
number of leaves in the pending subtree rooted at \eqn{v}), and
b) \eqn{v} has exactly two children.\cr\cr
Then, we can calculate the \eqn{I_v} value as follows:
\deqn{I_v=\frac{k_{v_a}-\left\lceil\frac{k_v}{2}\right\rceil}{k_v-1-\left\lceil\frac{k_v}{2}\right\rceil}}{I_v=(k_va-ceiling(k_v/2))/(k_v-1-ceiling(k_v/2))}
in which \eqn{k_{v_a}}{k_va} denotes the number of descendant (terminal) species
in the bigger one of the two pending subtrees rooted at \eqn{v}.\cr\cr
As the expected value of \eqn{I_v} under the Yule model depends on \eqn{k_v},
Purvis et al. (2002) suggested to take the corrected values \eqn{I'_v} or \eqn{I_v^w} instead.\cr
The \eqn{I'_v} value of \eqn{v} is defined as follows: \eqn{I'_v=I_v} if \eqn{k_v} is odd and \eqn{I'_v=\frac{k_v-1}{k_v}\cdot I_v}{(k_v-1)/k_v*I_v}
if \eqn{k_v} is even.\cr
The \eqn{I_v^w} value of \eqn{v} is defined as follows: \deqn{I_v^w=\frac{w(I_v)\cdot I_v}{mean_{V'(T)} w(I_v)}}{I_v^w=\frac{w(I_v)\cdot I_v}{mean_{V'(T)} w(I_v)}}
where \eqn{V'(T)} is the set of inner vertices of \eqn{T} that have precisely
two children and \eqn{k_v\geq 4}{k_v>=4}, and \eqn{w(I_v)} is a weight function with
\eqn{w(I_v)=1} if \eqn{k_v} is odd and \eqn{w(I_v)=\frac{k_v-1}{k_v}} if \eqn{k_v}
is even and \eqn{I_v>0}, and \eqn{w(I_v)=\frac{2\cdot(k_v-1)}{k_v}}{w(I_v)=2*(k_v-1)/k_v}
if \eqn{k_v} is even and \eqn{I_v=0}. \cr\cr
The \eqn{I}-based index of \eqn{T} can now be calculated using different methods.
Here, we only state the version for the \eqn{I'} correction method, but the non-corrected
version or the \eqn{I_v^w} corrected version works analoguously.
1) root: The \eqn{I'} index of \eqn{T} equals the \eqn{I'_v} value of the root of
\eqn{T}, i.e. \eqn{I'(T)=I'_{\rho}}{I'(T)=I'_\rho}, provided that the root fulfills the two
criteria. Note that this method does not fulfil the definition of an (im)balance index.
2) median: The \eqn{I'} index of \eqn{T} equals the median \eqn{I'_v} value of all
vertices \eqn{v} that fulfill the two criteria.
3) total: The \eqn{I'} index of \eqn{T} equals the summarised \eqn{I'_v} values of all
vertices \eqn{v} that fulfill the two criteria.
4) mean: The \eqn{I'} index of \eqn{T} equals the mean \eqn{I'_v} value of all
vertices \eqn{v} that fulfill the two criteria.
5) quartile deviation: The \eqn{I'} index of \eqn{T} equals the quartile
deviation (half the difference between third and first quartile) of the \eqn{I'_v} values of all
vertices \eqn{v} that fulfill the two criteria. \cr\cr
For details on the family of I-based indices, see 
also Chapter 17 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_17).
}
\examples{
tree <- ape::read.tree(text="(((((,),),),),);")
IbasedI(tree, method="mean")
IbasedI(tree, method="mean", correction="prime", specnum=c(1,1,2,1,1,1))

}
\references{
G. Fusco and Q. C. Cronk. A new method for evaluating the shape of large phylogenies. Journal of Theoretical Biology, 1995. doi: 10.1006/jtbi.1995.0136. URL https://doi.org/10.1006/jtbi.1995.0136.

A. Purvis, A. Katzourakis, and P.-M. Agapow. Evaluating Phylogenetic Tree Shape: Two Modifications to Fusco & Cronks Method. Journal of Theoretical Biology, 2002. doi: 10.1006/jtbi.2001.2443. URL https://doi.org/10.1006/jtbi.2001.2443.
}
\author{
Luise Kuehn and Sophie Kersting
}
