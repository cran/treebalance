% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/collesslikeI.R
\name{collesslikeI}
\alias{collesslikeI}
\title{Calculation of the Colless-like indices for rooted trees}
\usage{
collesslikeI(tree, f.size, dissim)
}
\arguments{
\item{tree}{A rooted binary tree in phylo format.}

\item{f.size}{A character string specifying the function \eqn{f} that shall be used to compute the f.size.
It can be one of the following: "exp", "ln" or the name of a function as a string.}

\item{dissim}{A character string specifying the dissimilarity that shall be
used. It can be one of the following: "mdm", "var", "sd" or the name of a function as a string.}
}
\value{
\code{collesslikeI} returns the Colless-like index of the given tree according to the chosen
function and dissimilarity.
}
\description{
This function calculates the Colless-like index for a given rooted
tree \eqn{T} according to the chosen weight function \eqn{f} and dissimilarity \eqn{D}.
The Colless-like index \eqn{CL(T)}
relative to \eqn{D} and \eqn{f} is the sum of the \eqn{(D,f)}-balance values
over all inner vertices of the tree. More precisely,
\deqn{CL(T)=\sum_{v\in V_{in}(T)} bal_{D,f}(v)}{CL(T)=\sum_{v in V_in(T)} bal_{D,f}(v)}
where \eqn{V_{in}(T)}{V_in(T)} is the
set of inner vertices of \eqn{T}. The \eqn{(D,f)}-balance value
of \eqn{v} with children \eqn{v_1,...,v_k}{v1,...,vk} is computed as
\deqn{bal_{D,f}(v)=D(fs(T_{v_1}),...,fs(T_{v_k}))}{bal_{D,f}(v)=D(fs(T_v1),...,fs(T_vk)))}
with \eqn{D} denoting the dissimilarity and \eqn{fs} denoting the f.size.\cr
The f.size \eqn{fs(T)} of a tree \eqn{T} uses the function \eqn{f}, which maps any
integer to a non-negative real number, to build a weighted sum of
the out-degrees of all vertices in \eqn{T}. More precisely,
\deqn{fs(T)=\sum_{v\in V(T)} f(deg+(v))}{fs(T)=\sum_{v in V(T)} f(deg+(v))}
where \eqn{V(T)} is the set of all
vertices of \eqn{T} and \eqn{deg+(v)} denotes the out-degree (i.e. the number of
children) of the vertex \eqn{v}. The \eqn{f}-functions that are already
implemented are \eqn{f(x)=e^x} and \eqn{f(x)=ln(x+e)}.\cr
The dissimilarity \eqn{D(x_1,...,x_k)} of a vector \eqn{x_1,...,x_k} assigns
a non-negative value to the vector, is independent of the order of the vector
entries and equals zero if and only if \eqn{x_1=...=x_k}. In this
implementation the following dissimilarity functions are already built-in:
mean deviation from the median (\eqn{mdm}),
the sample variance (\eqn{var}) and the sample standard deviation (\eqn{sd}).\cr
\code{collesslikeI} also allows the use of other functions for the weight function \eqn{f}
and the dissimilarity \eqn{D}.\cr\cr
Special cases: For \eqn{n=1} the function returns \eqn{CL(T)=0} and a warning. \cr\cr
For details on the family of Colless-like indices, see 
also Chapter 16 in "Tree balance indices: a comprehensive survey" (https://doi.org/10.1007/978-3-031-39800-1_16).
}
\examples{
tree <- ape::read.tree(text="((((,),),(,)),(((,),),(,)));")
collesslikeI(tree, f.size="exp", dissim="mdm")
collesslikeI(tree, f.size="exp", dissim="var")
collesslikeI(tree, f.size="ln", dissim="sd")
myfsize <- function(x) return(x+1)
mydissim <- function(x) return (var(x))
collesslikeI(tree, f.size="myfsize",dissim = "mydissim")

}
\references{
A. Mir, L. Rotger, and F. Rossello. Sound Colless-like balance indices for multifurcating trees. PLOSONE, 13(9):e0203401, 2018. doi: 10.1371/journal.pone.0203401
}
\author{
Luise Kuehn, Sophie Kersting
}
