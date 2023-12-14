#' Auxiliary functions
#'
#' \code{getDescMatrix} - Creates a matrix that contains the descendants of
#' node \eqn{i} in row \eqn{i}.
#'
#' @author Sophie Kersting, Luise Kuehn and Lina Herbst
#' @param tree A rooted tree in phylo format, >= 2 leaves
#' @return \code{desc_mat} numeric matrix
#' @export
#' @rdname auxFuncs
#' @examples
#' mat <- cbind(c(7,7,6,5,5,6),c(1,2,3,4,6,7))
#' tree <- list(edge=mat, tip.label=c("","","",""), Nnode=3)
#' getDescMatrix(tree)
#' mat <- cbind(c(5,5,5,5),c(1,2,3,4))
#' tree <- list(edge=mat, tip.label=c("","","",""), Nnode=1)
#' getDescMatrix(tree)
getDescMatrix <- function(tree){
  n <- length(tree$tip.label)
  m <- tree$Nnode
  descs <-  rep(NA,m+n-1)
  descs_index <-  rep(1,m+1) #where their descs start
  for(i in 1:nrow(tree$edge)){
    source <- tree$edge[i,1]
    descs_index <- descs_index + c(rep(0,source-n),rep(1,m+n-source+1))
  }
  cur_index <- rep(0,m)
  for(i in 1:nrow(tree$edge)){ # transfer every row (edge) to descs
    source <- tree$edge[i,1]
    descs[descs_index[source-n]+cur_index[source-n]] <- tree$edge[i,2]
    cur_index[source-n] <- cur_index[source-n]+1
  }
  numb_descs <- descs_index[2:(m+1)]-descs_index[1:m]
  desc_mat <- matrix(rep(NA,max(numb_descs)*(tree$Nnode+n)),
                     ncol=max(numb_descs))
  for(i in 1:m){ # transfer descs into desc_mat
    desc_mat[n+i,1:numb_descs[i]] <- descs[descs_index[i]:(descs_index[i+1]-1)]
  }
  return(desc_mat)
}
#' Auxiliary functions
#'
#' \code{getAncVec} - Creates a vector that contains the parent (direct ancestor) of
#' node \eqn{i} at position \eqn{i}.
#'
#' @return \code{anc_vec} numeric vector
#' @rdname auxFuncs
#' @export
#' @examples
#' getAncVec(tree)
getAncVec <- function(tree){
  n <- length(tree$tip.label)
  anc_vec <- rep(NA,tree$Nnode+n)
  for(i in 1:nrow(tree$edge)){ # transfer every row (edge) to anc_vec
    anc_vec[tree$edge[i,2]] <- tree$edge[i,1]
  }
  return(anc_vec)
}
#' Auxiliary functions
#'
#' \code{getNodesOfDepth} - Creates a matrix that contains the nodes of
#' depth \eqn{i} in row \eqn{i}.
#'
#' @param mat Descendants matrix from \code{getDescMatrix}
#' @param root Number (label) of the root of the tree
#' @param n Number of leaves of the tree
#' @return \code{nodes_of_depth} numeric matrix
#' @rdname auxFuncs
#' @export
#' @examples
#' getNodesOfDepth(mat=getDescMatrix(tree),root=length(tree$tip.label)+1,
#' n=length(tree$tip.label))
getNodesOfDepth <- function(mat,root,n){
  nodesOfDep <- matrix(rep(NA,n*n), ncol = n) #maxdepth=n-1
  depthOfNodes <- rep(NA,2*n-1)
  lastNodes <- root
  row_lengths <- rep(NA,n)
  row_lengths[1] <- 1
  current_depth <- 0
  while(length(lastNodes)>0){
    nodesOfDep[current_depth+1,1:length(lastNodes)] <- lastNodes
    depthOfNodes[lastNodes] <- current_depth
    lastNodes <- stats::na.omit(as.vector(mat[lastNodes,]))
    current_depth <- current_depth +1
    row_lengths[current_depth] <- length(lastNodes)
  }
  nodesOfDep <- nodesOfDep[1:(current_depth),1:max(row_lengths, na.rm = TRUE)]
  return(list(nodesOfDepth=nodesOfDep, maxdepth=current_depth-1,
              nodeDepths=depthOfNodes))
}
#' Auxiliary functions
#'
#' \code{symBucketLexicoSort} - Sorts the pairs of numbers lexicographically and
#' returns ranking. Uses bucket sort.
#'
#' @param workLabs numeric matrix (2 columns)
#' @return \code{ranking} numeric vector
#' @export
#' @rdname auxFuncs
#' @examples
#' myWorkLabs <- cbind(c(0,1,2,3,1,0),c(0,2,2,4,1,0))
#' symBucketLexicoSort(myWorkLabs)
symBucketLexicoSort <- function(workLabs){
  rows <- nrow(workLabs)
  if(is.null(rows)){
    return(1)
  }
  bigBuckets_numb <- max(workLabs[,1])+1
  tinyBuckets_numb <- max(workLabs[,2])+1
  ranks <- rep(NA,nrow(workLabs))
  #-----------
  bigBucket <- matrix(rep(NA,bigBuckets_numb*rows),ncol = rows)
  bBFill <- rep(1,bigBuckets_numb)
  for(i in 1:rows){
    firstEntry <- workLabs[i,1]
    bigBucket[firstEntry+1,bBFill[firstEntry+1]] <- i
    bBFill[firstEntry+1] <- bBFill[firstEntry+1] + 1
  }
  #----------
  current_rank <- 1
  for(buck in 1:bigBuckets_numb){
    tinyBucket <- matrix(rep(NA,tinyBuckets_numb*rows),ncol = rows)
    tBFill <- rep(1,tinyBuckets_numb)
    for(i in stats::na.omit(bigBucket[buck,])){
      secEntry <- workLabs[i,2]
      tinyBucket[secEntry+1,tBFill[secEntry+1]] <- i
      tBFill[secEntry+1] <- tBFill[secEntry+1] + 1
    }
    for(tinybuck in 1:tinyBuckets_numb){
      currentIndices <- stats::na.omit(tinyBucket[tinybuck,])
      if(length(currentIndices)>0){
        ranks[currentIndices] <- current_rank
        current_rank <- current_rank+1
      }
    }
  }
  return(ranks)
}
#' Auxiliary functions
#'
#' \code{getAllAncestors} - Returns all ancestors of \eqn{v} including \eqn{v} itself.
#'
#' @param v A vertex of the tree.
#'
#' @return \code{vectorWithAncs} numeric vector
#' @rdname auxFuncs
#' @export
#' @examples
#' getAllAncestors(tree,v=6)
getAllAncestors <- function(tree,v){
  Ancs <- getAncVec(tree)
  #v is an ancestor of itself
  vectorWithAncs <- v
  #get the parent of the current node as long as it exists
  while(!is.na(Ancs[v])){
    v <- Ancs[v]
    vectorWithAncs <- c(vectorWithAncs, v)
  }
  return (vectorWithAncs)
}
#' Auxiliary functions
#'
#' \code{cPL_inv} - Returns the binary tree that belongs to the input label in an incomplete
#' Newick format.
#'
#' @param label A Colijn-Plazotta label of desired tree, a positive integer.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' cPL_inv(label=6)
cPL_inv <- function(label)
{
  if(label%%1!=0 || label<1) stop("Invalid Colijn Plazzotta label (must be integer > 0).")

  if(label == 1) { # phylo format has no tree with a single node (only single edge)
    return(" ")
  } else { # the vertex has two direct descendants
    i <- ceiling((1+sqrt(8*label-7))/2)-1
    j <- label-1-i*(i-1)/2
  }
  return(paste0("(",cPL_inv(i),",",cPL_inv(j),")"))
}
#' Auxiliary functions
#'
#' \code{maxDepthLeaf} - Returns the maximumy< depth of a leaf in the subtree that
#' is rooted at \eqn{v}.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' maxDepthLeaf(tree,v=6)
maxDepthLeaf<- memoise::memoise(function(tree, v=length(tree$tip.label)+1)
{
  # if v is a leaf return 0
  if(!(v %in% tree$edge[,1])) {
    return (0)
  } else {
    vectorWithChildren <- tree$edge[tree$edge[,1]==v,2]
    maxDepthChild <- sapply(vectorWithChildren, function(x) maxDepthLeaf(tree,x)) # compare children of v
    return (max(maxDepthChild)+1)
  }
})
#' Auxiliary functions
#'
#' \code{get.subtreesize} - Creates a vector that contains at the \eqn{i}-th position
#' the number of leaves in the pending subtree rooted at \eqn{i}.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' get.subtreesize(tree)
get.subtreesize <- function(tree)
{
  n <- length(tree$tip.label)
  Descs <- getDescMatrix(tree)
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  #--------------------------
  nv <- rep(NA,n+tree$Nnode)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      nv[v] <- 1 #if leaf, then nv=1
    }else{
      nv[v] <- sum(nv[stats::na.omit(Descs[v,])])
    }
  }
  return(nv)
}
#' Auxiliary functions
#'
#' \code{getlca} - Returns the name of the lowest common ancestor of the two
#' input vertices \eqn{v} and \eqn{w}.
#'
#' @param w A vertex of the tree.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' getlca(tree,1,2)
getlca <- function(tree,v,w)
{
  lca <- 0
  v_par <- getAllAncestors(tree,v) # vector: node then direct parent to root
  w_par <- getAllAncestors(tree,w)
  for(i in 0:(min(length(v_par),length(w_par))-1))
  {
    if(v_par[length(v_par)-i] == w_par[length(w_par)-i]) {
      lca <- v_par[length(v_par)-i]
    } else {
      break
    }
  }
  if(lca == 0) {
    stop("No common ancestor found.")
  } else {
    return(lca)
  }
}
#' Auxiliary functions
#'
#' \code{we_eth} - Returns the Wedderburn-Etherington number \eqn{we(n)}
#' for a given non-negative integer \eqn{n}.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' we_eth(5)
we_eth <- function(n){
  if(n < 0 || n%%1 != 0)     stop("Input must be non-negative integer.")
  if(n == 0)                 return(0)
  return(treebalance::wedEth[n])
}

#' Auxiliary functions
#'
#' \code{getfurranks} - Returns for each vertex \eqn{i} the Furnas rank of the
#' subtree rooted at \eqn{i}.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' getfurranks(tree)
getfurranks <- function(tree){
  n <- length(tree$tip.label)
  # get the correct node order (from leaves towards root)
  Descs <- getDescMatrix(tree)
  Ancs <- getAncVec(tree)
  depthResults <- getNodesOfDepth(mat=Descs,root=n+1,n=n)
  nodeorder <- rev(stats::na.omit(as.vector(t(depthResults$nodesOfDepth))))
  # compute the number of descendant leaves for each vertex
  subleafnum <- rep(NA,n+tree$Nnode)
  for(v in nodeorder){
    if(is.na(Descs[v,1])){
      subleafnum[v] <- 1 #if leaf, then subleafnum=1
    }else{
      subleafnum[v] <- sum(subleafnum[Descs[v,1:2]])
    }
  }
  # get Wedderburn-Etherington numbers for 1:n
  we <- gmp::as.bigz(treebalance::wedEth[1:n])
  
  # get for each vertex i the rank of the subtree rooted at i
  subranks <- c(rep(gmp::as.bigz(1),n), rep(NA,tree$Nnode))
  for(i in nodeorder) {
    if(i <= n) next
    descs <- Descs[i,1:2]
    nLnR <- subleafnum[descs]
    rLrR <- subranks[descs]
    # ensure left-light rooted ordering
    if(nLnR[1]>nLnR[2]) {
      nLnR <- c(nLnR[2],nLnR[1])
      rLrR <- c(rLrR[2],rLrR[1])
    }
    if(nLnR[1]==nLnR[2] && rLrR[1]>rLrR[2]) {
      nLnR <- c(nLnR[2],nLnR[1])
      rLrR <- c(rLrR[2],rLrR[1])
    }
    nL <- nLnR[1]
    nR <- nLnR[2]
    rL <- rLrR[1]
    rR <- rLrR[2]
    # define auxiliary variable
    h <- gmp::as.bigz(0)
    if(nL > 1) {
      for (j in 1:(nL-1)) {
        h <- gmp::add.bigz(h, gmp::mul.bigz(we[j],we[subleafnum[i]-j]))
      }
    }
    # if left ('light') and right ('heavy') subtree do not have the same number of leaves
    if(nL < nR) {
      subranks[i] <- gmp::add.bigz(h, gmp::as.bigz(rL-1)*we[nR]) + rR - 1 + 1
    }
    # if left ('light') and right ('heavy') subtree have the same number of leaves
    if(nL == nR) {
      subranks[i] <- gmp::sub.bigz(gmp::add.bigz(h, gmp::as.bigz(we[nL]*(we[nL]+1)/2)), (we[nL]-rL+1)*(we[nL]-rL+2)/2) + rR - rL + 1
    }
  }
  
  # return final vector
  return(subranks)
}

#' Auxiliary functions
#'
#' \code{getsubtree} - Returns the pending subtree (in phylo format) that is
#' rooted at the input vertex. If the input vertex is a leaf, the function returns
#' the standard tree for \eqn{n=1} (with 1 edge).
#'
#' @param subroot A vertex of the tree. It is not recommended to use
#' leaves as subroots.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' getsubtree(tree,4)
getsubtree <- function(tree, subroot) {
  numLeaves <- length(tree$tip.label)
  if(subroot%%1 != 0 || subroot>numLeaves+tree$Nnode) stop("subroot must be the number of an inner vertex")

  # if the input vertex is a leaf
  if(sum(tree$edge[,1]==subroot) == 0) {
    warning("The function may not return the desired result if the input subroot is a leaf.")
    return(ape::read.tree(text=paste0("\"(",tree$tip.label[subroot],");\"")))
  }

  # else get all vertices and edges that belong to the wanted subtree
  ancestor    <- subroot
  next.anc    <- NULL
  sub.nodes   <- subroot
  sub.edge    <- NULL
  sub.elength <- NULL
  repeat {
    # for each vertex of the subtree
    next.anc <- NULL
    for(i in 1:length(ancestor)) {
      # get its direct descendants and save the corresponding edges (and their lengths if available)
      pos.descendants <- which(tree$edge[,1]==ancestor[i])
      sub.nodes <- c(sub.nodes,tree$edge[,2][pos.descendants])
      if(length(pos.descendants) >= 1) {
        for(j in 1:length(pos.descendants)) {
          sub.edge <- c(sub.edge, ancestor[i], tree$edge[,2][pos.descendants[j]])
          if("edge.length" %in% names(tree)) {sub.elength <- c(sub.elength, tree$edge.length[pos.descendants[j]])}
        }
      }
      next.anc <- c(next.anc, tree$edge[,2][pos.descendants])
    }
    ancestor <- next.anc
    if(length(ancestor) == 0) break
  }

  # rename the vertices in the subtree and hence the edges
  sub.nodes <- sort(sub.nodes)
  h1 <- which(sub.nodes==subroot)
  h2 <- which(sub.nodes==min(sub.nodes[sub.nodes>numLeaves]))
  sub.nodes[c(h1,h2)] <- sub.nodes[c(h2,h1)] # ensure that the root of the subtree gets the correct number
  for(i in 1:length(sub.edge)) {sub.edge[i] <- which(sub.nodes==sub.edge[i])}

  # get the labels of the tips in the subtree
  sub.tiplabs <- tree$tip.label[seq(1,numLeaves) %in% sub.nodes]

  # get the labels of the inner vertices in the subtree (if they exist)
  if("node.label" %in% names(tree)) {
    sub.nodelabs <- tree$node.label[sub.nodes[sub.nodes>numLeaves]-numLeaves]
  }

  # put everything together
  subtree                <- list()
  class(subtree)         <- "phylo"                               # set class to phylo (mandatory)
  subtree[["Nnode"]]     <- sum(sub.nodes>length(tree$tip.label)) # save number of inner vertices (mandatory)
  subtree[["tip.label"]] <- sub.tiplabs                           # save tip labels (mandatory)
  subtree[["edge"]]      <- matrix(sub.edge, ncol=2, byrow=TRUE)  # save edges (mandatory)
  if("node.label" %in% names(tree))  {subtree[["node.label"]]  <- sub.nodelabs} # save node labels (optional)
  if("edge.length" %in% names(tree)) {subtree[["edge.length"]] <- sub.elength}  # save edge lengths (optional)

  # return subtree
  return(subtree)
}
#' Auxiliary functions
#'
#' \code{is_binary} - Returns TRUE if the input tree is binary and FALSE otherwise.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' is_binary(ape::read.tree(text="((((,),),(,)),(((,),),(,)));"))
is_binary <- function(tree)
{
  #check for errors in input
  if (!inherits(tree, "phylo")) stop("The input tree must be in phylo-format.")
  # if n=1
  n <- length(tree$tip.label)
  if(n==1 && tree$edge[,1]==2) return(TRUE)
  # else: check if each node has exactly 0 or 2 direct descendants
  if(sum(table(tree$edge[,1])!=2)==0) return(TRUE)
}
#' Auxiliary functions
#'
#' \code{is_phylo} - Tests all requirements of the phylo format, and returns TRUE
#' if the tree is correctly formatted, else FALSE with detailed feedback on the
#' features that are not met.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' is_phylo(ape::read.tree(text="((((,),),(,)),(((,),),(,)));"))
is_phylo <- function(tree){
  IS_PHYLO <- TRUE

  #check class
  if (!inherits(tree, "phylo")) {
    IS_PHYLO <- FALSE
    message("The input tree must have class phylo: class(tree) <- \"phylo\"\n")
  }

  #check list elements
  if (sum(c("edge","Nnode","tip.label") %in% names(tree)) != 3) {
    IS_PHYLO <- FALSE
    message("The input tree has to be a list: list(edge,tip.label,Nnode,...).\n")
    return(IS_PHYLO)
  }

  # check tree structure------------------------------------------------------
  n <- length(tree$tip.label)
  m <- tree$Nnode

  if (ncol(tree$edge) != 2) {
    IS_PHYLO <- FALSE
    message("The edge matrix needs to have exactly two columns.\n")
    return(IS_PHYLO)
  }

  # does input fulfill one necessary criterion for a tree?
  if ((m+n-1) != nrow(tree$edge)){
    IS_PHYLO <- FALSE
    message("The input must fulfill |V|-1=|E| to be a tree.\n")
    return(IS_PHYLO)
  }

  # check node labeling-------------------------------------------------------
  # are the nodes labelled with 1,...,|V|?
  if (!identical(seq(1,(m+n)), as.integer(unique(sort(tree$edge))))) {
    IS_PHYLO <- FALSE
    message("Nodes are labeled with other values than 1,...,|V|.\n")
  }

  # are nodes 1,...,n the leaves?
  if (!identical(as.integer(sort(setdiff(tree$edge[,2], tree$edge[,1]))), seq(1,n))) {
    IS_PHYLO <- FALSE
    message("Leaves are labeled with other values than 1,...,n.\n")
  }

  # is (n+1) the root?
  if (sum(tree$edge[,2]==(n+1))!=0){
    IS_PHYLO <- FALSE
    message("Node (n+1) is not the root (n=number of leaves).\n")
  }

  # does each vertex except the root have exactly one direct ancestor?
  if (!identical(seq(1,(n+m))[-(n+1)], as.integer(sort(tree$edge[,2])))) {
    IS_PHYLO <- FALSE
    message("There are vertices (apart from the root) that do not have exactly one direct ancestor.\n")
  }

  if (!IS_PHYLO) return(IS_PHYLO) #break if anything was incorrect

  # is the tree connected?
  nodes <- c(n+1)
  i <- 1
  while(i <= length(nodes)) {
    nodes <- c(nodes, tree$edge[tree$edge[,1]==nodes[i],2])
    i <- i+1
  }
  if (!identical(seq(1,(n+m)), as.integer(sort(nodes)))) {
    IS_PHYLO <- FALSE
    message("Not all vertices are connected.")
    return(IS_PHYLO)
  }

  # check nodes and if tree is binary-----------------------------------------
  # interior nodes must have at least two children
  for(i in (n+1):(m+n)){
    if(sum(tree$edge[,1]==i) < 2) {
      IS_PHYLO <- FALSE
      message(paste("Non-leaves with out-degree <2 are not allowed, e.g.",i,"\n"))
      break
    }
  }

  #check if binary (does not change if phylo)
  if (m != n - 1) {
    message("The tree is not binary.\n")
  }

  # check optional arguments--------------------------------------------------
  if ("edge.length" %in% names(tree)) {
    if (length(tree$edge.length) != nrow(tree$edge)){
      IS_PHYLO <- FALSE
      message("Edge length vector has the wrong length.\n")
    }
    if (!is.numeric(tree$edge.length) || any(is.na(tree$edge.length))){
      IS_PHYLO <- FALSE
      message("Edge length vector contains lengths that are NA or not numeric.")
      return(IS_PHYLO)
    }
    if (any(tree$edge.length<=0)){
      message("Edge length vector has edges of negative or zero length.\n")
    }
  }

  if ("node.label" %in% names(tree)) {
    if (length(tree$node.label) != m){
      IS_PHYLO <- FALSE
      message("Length of node label vector not appropriate.\n")
    }
  }

  if ("root.edge" %in% names(tree)) {
    if (length(tree$root.edge)!=1 || !is.numeric(tree$root.edge)){
      IS_PHYLO <- FALSE
      message("Root edge length not appropriate.\n")
    }
  }

  #---------------------------------------------------------------------------
  return (IS_PHYLO)
}
#' Auxiliary functions
#'
#' \code{tree_decomposition} - Returns a list of length two, which
#' contains the two pending subtrees that are rooted at the children of the root
#' of the input tree. The
#' smaller one (according to the number of leaves) is stated first.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' tree_decomposition(ape::read.tree(text="((((,),),(,)),(((,),),(,)));"))
tree_decomposition <- function(tree)
{
  # number of leaves
  num.leaves <- length(tree$tip.label)

  # non-binary trees
  if(!is_binary(tree)) stop("The input tree is not binary.")

  # trees with one leaf
  if(num.leaves == 1) stop("A decomposition of the input tree is not possible since it contains just one leaf.")

  # trees with two leaves
  if(num.leaves == 2)
  {
    first.subtree  <- ape::read.tree(text="();")
    second.subtree <- ape::read.tree(text="();")
    maxpensub <- list(first.subtree, second.subtree) #contains maximal pending subtrees
    return(maxpensub)
  }

  # trees with more than two leaves
  # subtreelist: list of all subtrees with at least two leaves
  subtreelist <- ape::subtrees(tree)
  first.subtree <- subtreelist[[2]]
  # if first subtree has num.leaves-1 leaves, the second subtree consists of one node
  if(first.subtree$Ntip == num.leaves-1)
  {
    second.subtree <- ape::read.tree(text="();")
  }
  # otherwise search for second subtree in subtreelist
  else
  {
    second.subtree <- subtreelist[[first.subtree$Nnode + 2]]
  }
  # order first and second subtree according to their number of leaves (smaller one is stated first)
  if(length(first.subtree$tip.label) <= length(second.subtree$tip.label))
  {
    maxpensub <- list(first.subtree, second.subtree) # contains maximal pending subtrees
  }
  else
  {
    maxpensub <- list(second.subtree, first.subtree) # contains maximal pending subtrees
  }
  return(maxpensub)
}

#' Auxiliary functions
#'
#' \code{tree_merge} - Returns a rooted tree \eqn{T} in phylo
#' format, which contains the input trees \eqn{tree1} and \eqn{tree2} as
#' "left" and "right" maximal pending subtrees.
#'
#' @param tree1 A rooted tree in phylo format.
#' @param tree2 A rooted tree in phylo format.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' treeA <- ape::read.tree(text="(((,),),(,));")
#' treeB <- ape::read.tree(text="((,),);")
#' tree_merge(treeA, treeB)
tree_merge <- function(tree1, tree2)
{
  # tree with 1+1=2 leaves
  if(length(tree1$tip.label)==1 & length(tree2$tip.label)==1)
  {
    tree <- ape::read.tree(text="(,);")
    return(tree)
  }
  # tree with 1+x leaves
  if(length(tree1$tip.label)==1 & length(tree2$tip.label)!=1)
  {
    r1 <- length(tree1$tip.label) + 1
    r2 <- length(tree2$tip.label) + 1
    tree1$root.edge <- 1
    tree2$root.edge <- 1
    # tree1 and tree2 are bound together
    tree <- ape::bind.tree(tree1, tree2)
    tree$root.edge <- NULL
    return(tree)
  }
  # tree with x+1 leaves
  if(length(tree1$tip.label)!=1 & length(tree2$tip.label)==1)
  {
    r1 <- length(tree1$tip.label) + 1
    r2 <- length(tree2$tip.label) + 1
    #tree1 and tree2 are bound together
    tree <- ape::bind.tree(tree1, tree2, where=r1, position=r2)
    tree$root.edge <- NULL
    return(tree)
  }
  # tree with x+y leaves
  else
  {
    r1 <- length(tree1$tip.label) + 1
    r2 <- length(tree2$tip.label) + 1
    tree2$root.edge <- 1
    # tree1 and tree2 are bound together
    tree <- ape::bind.tree(tree1, tree2, where=r1, position=r2)
    tree$root.edge <- NULL
    return(tree)
  }
}

#' Auxiliary functions
#'
#' \code{treenumber} - Returns the unique tree number \eqn{tn(T)} of the given tree.
#' \eqn{tn(T)} is the rank of the tree \eqn{T} among all
#' rooted binary trees in the left-light rooted ordering. It can
#' be calculated as follows: \deqn{tn(T)=F(T) + \sum_{i=1}^{n-1} we(i)}{
#' tn(T)=F(T) + \sum we(i) over 1=i<=n-1}
#' in which \eqn{n} is the number of leaves in \eqn{T}, \eqn{F(T)} is the Furnas
#' rank of \eqn{T}, i.e. the rank of \eqn{T} in the left-light rooted ordering
#' of all rooted binary trees with \eqn{n} leaves, and \eqn{we(i)} is the
#' Wedderburn-Etherington number of \eqn{i}.
#' The concept of assigning each rooted binary tree a unique tree number allows
#' to store many trees with minimal storage use.
#' For \eqn{n=1} the function returns \eqn{tn(T)=1} and a warning.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' treenumber(ape::read.tree(text="((((,),),(,)),(((,),),(,)));"))
treenumber <- function(tree)
{
  # number of leaves
  num.leaves <- length(tree$tip.label)

  if(num.leaves == 1)
  {
    warning("The function might not deliver accurate results for n=1.")
    return(1)
  }

  # calculate the tree number of the given tree by adding the Furnas rank to
  # the sum of Wedderburn-Etherington numbers from 1 to num.leaves-1
  treenum <- sum(sapply(seq(1,num.leaves-1),we_eth)) + furnasI(tree)

  # return the tree number
  return(treenum)
}

#' Auxiliary functions
#'
#' \code{treenumber_inv} - Returns the unique tree (in phylo format) for
#' the given tree number.
#'
#' @param treenum An integer denoting the tree number of the sought tree.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' treenumber_inv(192)
treenumber_inv <- function(treenum)
{
  # check for errors in input
  if(treenum%%1!=0 | treenum<1) stop("Tree cannot be calculated, because tree number is not valid.")

  # initial conditions
  if(treenum == 1) {return(ape::read.tree(text="();"))}
  if(treenum == 2) {return(ape::read.tree(text="(,);"))}

  # calculate the number of leaves n of the tree
  we_sum <- 0
  nT <- 1
  while(we_sum < treenum)
  {
    we_sum <- we_sum + we_eth(nT)
    nT <- nT + 1
  }
  nT <- nT - 1
  we_sum <- we_sum - we_eth(nT)

  # calculate the rank of the tree among all rooted binary trees with n leaves
  rankT <- treenum - we_sum

  # calculate the tree and return
  return(furnasI_inv(rank=rankT, n=nT))
}

#' Auxiliary functions
#'
#' \code{auxE_l_X} - Returns the sum of all products of l different values in X.
#'
#' @param subX integer >=1, size of the subsets of X.
#' @param Xset Vector (multiset) of numeric values.
#'
#' @rdname auxFuncs
#' @export
#' @examples
#' auxE_l_X(subX=3,Xset=c(1,1,2,2))
auxE_l_X <- function(subX,Xset){
  if(subX>length(Xset)) return(0)
  if(subX==length(Xset)) return(prod(Xset))
  if(subX==1) return(sum(Xset))
  P_i <- sapply(1:subX, function(x) sum(Xset^x))
  mat <- matrix(0, nrow = subX, ncol = subX)
  for(row in 1:subX){
    mat[row,1:row] <- rev(P_i[1:row])
  }
  for(row in 1:(subX-1)){
    mat[row,row+1] <- row
  }
  E_result <- det(mat)/factorial(subX)
  return(E_result)
}
