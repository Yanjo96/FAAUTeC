#!/usr/bin/R
library(kdetrees)
library(ape)
library(ggplot2)
library(plotly)

ks <- seq(0,2,0.1)

tree_data <- read.table("SUMMARY/pylogeny_hypoTreeShortestDistUnconstTree.tre",sep=" ", header = FALSE, stringsAsFactors = FALSE)
multitrees <- ape::read.tree(text=tree_data$V2, tree.names = tree_data$V1)
multitrees <- multi2di(multitrees)

klen <- as.numeric(length(ks))
tlen <- as.numeric(length(tree_data$V1))
mat_geo_pylogeny <- matrix(rep(0,klen*tlen),nrow=tlen,ncol=klen, dimnames = list(tree_data$V1, ks))
mat_dis_pylogeny <- matrix(rep(0,klen*tlen),nrow=tlen,ncol=klen, dimnames =list(tree_data$V1, ks))
    
# See Advanced options (outlier detection [k] & distance method [distance])
for(k in ks){
  result <- kdetrees(multitrees, k=k, distance="geodesic")
  for(out in names(result$outliers)){
    mat_geo_pylogeny[out,as.character(k)] <- 1
  }
  result <- kdetrees(multitrees, k=k, distance="dissimilarity")
  for(out in names(result$outliers)){
    mat_dis_pylogeny[out,as.character(k)] <- 1
  }
}

tree_data <- read.table("SUMMARY/iqtree_hypoTreeShortestDistUnconstTree.tre",sep=" ", header = FALSE, stringsAsFactors = FALSE)
multitrees <- ape::read.tree(text=tree_data$V2, tree.names = tree_data$V1)
multitrees <- multi2di(multitrees)

klen <- as.numeric(length(ks))
tlen <- as.numeric(length(tree_data$V1))
mat_geo_iqtree <- matrix(rep(0,klen*tlen),nrow=tlen,ncol=klen, dimnames = list(tree_data$V1, ks))
mat_dis_iqtree <- matrix(rep(0,klen*tlen),nrow=tlen,ncol=klen, dimnames = list(tree_data$V1, ks))

# See Advanced options (outlier detection [k] & distance method [distance])
for(k in ks){
  result <- kdetrees(multitrees, k=k, distance="geodesic")
  for(out in names(result$outliers)){
    mat_geo_iqtree[out,as.character(k)] <- 1
  }
  result <- kdetrees(multitrees, k=k, distance="dissimilarity")
  for(out in names(result$outliers)){
    mat_dis_iqtree[out,as.character(k)] <- 1
  }
}

mat_dis_iqtree <- reshape2::melt(mat_dis_iqtree)
mat_geo_iqtree <- reshape2::melt(mat_geo_iqtree)

mat_dis_pylogeny <- reshape2::melt(mat_dis_pylogeny)
mat_geo_pylogeny <- reshape2::melt(mat_geo_pylogeny)

ggplot(mat_dis_iqtree, aes(Var2, Var1)) + geom_tile(aes(fill = value), colour = "white", show.legend=F) + scale_fill_gradient(low = "darkgray", high = "black") + ggtitle("IQTree dissimilarity") +
  xlab("k") + ylab("Best tree for gene")
ggplot(mat_geo_iqtree, aes(Var2, Var1)) + geom_tile(aes(fill = value), colour = "white", show.legend=F) + scale_fill_gradient(low = "darkgray", high = "black") + ggtitle("IQTree geodesic") +
  xlab("k") + ylab("Best tree for gene")
ggplot(mat_dis_pylogeny, aes(Var2, Var1)) + geom_tile(aes(fill = value), colour = "white", show.legend=F) + scale_fill_gradient(low = "darkgray", high = "black") + ggtitle("Pylogeny dissimilarity") +
  xlab("k") + ylab("Best tree for gene")
ggplot(mat_geo_pylogeny, aes(Var2, Var1)) + geom_tile(aes(fill = value), colour = "white", show.legend=F) + scale_fill_gradient(low = "darkgray", high = "black")+ ggtitle("IQTree geodesic") +
  xlab("k") + ylab("Best tree for gene")