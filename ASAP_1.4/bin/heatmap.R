#! R
args<-commandArgs(TRUE)

#install.packages("pheatmap")
library("pheatmap")

data<-read.csv(args[1],header=T,sep="\t",row.names=1)
x<-as.matrix(data)
x<-sqrt(x)
colfunc <- colorRampPalette(c("white", "blue"))
pdf("heatmap.pdf", width=(length(x[1,])*0.5+15), height=(length(x[,1])*0.5)+15)
pheatmap(x, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="average", color=colfunc(15), fontsize_row=(10/length(x[,1]))+20, fontsize_col=(10/length(x[,1]))+20, )
dev.off()
