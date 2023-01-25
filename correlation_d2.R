intpsc2 <- readRDS("/Users/qinluo/Downloads/intpsc2.RDS")
epi_pgc <- readRDS("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/epi_pgc.rds")
torev = subset(epi_pgc,subset=batch=="Tore_bioRxiv_2020") 
torev = subset(torev,subset=lineage!="EB_d2")
torev = subset(torev,subset=lineage!="EB_d4")
torev = subset(torev,subset=lineage!="ESC_SL")
torev = subset(torev,subset=lineage!="EpiSC")
torev <- subset(torev, subset = nCount_RNA < 75000  & percent.mt < 5)

setwd("/Users/qinluo/Desktop/revision")
library(Seurat)
d2 = subset(torev,subset = lineage=="EpiLC_d2_A")
DefaultAssay(intpsc2)= "RNA"
genes = intersect(rownames(d2),rownames(intpsc2))

   
cc = AggregateExpression(intpsc2,features = genes,
                         group.by = "id",
                         slot = "counts",assays = "RNA")

cc = cc$RNA
cc = log2(cc+1)
d2 = as.data.frame(d2@assays$RNA@counts)
d2 = d2[genes,]
d2=log(d2+1)
dim(d2)
#cor = data.frame(cor=0, clus=0)
cor = data.frame(cor1=0, cor2=0,cor3=0)

for (i in c(1:263)){
    test1=cor.test(d2[,i],
                   cc[,1])
    cor1  =test1$estimate
    test2=cor.test(d2[,i],
                   cc[,2])
    cor2  =test2$estimate
    test3=cor.test(d2[,i],
                   cc[,3])
    cor3  =test3$estimate
    
   # cor[i,1] = max(cor1,cor2,cor3)
    #cor[i,2] = which.max(c(cor1,cor2,cor3))
    
    cor[i,] = c(cor1,cor2,cor3)
}
#cor = rbind(cor,c(0,3))
#cor$clus= as.factor(cor$clus)
pheatmap::pheatmap(cor,cluster_cols = F,
                   cluster_rows = T,
                   scale = "row",show_rownames = F)

count1 = as.data.frame(intpsc2@assays$RNA@counts)
count1 = count1[genes,]
cor2 = data.frame(row.names = 1:528,c=rep(0,528))

for (i in c(1:263)){
    r1 = cbind(d2[,i],count1)
    cor1 = cor(r1)
    cor1 = cor1[1,-1]
    cor2[,i]=cor1
    print(i)
}



rownames(cor2) = colnames(intpsc2)
cor2 = t(cor2)
annotation_col = data.frame(
    CellType = intpsc2$id)
pheatmap::pheatmap(cor2,
 annotation_col  =annotation_col ,
 show_colnames = F,show_rownames = F)
