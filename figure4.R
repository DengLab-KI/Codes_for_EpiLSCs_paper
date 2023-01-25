library(Seurat)
emb_new <- readRDS("/Users/qinluo/OneDrive - KI.SE/Mac/Desktop/velocity/emb_new.rds")
formativeCell = subset(emb_new,subset=batch=="formative")
formativeCell = subset(formativeCell,subset = lineage2!="fPSC")
formativeCell$lineage2[
  which(formativeCell$lineage2=="INTPSC")]="dyEpiLC"
DefaultAssay(formativeCell)="RNA"
clustercol=c("#eeac99", "#8ca3a3", "#77AADD", "orange", 
             "#99CCBB", "#117755","#777711",  
             "#ffef96",  "#c1946a", "#82b74b", "#034f84", 
             "#b8a9c9", "#F5DEB3", "#7e4a35","#778822") 

formativeCell$id = formativeCell$seurat_clusters
formativeCell$id = as.numeric(formativeCell$id)
formativeCell$id [which(formativeCell$id==2)]="C1"
formativeCell$id [which(formativeCell$id==6)]="C2"
formativeCell$id [which(formativeCell$id==3)]="C3"
formativeCell$id [which(formativeCell$id==1)]="C4"
formativeCell$id [which(formativeCell$id==5)] ="C5"
formativeCell$id [which(formativeCell$id==4)] ="C6"

all.genes <- rownames(formativeCell)
#formativeCell = NormalizeData(formativeCell)

deg = FindAllMarkers(emb_new, #test.use="MAST",
                     #only.pos=T,
                     #logfc.threshold=1,
                     #group.by = "id",
                     assay = "integrated",slot="data")
deg2 = deg[which(abs(deg$avg_log2FC)> 0.5849625),]

genes = unique(deg$gene)
genes2 = unique(deg$gene)[1:3000]
#formativeCell <- FindVariableFeatures(formativeCell,
#                                     selection.method = "vst",
#                                     nfeatures = 10000)

formativeCell <- ScaleData(formativeCell, 
                           features = genes,
                           vars.to.regress =c("lineage","lineage2" ))

formativeCell <- RunPCA(formativeCell, 
                        features = genes2,
                        npcs = 50,assay = "RNA")
ElbowPlot(formativeCell,ndims = 50)
library(harmony)
#formativeCell<- RunHarmony(formativeCell,
#                           group.by.vars = c("lineage","lineage2"))
#ElbowPlot(formativeCell,reduction = "harmony")
formativeCell<- FindNeighbors(formativeCell, 
                              #reduction = "harmony",
                              dims = 1:10)
formativeCell <- FindClusters(formativeCell,
                              resolution = 0.5)
formativeCell <- RunUMAP(formativeCell, 
                         #reduction = "harmony",
                         dims = 1:20,
                         metric = "cosine",
                         spred=1,min.dist = 0.1,n.neighbors = 30)

DimPlot(formativeCell , reduction = "umap",
        group.by = "id",cols =  c("#6b5b95","#feb236","#82b74b",
                                         "#dac292","#92a8d1","#034f84"))
DimPlot(formativeCell , reduction = "umap",
        group.by = "lineage2",
        cols = clustercol[10:15])
setwd("/users/qinluo/Downloads")
pdf("formativecell.pdf",width=24,height=24)
FeaturePlot(formativeCell,features = 
              c("Esrrb","Klf4","Klf5","Tdh",
                "Stat3","Fgf4","Dnmt3l","Dnmt3a",
                "Otx2","Pou3f1","Eomes","Hes1","T","Krt8"),
            split.by = "lineage2")
dev.off()

FeaturePlot(formativeCell,features = 
              c("Esrrb","Klf4","Stat3",
                "Dnmt3b","Otx2","T","Krt8"))
