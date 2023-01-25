intpsc2 <- readRDS("/Users/qinluo/Downloads/intpsc2.RDS")

epi_pgc <- readRDS("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/epi_pgc.rds")

torev = subset(epi_pgc,subset=batch=="Tore_bioRxiv_2020") 
torev = subset(torev,subset=lineage!="EB_d2")
torev = subset(torev,subset=lineage!="EB_d4")
torev = subset(torev,subset=lineage!="ESC_SL")
torev = subset(torev,subset=lineage!="EpiSC")



torev <- subset(torev, subset = nCount_RNA < 75000  & percent.mt < 5)
VlnPlot(torev,features = c("nCount_RNA","nCount_Feature","percent.mt"))
formative <- readRDS("/Users/qinluo/Downloads/intpsc_wgcna_3clusters/formative.rds")
EpiLSC = subset(formative ,subset=lineage2=="dyEpiLC")

#saveRDS(EpiLSC,"./Desktop/EpiLSC.rds")
#saveRDS(torev,"./Desktop/torev.rds")

suppressMessages(library(Seurat))
suppressMessages(library(slingshot))
suppressMessages(library(splatter))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tradeSeq))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(harmony))
suppressMessages(library(scater))
suppressMessages(library(mclust, quietly = TRUE))
suppressMessages(library(MLmetrics))
suppressMessages(library(scStalt))

clustercol=c("#eeac99", "#77AADD", "#c1946a", "#82b74b", "#034f84",  "#b8a9c9", "#F5DEB3", "#7e4a35") 

clustercol=c("#eeac99", "#8ca3a3", "#77AADD", "orange", 
             "#99CCBB", "#117755","#777711",  
            "#ffef96",  "#c1946a", "#82b74b", "#034f84", 
            "#b8a9c9", "#F5DEB3", "#7e4a35","#778822") 
            
torev <- NormalizeData(torev, normalization.method = "LogNormalize")
torev <- FindVariableFeatures(torev,selection.method = "vst",
                              nfeatures = 3000)
torev = ScaleData(torev)
torev = RunPCA(torev)
DimPlot(torev,group.by = "lineage",reduction = "pca")
torev = FindNeighbors(torev,dims = 1:20)
torev = FindClusters (torev)
DimPlot(torev,reduction = "pca")

data = torev@reductions$pca@cell.embeddings[,1:2]
#pca <- prcomp(t(torev@assays$RNA@data[VariableFeatures(object = torev),]), 
#              scale. = FALSE)
#data =  pca$x[,1:2]

colnames(data) = c("PC_1" ,"PC_2")

#torev@reductions$pca@cell.embeddings = data
DimPlot(torev,reduction = "pca",
        group.by = "lineage")

#cl1 <- Mclust( pca$x[,1:2])
#cl1 = cl1$classification
#torev[["clus"]] = cl1
#DimPlot(torev,reduction = "pca",
#        group.by = "clus",label = T)


#torev[["mclust"]] = cl1
torev.sce = as.SingleCellExperiment(torev)
reducedDims(torev.sce) = SimpleList(PCA = data)
colData(torev.sce)$GMM = torev$lineage

#VlnPlot(torev,features = "nCount_RNA",
#        group.by = "clus")
#torev = subset(torev,subset=clus!="5")
torev.sce = as.SingleCellExperiment(torev)
reducedDims(torev.sce) = SimpleList(PCA = torev@reductions$pca@cell.embeddings)
colData(torev.sce)$GMM = torev$lineage
torev.sce = slingshot(torev.sce, 
                      clusterLabels = 'GMM', 
                      reducedDim = 'PCA',
                      start.clus="ES_2i_A",
                      end.clus="EpiLC_d3")
pseudotime = slingPseudotime(torev.sce,na=F)
cellWeights = slingCurveWeights(torev.sce)
#data =  torev@reductions$pca@cell.embeddings
#lin1 <- getLineages(data, torev$clus)
#plot(data, col = brewer.pal(4,"Set1")[cl1], asp = 1, pch = 16)
#lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

torev$time = pseudotime
FeaturePlot(torev,features = "time")

library(scStalt)
ref = doKernel(count=torev@assays$RNA@counts,
               pseudotime = pseudotime,
               cellWeights = cellWeights, 
               nGenes = 3000, nKnots = 6,
               nPoints = 200, fdr = 0.05)

DefaultAssay(EpiLSC) = "RNA"
EpiLSC = NormalizeData(EpiLSC)
EpiLSC = ScaleData(EpiLSC)
#query_b2 = EpiLSC@assays$RNA@data[rownames(ref),]
#query_b2 = as.data.frame(query_b2)
#query_id = stageAssign(count = query_b2,kernel = log2(ref+1))
#query_id = as.data.frame(query_id)
#query_id$lineage_cosine = as.numeric(query_id$lineage_cosine)
#query_id$score_cosine = as.numeric(query_id$score_cosine)

count1= as.data.frame(torev@assays$RNA@counts)
count2=as.data.frame(EpiLSC@assays$RNA@counts)
genes = intersect(rownames(count1),rownames(count2))
count1 = count1[genes,]
count2 = count2[genes,]
bm = merge(torev[genes,],EpiLSC[genes,])

bm = NormalizeData(bm,verbose = F)
bm = FindVariableFeatures(bm,verbose = F)
bm = ScaleData(bm ,verbose = F)
bm = RunPCA(bm ,verbose = F)
DimPlot(bm,group.by = "lineage2")

g = c(grep("mt-",rownames(ref)),grep("Rsp",rownames(ref)),
      grep("Rpl",rownames(ref)),
      grep("Gm",rownames(ref)))
ref = ref[-g,]


all = merge(torev[genes,],EpiLSC[genes,])
all = FindVariableFeatures(all,nfeatures = 3000)

count1= as.data.frame(torev@assays$RNA@counts)
count2=as.data.frame(EpiLSC@assays$RNA@counts)
genes = intersect(rownames(count1),rownames(count2))
count1 = count1[genes,]
count2 = count2[genes,]

bm_stalt_corrected =staltIntegration2(count1,
                                     count2,  
                                     nPoints = 200,
                                     refk=ref,
                                     ref=ref)


bm_stalt_corrected = CreateSeuratObject(bm_stalt_corrected)
bm_stalt_corrected$batch = bm$batch
bm_stalt_corrected$col = bm$lineage2
VlnPlot(bm_stalt_corrected,
        features = c("nCount_RNA","nFeature_RNA"),group.by = "batch")
VlnPlot(bm_stalt_corrected,
        features = c("nCount_RNA","nFeature_RNA"),group.by = "col")

#bm_stalt_corrected$step = c(b1_sub$Step,b2_sub$Step)
#bm_stalt_corrected$lineage =c(b1_id$lineage_cosine,b2_id$lineage_cosine)
bm_stalt_corrected = NormalizeData(bm_stalt_corrected)
bm_stalt_corrected = FindVariableFeatures(bm_stalt_corrected)
#bm_stalt_corrected@assays$RNA@data = bm_stalt_corrected@assays$RNA@counts
bm_stalt_corrected = ScaleData(bm_stalt_corrected )
bm_stalt_corrected = RunPCA(bm_stalt_corrected,
             features = rownames(bm_stalt_corrected) ) 
bm_stalt_corrected[["id"]] = as.vector(bm_stalt_corrected$col)
bm_stalt_corrected$id[which(bm_stalt_corrected$id=="dyEpiLC")]=as.vector(intpsc2$id)
DimPlot(bm_stalt_corrected,group.by = "id")
bm_stalt_corrected = RunUMAP(bm_stalt_corrected,reduction = "pca",
                             dims = 1:10)
DimPlot(bm_stalt_corrected,group.by = "id")

DimPlot(bm_stalt_corrected,group.by = "id",
        reduction = "pca",
       split.by="batch")

DimPlot(bm_stalt_corrected,group.by = "id",
        reduction = "pca",cols = c("#82b74b","#6b5b95",
            "#dac292",clustercol[c(1,2,11,8)]),
        shuffle = T)+
  theme_void()
saveRDS(bm_stalt_corrected,"/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Desktop/0909/0306/velocity/bm_stalt_corrected.rds")

############################################################
############################################################
#2 prepare velo files
setwd("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Desktop/0909/0306/velocity")
library(SeuratDisk)
library(anndata)

c1 = as.data.frame(torev@assays$RNA@counts)
c2 = as.data.frame(EpiLSC@assays$RNA@counts)
genes = intersect(rownames(c1),rownames(c2))
call = cbind(c1[genes,],c2[genes,])
all_new = CreateSeuratObject(call)
all_new$id = bm_stalt_corrected$id
all_new$batch = bm_stalt_corrected$batch
all_new = NormalizeData(all_new)
all_new = ScaleData(all_new)
all_new = FindVariableFeatures(all_new)
all_new = RunPCA(all_new)
all_new@reductions$pca@cell.embeddings = 
  bm_stalt_corrected@reductions$pca@cell.embeddings
#VariableFeatures(all_new) = rownames(ref)

dy_scvelo = subset(all_new,subset=batch=="formative")
tore_scvelo = subset(all_new,subset=batch!="formative")
DimPlot(dy_scvelo ,reduction = "pca",group.by = "id")
SaveH5Seurat(dy_scvelo , 
             filename = "EpiLSC.h5Seurat")
Convert("EpiLSC.h5Seurat", dest = "h5ad")

SaveH5Seurat(tore_scvelo , 
             filename = "tore_new.h5Seurat")
Convert("tore_new.h5Seurat", dest = "h5ad")


############################################################
############################################################
#3 prepare e/nmat
######
setwd("/Users/qinluo/Downloads")
library(velocyto.R)
library(sccore)
ldat1 <- read.loom.matrices("ES.loom")
ldat2 <- read.loom.matrices("d1.loom")
ldat3 <- read.loom.matrices("d2.loom")
ldat4 <- read.loom.matrices("d3.loom")
ldat5 <- read.loom.matrices("EpiSC.loom")

#source("final_code_chir_compare_...")

emat1 <- ldat1$spliced
emat2 <- ldat2$spliced
emat3 <- ldat3$spliced
emat4 <- ldat4$spliced
emat5 <- ldat5$spliced
#library(pagoda2)
#r <- Pagoda2$new(emat,modelType='plain',trim=10,log.scale=T)
#r$adjustVariance(plot=T,do.par=T,gam.k=10)
#r$calculatePcaReduction(nPcs=100,n.odgenes=3e3,maxit=300)
#r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine');
#r$getKnnClusters(method=multilevel.community,type='PCA',name='multilevel')
#r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T)

#par(mfrow=c(1,2))
#r$plotEmbedding(type='PCA',embeddingType='tSNE',
#show.legend=F,mark.clusters=T,min.group.size=10,
#shuffle.colors=F,mark.cluster.cex=1,alpha=0.3,
#main='cell clusters')

emat1 <- ldat1$spliced; 
nmat1 <- ldat1$unspliced
emat2 <- ldat2$spliced; 
nmat2 <- ldat2$unspliced
emat3 <- ldat3$spliced; 
nmat3 <- ldat3$unspliced
emat4 <- ldat4$spliced; 
nmat4 <- ldat4$unspliced
nmat5 <- ldat5$unspliced
genes <- rownames(emat1)

colnames(emat1) = gsub( "ES:","",colnames(emat1))
colnames(emat2) = gsub( "d1:","",colnames(emat2))
colnames(emat3) = gsub( "d2:","",colnames(emat3))
colnames(emat4) = gsub( "d3:","",colnames(emat4))
colnames(emat5) = gsub( "EpiSC:","",colnames(emat5))
colnames(nmat1) = gsub( "ES:","",colnames(nmat1))
colnames(nmat2) = gsub( "d1:","",colnames(nmat2))
colnames(nmat3) = gsub( "d2:","",colnames(nmat3))
colnames(nmat4) = gsub( "d3:","",colnames(nmat4))
colnames(nmat5) = gsub( "EpiSC:","",colnames(nmat5))

colnames(emat1) = gsub( "x","",colnames(emat1))
colnames(emat2) = gsub( "x","",colnames(emat2))
colnames(emat3) = gsub( "x","",colnames(emat3))
colnames(emat4) = gsub( "x","",colnames(emat4))
colnames(emat5) = gsub( "x","",colnames(emat5))
colnames(nmat1) = gsub( "x","",colnames(nmat1))
colnames(nmat2) = gsub( "x","",colnames(nmat2))
colnames(nmat3) = gsub( "x","",colnames(nmat3))
colnames(nmat4) = gsub( "x","",colnames(nmat4))
colnames(nmat5) = gsub( "x","",colnames(nmat5))
colnames(emat1) = paste(colnames(emat1),"1",sep=".")
colnames(nmat1) = paste(colnames(nmat1),"1",sep=".")
colnames(emat2) = paste(colnames(emat2),"3",sep=".")
colnames(nmat2) = paste(colnames(nmat2),"3",sep=".")
colnames(emat3) = paste(colnames(emat3),"2",sep=".")
colnames(nmat3) = paste(colnames(nmat3),"2",sep=".")
colnames(emat4) = paste(colnames(emat4),"8",sep=".")
colnames(nmat4) = paste(colnames(nmat4),"8",sep=".")
colnames(emat5) = paste(colnames(emat5),"5",sep=".")
colnames(nmat5) = paste(colnames(nmat5),"5",sep=".")

emat_invitro = cbind(emat1,emat2,emat3,emat4,emat5)
nmat_invitro = cbind(nmat1,nmat2,nmat3,nmat4,nmat5)
emat_invitro = emat_invitro[!duplicated(rownames(emat_invitro)),]
nmat_invitro = nmat_invitro[!duplicated(rownames(nmat_invitro)),]
#emat <- emat[,rownames(r$counts)];
#nmat <- nmat[,rownames(r$counts)]; 
emat <- readRDS("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Desktop/0909/0306/velocity/emat.rds")
nmat <- readRDS("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Desktop/0909/0306/velocity/nmat.rds")

gene = intersect(rownames(emat),rownames(emat_invitro))

emat <- emat[gene,]
nmat <- nmat[gene,] 
emat_invitro <- emat_invitro[gene,]
nmat_invitro <- nmat_invitro[gene,] 

emat_all = cbind(emat,emat_invitro)
nmat_all = cbind(nmat,nmat_invitro)

emat_all = emat_all[,intersect(colnames(bm_stalt_corrected),colnames(emat_all))]
nmat_all = nmat_all[,intersect(colnames(bm_stalt_corrected),colnames(nmat_all))]

saveRDS(emat_all,"emat_all.rds")
saveRDS(nmat_all,"nmat_all.rds")
############################################################
############################################################
#4 velo