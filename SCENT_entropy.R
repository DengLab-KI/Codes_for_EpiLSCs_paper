library('biomaRt')
library(R.matlab)
library(ggplot2)
library(SCENT)
library(LandSCENT)
library(Seurat)
hmo <- readMat("~/Downloads/human_mouse_orthologue_pairs_20171018.mat")
hmo = data.frame(unlist(hmo))
hmo = cbind(hmo[1:18501,],hmo[18502:37002,])
colnames(hmo)=c("hs","mm")

library(SCENT);
data(net17Jan16)

emb$lineage[which(emb$lineage=="E7.0")]="ecto_early"

#need to use log normalized data
merge = emb
merge$lineage[which(merge$lineage=="E7.0")]="ecto_early"

#merge = SCTransform(merge,vars.to.regress = c("batch","split"))
rt = merge@assays$integrated@data
#rt = merge@assays$RNA@counts
#rt = log2(rt+1.1)
rt = as.data.frame(rt)
genes = rownames(rt)
genes= intersect(genes,hmo[,2])
rt = rt[ genes,]
genes = rownames(rt)
genes_match = match(genes,hmo[,2])
rt =cbind(rt,id= hmo[ genes_match,1])
rt = na.omit(rt)
rt = rt[!duplicated(rt$id),]
rownames(rt) = rt$id
rt = rt[,-dim(rt)[2]]





mart <- useDataset("hsapiens_gene_ensembl",
                   useMart("ensembl"))

genes = rownames(rt)
id = getBM(attributes = c('entrezgene_id', 'hgnc_symbol'),
           filters = 'hgnc_symbol',
           values = genes,
           mart = mart)


id =  id[complete.cases(id), ]
id = id[!duplicated(id$hgnc_symbol),]
id = id[!duplicated(id$entrezgene_id),]

x = match(id$hgnc_symbol, genes)
rt = rt[x,]

identical(rownames(rt),id$hgnc_symbol)
rt$id = id$entrezgene_id

rownames(rt ) = rt$id
rt = rt[,-dim(rt)[2]]

rt = as.matrix(rt)

#ccat.v = CompCCAT(exp.m  = rt, ppiA = net17Jan16.m)
library(SCENT)
library(LandSCENT)
Integration.l <- DoIntegPPI(exp.m = rt, ppiA.m = net17Jan16.m)
Integration.l = CCAT(  Integration.l = Integration.l)
Integration.l$coordinates = merge@reductions$umap@cell.embeddings



potency_state.v <- 0



#merge
dt = data.frame(lineage = merge$lineage2,
                cluster=merge$batch,
                entropy=Integration.l$CCAT ,umap1=merge@reductions$umap@cell.embeddings[,1],
                umap2=merge@reductions$umap@cell.embeddings[,2],
                id=merge$seurat_clusters)

saveRDS(dt,"entrop.rds")
dt$category <- cut(dt$entropy, 
                   breaks=c(0.1,0.4, 0.45, 0.46, 0.47 ,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.6), 
                   labels=c(1:13))

dt$formative = 0
dt$formative[which(dt$id==3|dt$id==2|dt$id==6)]=1

library(viridis)
p2=ggplot(dt)+geom_point(aes(x=umap1,y=umap2,col=category))+
  scale_color_manual(values =  viridis(13))


ggplot(dt)+geom_violin(aes(x=lineage,y=entropy),trim = T,scale="width")+
  geom_jitter(aes(x=lineage,y=entropy,col=cluster))+
  scale_color_manual(values = Myrainbow)+facet_wrap(~formative)


dt$lineage = factor(dt$lineage,levels = c("ES 2i", "Day1 EpiLC",
                                          "Day2 EpiLC",
                                          "EpiSC","E3.5","E4.5","E5.5",
                                          "E6.5","ecto_early","streak_pre","PS_early","INTPSC",
                                          "rosetteStage","FS","fPSC"))
pdf("fs_entropy.pdf")
ggplot(dt,aes(x=lineage,y=entropy,fill=lineage))+
  facet_wrap(~formative)+
  geom_violin(trim = T,scale="width")+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = .6,
               colour = "black")+
  scale_fill_manual(values = Myrainbow)

dev.off()
