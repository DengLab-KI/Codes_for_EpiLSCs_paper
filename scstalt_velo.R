######
setwd("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Desktop/0909/0306/velocity")

library(SeuratDisk)
library(anndata)
library(Seurat)
all.ad <- read_h5ad("EpiLSC.h5ad")
dy_scvelo <- LoadH5Seurat("EpiLSC.h5Seurat")

emat_all = as.data.frame(emat_all)
nmat_all = as.data.frame(nmat_all)

gene =intersect(rownames(emat_all),rownames(dy_scvelo))
ad_new <- import("anndata", convert = FALSE)
ad_new <- ad_new$AnnData(
    X= as.data.frame(all.ad$X[,gene]),
    obs=all.ad$obs,
    var=all.ad$var[gene,],
    layers=list('spliced'=t(emat_all[gene,colnames(dy_scvelo)]), 
                'unspliced'=t(nmat_all[gene,colnames(dy_scvelo)])),
    obsm=list('X_umap'=dy_scvelo@reductions$pca@cell.embeddings[,1:2]) 
)

####
adata = ad_new
#library(reticulate)
#use_condaenv("r-velocity", required = TRUE)
#scv <- import("scvelo")
#saveRDS(adata,"adata.rds")
#adata <- scv$datasets$pancreas()
scv$pp$filter_genes(adata)
scv$pp$filter_genes(adata,min_cells=10,
                    min_counts_u = 1,
                    min_cells_u =5) ## filter

scv$pp$moments(adata) ## normalize and compute moments
scv$tl$recover_dynamics(adata,n_jobs=1) ## model
#saveRDS(adata,"adata.rds")
## plot (creates pop up window)
scv$tl$velocity(adata, mode='dynamical')
#saveRDS(adata,"adata.rds")
scv$tl$velocity_graph(adata)
scv$pl$velocity_embedding_stream(adata,
                                 basis='umap',
                                 color="royalblue",size=100)


scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color=c("id"),size=100,
                                 palette = c(c1="#82b74b",c2="#6b5b95",
                                             c3="#dac292"))

scv$pl$velocity_embedding_stream(adata, basis='umap',
                                 color="id",
                                 save = "EpiLC_new.svg",
                                 size=100,
                                 alpha=0.9,
                                 #dpi=300,
                                 palette = c(c1="#82b74b",c2="#6b5b95",
                                             c3="#dac292"))
topgenes <- adata$var["fit_likelihood"]
topgenes =scv$DataFrame(topgenes)

topgenes = data.frame(gene=rownames(topgenes)[
    order(topgenes$fit_likelihood,decreasing = T)],
    ll=topgenes[order(
        topgenes$fit_likelihood,decreasing = T),])
write.table(topgenes,"tore_dy.txt",sep="\t")

###########gnees

##############
all.ad2 <- read_h5ad("tore_new.h5ad")
tore_scvelo <-  LoadH5Seurat("tore_new.h5Seurat")
gene =intersect(rownames(emat_all),rownames(tore_scvelo))

ad_new2 <- import("anndata", convert = FALSE)
ad_new2 <- ad_new2$AnnData(
    X= as.data.frame(all.ad2$X[,gene]),
    obs=all.ad2$obs,
    var=all.ad2$var[gene,],
    layers=list('spliced'=t(emat_all[gene,colnames(tore_scvelo)]), 
                'unspliced'=t(nmat_all[gene,colnames(tore_scvelo)])),
    obsm=list('X_umap'=tore_scvelo@reductions$pca@cell.embeddings[,1:2]) 
)

####
adata2 = ad_new2
#saveRDS(adata,"adata.rds")
#adata <- scv$datasets$pancreas()
scv$pp$filter_genes(adata2) ## filter
scv$pp$filter_genes(adata,min_cells=10,
                    min_counts_u = 1,
                    min_cells_u =5) 
scv$pp$moments(adata2) ## normalize and compute moments
scv$tl$recover_dynamics(adata2,n_jobs=1) ## model
#saveRDS(adata,"adata.rds")
## plot (creates pop up window)
scv$tl$velocity(adata2, mode='dynamical')

topgenes <- adata2$var["fit_likelihood"]
topgenes =scv$DataFrame(topgenes)

topgenes = data.frame(gene=rownames(topgenes)[
    order(topgenes$fit_likelihood,decreasing = T)],
    ll=topgenes[order(
        topgenes$fit_likelihood,decreasing = T),])

write.table(topgenes,"tore_topgenes.txt",sep="\t")

#saveRDS(adata,"adata.rds")
scv$tl$velocity_graph(adata2)
scv$pl$velocity_embedding_stream(adata2, basis='umap',
                                 color="id",size=100,
                                 palette=clustercol[c(1,2,12,9)])
scv$pl$velocity_embedding_stream(adata2, basis='umap',
                                 color="id",size=100,
                                 palette=clustercol[c(8,1,2,11)],alpha = 1)
scv$pl$velocity_embedding_stream(adata2, basis='umap',
                                 color="id",
                                 save = "tore_new_velocity.svg",
                                 size=100,
                                 alpha=0.9,
                                 #dpi=300,
                                 palette=clustercol[c(8,1,2,11)])
clustercol=c("#eeac99", "#8ca3a3", "#77AADD", "orange", 
             "#99CCBB", "#117755","#777711",  
             "#ffef96",  "#c1946a", "#82b74b", "#034f84", 
             "#b8a9c9", "#F5DEB3", "#7e4a35","#778822") 
scv$pl$velocity(adata2,c("Esrrb","Nanog","Lef1"),
                ncols=2,color ="lineage2",
                palette=clustercol[c(1,2,12,9)],
                alpha=0.8)
scv$pl$velocity(adata,c("Esrrb","Nanog","Otx2"),
                ncols=2,color ="id",
                palette=c("#6b5b95","#feb236","#82b74b",
                          "#dac292","#92a8d1","#034f84"),
                alpha=0.8)

scv$pl$velocity(adata,names(topgenes_vals[1:20]),
                ncols=4,color ="id",
                palette=c("#6b5b95","#feb236","#82b74b",
                          "#dac292","#92a8d1","#034f84"),
                alpha=0.8,size=30)

scv$pl$velocity(adata2,c("Klf5","Tfcp2l1","Nanog","Esrrb","Dnmt3l",
                         "Usp9x","Fgf4","Lin28a","Dnmt3a","Hes1",
                         "Fgf8","T","Krt8"),
                ncols=2,color ="lineage2",
                palette=clustercol[c(1,2,11,8)],
                alpha=0.8,
                save = "tore_genes.png",
                size=30)

scv$pl$velocity(adata,c("Klf5","Tfcp2l1","Nanog","Esrrb","Dnmt3l",
                        "Usp9x","Fgf4","Lin28a","Dnmt3a","Hes1",
                        "Fgf8","T","Krt8"),
                ncols=2,color ="id",
                palette=c("#6b5b95","#feb236","#82b74b",
                          "#dac292","#92a8d1","#034f84"),
                save = "dy_genes.png",
                size=30,
                alpha=0.8)



scv$pl$velocity(adata,c("Tdh","Nr5a2","Nanos","Esrrb",
                        "Hmces",
                        "Lin28a","Dnmt3l","pim2",
                        "Krt8"),
                ncols=2,color ="id",
                palette=c("#6b5b95","#feb236","#82b74b",
                          "#dac292","#92a8d1","#034f84"),
                save = "dy_genes.png",
                size=30,
                alpha=0.8)
scv$pl$velocity(adata2,c("Tdh","Nr5a2","Nanos","Esrrb",
                         "Hmces",
                         "Lin28a","Dnmt3l","pim2",
                         "Krt8"),
                ncols=2,color ="lineage2",
                palette=clustercol[c(1,2,11,8)],
                alpha=0.8,
                save = "tore_genes.png",
                size=30)

## top dynamic genes,




all.ad3 <- read_h5ad("/Users/qinluo/Downloads/chir_scvelo.h5ad")
chir_velo <- LoadH5Seurat("/Users/qinluo/Downloads/chir_scvelo.h5Seurat")

ad_new3 <- import("anndata", convert = FALSE)
ad_new3 <- ad_new3$AnnData(
    X= as.data.frame(all.ad3$X[,gene]),
    obs=all.ad3$obs,
    var=all.ad3$var[gene,],
    layers=list('spliced'=t(emat_all[,colnames(chir_velo)]), 
                'unspliced'=t(nmat_all[,colnames(chir_velo)])),
    obsm=list('X_umap'=chir_velo@reductions$harmony@cell.embeddings[,1:2]) 
)
####
adata3 = ad_new3
#saveRDS(adata,"adata.rds")
#adata <- scv$datasets$pancreas()
scv$pp$filter_genes(adata3) ## filter
scv$pp$moments(adata3) ## normalize and compute moments
scv$tl$recover_dynamics(adata3,n_jobs=1) ## model
#saveRDS(adata,"adata.rds")
## plot (creates pop up window)
scv$tl$velocity(adata3, mode='dynamical')
#saveRDS(adata,"adata.rds")
scv$tl$velocity_graph(adata3)
scv$pl$velocity_embedding_stream(adata3, basis='umap',
                                 color = "orange",size=100)

scv$pl$velocity(adata3,c("Esrrb"),
                ncols=2,color ="id",
                palette=clustercol[c(1,2,11,8)],
                alpha = .8)




## top dynamic genes
topgenes <- adata$var["fit_likelihood"]
topgenes =scv$DataFrame(topgenes)
write.table(topgenes,"./Downloads/chir_topgenes.txt",sep="\t")

topgenes_vals <- topgenes[,1]
names(topgenes_vals) <- rownames(topgenes)
write.table(topgenes,"./Downloads/chir_topgenes.txt",sep="\t")

topgenes_vals <- sort(topgenes_vals, decreasing=TRUE)
head(topgenes_vals)

scv$pl$velocity(adata,names(topgenes_vals[1:5]),ncols=2)
scv$pl$velocity(adata2,names(topgenes_vals[1:5]),ncols=2)


scv$pl$velocity(adata,c("Esrrb","Nanog","Lef1"),
                ncols=2,color ="id",
                palette=clustercol[c(1,2,11,8)],
                alpha = .8)
scv$pl$velocity(adata2,c("Esrrb","Nanog","Lef1"),
                ncols=2,color ="lineage2",
                palette=clustercol[c(1,2,11,8)],
                alpha = .8)



