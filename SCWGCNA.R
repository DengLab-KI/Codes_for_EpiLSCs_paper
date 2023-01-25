# load scWGCNA package
library(scWGCNA)
library(Seurat)
emb_new <- readRDS("/Users/qinluo/OneDrive - KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/emb_new.rds")
emb_subset = subset(emb_new,subset=lineage2=="fPSC"|lineage2=="FS"|lineage2=="FTW"|lineage2=="INTPSC")
# construct metacells by ODC group, condition
seurat_odc = emb_subset
#seurat_odc$metacell_group <- seurat_odc$split

seurat_odc$metacell_group <- paste0(
  as.character(seurat_odc$split), '_',
  as.character(seurat_odc$lineage2))
seurat_odc = FindVariableFeatures(seurat_odc,nfeatures = 21627)
genes.keep <- VariableFeatures(seurat_odc)
# loop through each group and construct metacells
seurat_list <- list()
for(group in unique(seurat_odc$metacell_group)){
  print(group)
  cur_seurat <- subset(seurat_odc, metacell_group == group)
  cur_seurat <- cur_seurat[genes.keep,]
  cur_metacell_seurat <- scWGCNA::construct_metacells(
    cur_seurat, name=group,
    k=20, #reduction='umap',
    assay='RNA', slot='count'
  )
  cur_metacell_seurat$Condition <- as.character(unique(cur_seurat$split))
  cur_metacell_seurat$odc_group <- as.character(unique(cur_seurat$lineage2))
  #cur_metacell_seurat$seurat_clusters <- as.character(unique(cur_seurat$seurat_clusters))
  
  seurat_list[[group]] <- cur_metacell_seurat
}

# merge all of the metacells objects
metacell_subset <- merge(seurat_list[[1]], seurat_list[2:length(seurat_list)])
#saveRDS(metacell_seurat, file='data/metacell_seurat.rds')
#metacell_subset <- ScaleData(metacell_subset, features = rownames(metacell_subset))
#metacell_subset <- RunPCA(metacell_subset, features=rownames(metacell_subset))
#library(harmony)
#metacell_subset <- RunHarmony(metacell_subset, group.by=c('Condition','odc_group'), dims=1:15)
#metacell_subset <- RunUMAP(metacell_subset, reduction='harmony', dims=1:15)
#DimPlot(metacell_subset, group.by="odc_group",        reduction='umap', label=TRUE,ncol=3) 
metacell_subset = SCTransform(metacell_subset,vars.to.regress = "Condition",return.only.var.genes = F)
metacell_subset = ScaleData(metacell_subset,features = rownames(metacell_subset))
metacell_subset = FindVariableFeatures(metacell_subset,nfeatures = 3000)
metacell_subset = RunPCA(metacell_subset,npcs=20)
ElbowPlot(metacell_subset)
metacell_subset <- RunHarmony(metacell_subset,dims.use=1:20, 
                   group.by.vars=c("Condition"),assay.use = "SCT")
metacell_subset <- RunUMAP(metacell_subset, 
                           reduction = "pca",
                dims=1:20,
                metric="cosine",
                local.connectivity=10,
                min.dist=.3,n.neighbors=10,
                spread=1)
metacell_subset <- FindNeighbors(metacell_subset  , 
                                 dims = 1:20,
                       reduction = "pca",
                       annoy.metric="cosine",
                       nn.method = "annoy")
metacell_subset   <- FindClusters(metacell_subset ,  
                                  resolution = .1,
                                  method="igraph",graph.name = "SCT_snn",
                       group.singletons=F,
                       algorithm = 2)
library(ggplot2)
p31=DimPlot(metacell_subset,label = T)+coord_fixed()+theme_void()
p32=DimPlot(metacell_subset, reduction = "umap",
            group.by = "odc_group",label = T,
            cols=Myrainbow,
            pt.size = 1)+ coord_fixed(ratio = 1)+
  theme_void()
p33=DimPlot(metacell_subset, reduction = "umap",
            group.by = "odc_group",label = T,
            cols=Myrainbow[-1],
            pt.size = 1,split.by = "Condition",ncol = 2)+ 
  coord_fixed(ratio = 1)+theme_void()
p31+p32+p33

library(tidyverse)
library(WGCNA)
enableWGCNAThreads()
metacell_seurat =metacell_subset
# how many groups are there
nclusters <- length(unique(metacell_seurat$odc_group))

# which genes are we using ?
genes.use <- rownames(metacell_seurat)

# cell meta-data table
targets <- metacell_subset@meta.data


# vector of cell conditions
group <- as.factor(metacell_seurat$Condition)

# format the expression matrix for WGCNA
datExpr <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[genes.use,])
datExpr <- as.data.frame(t(datExpr))

# only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]

#Next we perform a parameter search to identify a good soft power threshold.



# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));

# Call the network topology analysis function for each set in turn
powerTable = list(
  data = pickSoftThreshold(
    datExpr,
    powerVector=powers,
    verbose = 100,
    networkType="signed",
    corFnc="bicor"
  )[[2]]
);

# Plot the results:
pdf("1_Power.pdf", height=10, width=18)

colors = c("blue", "red","black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "mean connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (col in 1:length(plotCols)){
  ylim[1, col] = min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
  ylim[2, col] = max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE);
}

# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;

for (col in 1:length(plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
       xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
       main = colNames[col]);
  addGrid();
  
  if (col==1){
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2],
         labels=powers,cex=cex1,col=colors[1]);
  } else
    text(powerTable$data[,1], powerTable$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[1]);
  if (col==1){
    legend("bottomright", legend = 'Metacells', col = colors, pch = 20) ;
  } else
    legend("topright", legend = 'Metacells', col = colors, pch = 20) ;
}
dev.off()

#Based on the soft power threshold that you have selected, now we can build the co-expression network. Consult the WGCNA documentation if you need help selecting a soft power value. Furthermore, you should read the function description for blockwiseConsensusModules carefully to select the different parameters, however the ones that I have chosen here have generally given good results on a variety of datasets.



softPower=20

nSets = 1
setLabels = 'ODC'
shortLabels = setLabels

multiExpr <- list()
multiExpr[['ODC']] <- list(data=datExpr)

checkSets(multiExpr) # check data size

# construct network
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson",
                              power = softPower,
                              consensusQuantile = 0.3,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0.8,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 4,
                              pamStage=FALSE,
                              detectCutHeight = 0.995, minModuleSize = 50,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")


consMEs = net$multiMEs;
moduleLabels = net$colors;

# Convert the numeric labels to color labels
moduleColors = as.character(moduleLabels)
consTree = net$dendrograms[[1]];

# module eigengenes
MEs=moduleEigengenes(multiExpr[[1]]$data, colors = moduleColors, nPC=1)$eigengenes
MEs=orderMEs(MEs)
meInfo<-data.frame(rownames(datExpr), MEs)
colnames(meInfo)[1]= "SampleID"

# intramodular connectivity
KMEs<-signedKME(datExpr, MEs,outputColumnName = "kME",corFnc = "bicor")

# compile into a module metadata table
geneInfo=as.data.frame(cbind(colnames(datExpr),moduleColors, KMEs))

# how many modules did we get?
nmodules <- length(unique(moduleColors))

# merged gene symbol column
colnames(geneInfo)[1]= "GeneSymbol"
colnames(geneInfo)[2]= "Initially.Assigned.Module.Color"

# save info
write.csv(geneInfo,file=paste0('geneInfoSigned.csv'))

PCvalues=MEs

#Next we will visualize the WGCNA dendrogram:
  
  
  
  pdf("figures/SignedDendro.pdf",height=5, width=8)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = paste0("ODC lineage gene dendrogram and module colors"))
dev.off()

#ODC WGCNA Dendrogram

#Visualizations
#Here I will show some more visualizations you can use to show the results of the co-expression analysis:
  
#  Module Trajectories:
#  We plot the distribution of module eigengenes in each of the differnt ODC subgroups and split by the different disease conditions.



plot_df <- cbind(select(targets, c(odc_group, Condition)), PCvalues,cluster=metacell_seurat$seurat_clusters)
plot_df <- reshape2::melt(plot_df, id.vars = c('odc_group', 'Condition','cluster'))
#plot_df$odc_group <- factor(plot_df$odc_group, levels=c('progenitor', 'early', 'intermediate', 'mature'))

colors <- sub('ME', '', as.character(levels(plot_df$variable)))
p <- ggplot(plot_df, aes(x=odc_group, y=value, fill=cluster)) +
  geom_boxplot() +
  RotatedAxis() + ylab('Module Eigengene')  
  

# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf('figures/ME_trajectory_Plot_celtype_condition.pdf',width=w,height=h,useDingbats=F)
p + facet_wrap(~variable, scales='free', ncol=3)
dev.off()

