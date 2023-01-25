# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(scWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)
DimPlot(intpsc, group.by='id', label=FALSE) +
 theme_void() + ggtitle('Cell Type')

# load scWGCNA package
library(scWGCNA)
library(Seurat)
seurat_odc = intpsc
#seurat_odc$metacell_group <- seurat_odc$split

seurat_odc$metacell_group <-intpsc$id
seurat_odc = FindVariableFeatures(seurat_odc,nfeatures = 5000)
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
metacell_subset = LogNormalize(metacell_subset)
metacell_subset = ScaleData(metacell_subset,features = rownames(metacell_subset))
metacell_subset = FindVariableFeatures(metacell_subset,nfeatures = 3000)
metacell_subset = RunPCA(metacell_subset,npcs=20)

library(tidyverse)
library(WGCNA)
enableWGCNAThreads()
metacell_seurat =metacell_subset
# how many groups are there
nclusters <- length(unique(metacell_seurat$orig.ident))

# which genes are we using ?
genes.use <- rownames(metacell_seurat)

# cell meta-data table
targets <- metacell_subset@meta.data

# vector of cell conditions
group <- as.factor(metacell_seurat$orig.ident)

# format the expression matrix for WGCNA
datExpr <- as.data.frame(GetAssayData(metacell_seurat, assay='RNA', slot='data')[genes.use,])
datExpr <- as.data.frame(t(datExpr))

# only keep good genes:
datExpr <- datExpr[,goodGenes(datExpr)]

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

wnt_related_genes = read.table("/Users/qinluo/Downloads/GO_term_summary_20220210_093945.txt",
                               header = T)
wnt_Genes_id  = match(wnt_Genes,geneInfo$GeneSymbol)
wnt_Genes_module = geneInfo[wnt_Genes_id,]
wnt = data.frame(table(wnt_Genes_module$Initially.Assigned.Module.Color))
wnt$Var1 = factor(wnt$Var1,levels = wnt$Var1)
ggplot(wnt) +
  geom_bar(aes(x="",y=Freq,fill=Var1),stat="identity", width=1,position = "fill") +
  coord_polar("y", start=2)+
  theme_void()+ scale_fill_manual(values=c("#633517","#e06000","#ffab00",	"#004d33",
                                           "#00477e","black"))
DoHeatmap(metacell_seurat,
          features = wnt_Genes_module$GeneSymbol)
wnt = wnt[c(2,4,7,3,6,1),]
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



plot_df <- cbind(select(targets, c(odc_group, Condition)), PCvalues,cluster=metacell_seurat$orig.ident)
plot_df <- reshape2::melt(plot_df, id.vars = c('odc_group', 'Condition','cluster'))
#plot_df$odc_group <- factor(plot_df$odc_group, levels=c('progenitor', 'early', 'intermediate', 'mature'))

colors <- sub('ME', '', as.character(levels(plot_df$variable)))
p <- ggplot(plot_df, aes(x=odc_group, y=value, fill=cluster)) +
  geom_boxplot() +
  RotatedAxis() + ylab('Module Eigengene')  


# plot width and height:
w=2*nmodules; h=2*nclusters;

pdf('figures/ME_trajectory_Plot_celtype_condition.pdf',width=w,height=h,useDingbats=F)
p + facet_wrap(~variable, scales='free', ncol=3)+
  scale_fill_manual(values= c("#6b5b95","#feb236","#82b74b",
                 "#dac292","#92a8d1","#034f84"))+theme_classic()
dev.off()

