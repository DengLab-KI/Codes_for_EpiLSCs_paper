library(Seurat)
#pseudo_seu <- readRDS("/Users/qinluo/OneDrive - KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/pseudo_seu.rds")
plate1 = subset(pseudo_seu, subset=batch=="Lab")
#plate2 = readRDS("/Users/qinluo/OneDrive - KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/INTPSC_2.rds")
#all = merge(plate1, plate2)

chir = subset(plate1,subset=lineage=="CHIR")
chir$id = chir$lineage
#all = merge(chir, intpsc2)
#DefaultAssay(all)="RNA"

DefaultAssay(chir) ="RNA"
DefaultAssay(intpsc2) ="RNA"
chir = NormalizeData(chir)
chir = ScaleData(chir,features = rownames(chir))
intpsc2 = NormalizeData(intpsc2)
intpsc2 = ScaleData(intpsc2,features = rownames(intpsc2))

count1= as.data.frame(chir@assays$RNA@counts)
count2=as.data.frame(intpsc2@assays$RNA@counts)
genes = intersect(rownames(count1),rownames(count2))
count1 = count1[genes,]
count2 = count2[genes,]

meta = intpsc2@meta.data
meta = as.data.frame(meta)
#ggplot(meta)+geom_bar(aes(x=lineage,fill=id),position = "fill")

meta$batch2 = 1
meta1 = meta[which(meta$lineage=="INTPSC"),]
meta2 = meta[which(meta$lineage=="INTPSC2"),]
set.seed(2266)
a=sample(1:252,120)
b=sample(1:276,140)
meta1$batch2[a] = "P17"
meta1$batch2[-a] = "P10"
meta2$batch2[b] = "P17"
meta2$batch2[-b] = "P10"

metanew = rbind(meta1,meta2)
metanew$batch2 = as.factor(metanew$batch2)
intpsc2$id2 = metanew$batch2
chir$id2 = chir$id

bm = merge(chir[genes,],intpsc2[genes,])
g = c(grep("mt-",rownames(ref)),grep("Rsp",rownames(ref)),
      grep("Rpl",rownames(ref)),
      grep("Gm",rownames(ref)))
ref = ref[-g,]
colnames(ref)=1:100
bm_stalt_corrected =staltIntegration2(count1,
                                      count2,  
                                      nPoints = 100,
                                      refk=ref,
                                      ref=ref)
bm_stalt_corrected = CreateSeuratObject(bm_stalt_corrected)
bm_stalt_corrected$id = bm$id
bm_stalt_corrected$id2 = bm$id2
bm_stalt_corrected = NormalizeData(bm_stalt_corrected)
bm_stalt_corrected = ScaleData(bm_stalt_corrected)
bm_stalt_corrected = FindVariableFeatures(bm_stalt_corrected)
bm_stalt_corrected = RunPCA(bm_stalt_corrected)
bm_stalt_corrected = RunUMAP(bm_stalt_corrected,dims = 1:15)
p1=DimPlot(bm_stalt_corrected,
        group.by = "id2",
        cols = c("#c1946a","grey","black"),
        pt.size = .8)+
    coord_fixed(ratio = .6)+theme_void()

p2=DimPlot(bm_stalt_corrected,
        group.by = "id",
        cols = c("#82b74b","#6b5b95",
                 "#dac292","#c1946a"),
        pt.size = .8)+
    coord_fixed(ratio = .6)+theme_void()

p1+p2
