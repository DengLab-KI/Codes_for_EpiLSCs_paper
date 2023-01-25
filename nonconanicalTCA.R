library(Seurat)
library(ggplot2)
formative <- readRDS("/Users/qinluo/Library/CloudStorage/OneDrive-KI.SE/Mac/Documents/Project_PGC_scRNA_data_integration/data/formative.rds")
DefaultAssay(formative)="RNA"
formative = NormalizeData(formative)
formative = ScaleData(formative,features = rownames(formative))

tca= c("Cs","Acly","Aco2","Aco1","Idh1","Idh2",
       "Idh3g","Idh3a","Idh3b","Ogdh","Oghd1",
       "Dlst","Dld","Suclg2","Sucla2","Sdha",
       "Sdhb","Sdhc","Sdhd","Fh1","Mdh1","Mdh2",
       "Pdhb","Pdha1",
       "Dlat")

DoHeatmap(formative,features = tca,
          group.by = "lineage2",
          group.colors = c("grey" ,"#3b3a30",
                           "#f9d5e5","#c2d4dd"))+
    scale_fill_gradientn(colors =
                             colorRampPalette(c("#2874A6", "white","#C0392B"))(100))

DoHeatmap(formative,features = c("Slc25a1","Acly","Mdh1","Fasn","Acaca"),
          group.by = "lineage2",
          group.colors = c("grey" ,"#3b3a30",
                           "#f9d5e5","#c2d4dd"))+
    scale_fill_gradientn(colors =
                             colorRampPalette(c("#2874A6", "white","#C0392B"))(100))

