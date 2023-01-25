atp = read.table("/Users/qinluo/Desktop/atp.txt",header = T)

DoHeatmap(intpsc2,
          features = atp$A,
          group.by = "id")
DoHeatmap(intpsc2,
          features = atp$T,
          group.by = "id")
DoHeatmap(intpsc2,
          features = atp$P,
          group.by = "id")

DoHeatmap(intpsc2,
          features = atp$P,
          group.by = "id")
DoHeatmap(intpsc2,
            features = c("Dppa4","Spp1","Cyp2b23",
                        "Trh","Ntrk2","Sprr4","Fos",
                        "Kdr","T","Mesp1",
                        "Klhl6","Wnt5a"),
          group.by = "id")
intpsc2$id2 = intpsc2$i
length(intpsc2$id2[which(intpsc2$id2=="C3")]) 


dat = intpsc2@assays$RNA@counts
shuffle = sample(1:ncol(dat))
dat = dat[,shuffle]
new = CreateSeuratObject(dat)
new@meta.data = intpsc2@meta.data[shuffle,]
new = LogNormalize(new)
new = FindVariableFeatures(new)
new = ScaleData(new,features = rownames(new))
DoHeatmap(new,
          features = c("Pou5f1","Sox2","Sall4","Tdgf1",
                       "Utf1","Dppa3","Zic3",
                       "Nanog","Esrrb","Nr0b1","Tfcp2l1",
                       "Prdm14","Klf2",
                       "Otx2","Dnmt3b","Lef1","Fgf5","Sox4",
                       "Pou3f1","Sox3",
                       "Dppa4","Aire","Fam25c","Dppa2",
                       "Trh","Ntrk2","Klhl6","Aplnr",
                       "T","Sp5"),
          group.by = "id",slot = "scale.data")+
    scale_fill_gradientn(colors =
                             colorRampPalette(c("#2874A6", "white","#C0392B"))(100))
 

DoHeatmap(new,
          features = c(atp$A, atp$T,atp$P),
          group.by = "id",slot = "scale.data",
          group.colors = c("#82b74b","#6b5b95","#dac292"))+
    scale_fill_gradientn(colors =
          colorRampPalette(c("#2874A6", "#E7EFF5",
                             "white","#C0392B"))(10)
          )
          