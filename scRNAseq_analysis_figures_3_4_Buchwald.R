library(Seurat)
library(ggplot2)
library(devtools)
library(scales)
library(RColorBrewer)
library(monocle3)
library(SeuratWrappers)
library(dplyr)
library(UCell)
library(Nebulosa)
library(readr)
library(dittoSeq)
library(ggpubr)
library(cowplot)



setwd("/Users/vdhere/Desktop/Buchwald/SeuratObjects")


#load in figure 3 seurat object
Fig3_v_final_Buchwald<-readRDS('Fig3_v_final_Buchwald.rds')

#figure 3a
DimPlot(Fig3_v_final_Buchwald, group.by = 'origin', label = T)

DimPlot(Fig3_v_final_Buchwald, group.by = 'celltypes3', label = T)


#figure 3c- 
dittoBarPlot(
  object = Fig3_v_final_Buchwald,
  var = "celltypes3", 
  group.by = "Batch", 
  theme=theme_half_open(),,
  main = "",xlab=NULL,  )+ RotatedAxis()+ theme(axis.text.x = element_text(size = 10,face = "bold"))+ theme(axis.text.y.left = element_text(size = 12,face = "bold")) 
table(DN2s_male_female$orig.ident)

Idents(Fig3_v_final_Buchwald)<-'celltypes3'
#retrieve proportion table for donut plots
cluster_membership_by_Study_ID<-as.data.frame(prop.table(table(Idents(Fig3_v_final_Buchwald), Fig3_v_final_Buchwald$Batch), margin = 2))
write.csv(cluster_membership_by_Study_ID,'cluster_membership_by_Study_ID.csv')


#FIGURE 3D

gyrdpu <- grDevices::colorRampPalette(c('#e5e7e9', RColorBrewer::brewer.pal(n=9, name="RdPu")))(100)

FeaturePlot(object = Fig3_v_final_Buchwald,ncol=5,
            features = c("Tcf7",'Il7r','Sell','Ccr7',
                         'Cd69','Junb','Fos','Jun',
                         'Tox','Ctla4','Pdcd1','Havcr2',
                         'Gzmb','Ifng','Ly6a','Mki67'
            ), cols=gyrdpu,pt.size = 2, min.cutoff = .1, order=T,
            raster=FALSE,label = FALSE, label.size = 5)&NoAxes()+ theme_void()&NoLegend()


#Figure 3E

Idents(Fig3_v_final_Buchwald)<-'celltypes3'
Fig1_woTcyc<-subset(Fig3_v_final_Buchwald,  idents = c("Tcyc"), invert=T)
Fig1_woTcyc <- NormalizeData(object = Fig1_woTcyc)

DotPlot(object = Fig1_woTcyc, dot.scale = 8,
        features = c(
          "Tcf7",'Il7r','Sell',
          'Fos','Jun','Cd69',
          'Isg15','Irf7','Ly6c2','Ly6a',
          'Tox','Lag3','Ctla4',
          'Pdcd1','Entpd1','Havcr2','Ifng','Gzmb'), 
        cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 8,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))

#figure 3F
library(monocle3)
cds <- as.cell_data_set(Fig3_v_final_Buchwald)
cds <- preprocess_cds(cds)
cds <- learn_graph(cds, use_partition = F)
cds <- order_cells(cds, reduction_method = 'UMAP' )

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           group_label_size = 0, 
           label_roots = F, 
           label_leaves = FALSE,
           label_branch_points = FALSE)

#figures 3G, 3H
markers <- list()

markers$stem<-c('Tcf7','Sell','Il7r','Fos','Junb','Cd69')
markers$exhaustion<-c('Lag3','Ctla4','Entpd1','Havcr2','Tigit','Pdcd1')

Fig3_v_final_Buchwald <- AddModuleScore_UCell(Fig3_v_final_Buchwald, 
                                              features = markers)
plot_density(
  Fig3_v_final_Buchwald,
  c('stem_UCell','exhaustion_UCell'),
  slot = NULL,
  joint = F,
  reduction = NULL,
  dims = c(1, 2),
  method = c("ks", "wkde"),
  adjust = 1,
  size = 1,
  shape = 16,
  combine = TRUE,
  pal = "viridis"
)& NoAxes()




#load in figure 4 seurat object
Fig4_v_final_Buchwald<-readRDS('Fig4_v_final_Buchwald.rds')

#figure 4A
DimPlot(Fig4_v_final_Buchwald, group.by = 'celltypes4', split.by = 'Treatment', label = T)

#figure 4A- retrieve proportion table for donut plots

Idents(Fig4_v_final_Buchwald)<-'celltypes4'
cluster_membership_by_Treatment<-as.data.frame(prop.table(table(Idents(Fig4_v_final_Buchwald), Fig4_v_final_Buchwald$Treatment), margin = 2))
write.csv(cluster_membership_by_Treatment,'cluster_membership_by_Treatment.csv')

#FIGURE 4B
Idents(Fig4_v_final_Buchwald)<-'Treatment'
DimPlot(Fig4_v_final_Buchwald,  split.by = 'Treatment', label = T)


#figure 4C- retrieve cell counts table for bar plots

Idents(Fig4_v_final_Buchwald)<-'celltypes4'
Fig4_v_final_Buchwald_wo_Tcyc_Tex<-subset(Fig4_v_final_Buchwald,  idents = c("Tcyc",'Tex'), invert=T)
cluster_membership_by_Treatment_stem_pops<-as.data.frame(table(Idents(Fig4_v_final_Buchwald_wo_Tcyc_Tex), Fig4_v_final_Buchwald_wo_Tcyc_Tex$Treatment), margin = 2)


PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGene(Fig4_v_final_Buchwald_wo_Tcyc_Tex ,genes =c('Tcf7')) 

#Find cells coexpressing a Gzmb, Ly6a, and Tcf7 for each treatment and for each CD8 T stem subpopulation

Tstem_control <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tprimed") & (Treatment == "Buchwald_control"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tstem_control, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tstem_control, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tstem_control, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tstem_control, vars = Gene2) > Gene2.cutoff 
                                & FetchData(Tstem_control, vars = Gene1) > Gene1.cutoff
                                & FetchData(Tstem_control, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tstem_control@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in control Tstems that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in control Tstems that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Tstems that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Tstems that co-express Gzmb + Ly6a + Tcf7

Teffstem_control <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Teffstem") & (Treatment == "Buchwald_control"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Teffstem_control, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Teffstem_control, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Teffstem_control, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Teffstem_control, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Teffstem_control, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Teffstem_control, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Teffstem_control@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in control Teffstems that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in control Teffstems that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Teffstems that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Teffstems that co-express Gzmb + Ly6a + Tcf7

Tpex_control <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tpex") & (Treatment == "Buchwald_control"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tpex_control, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tpex_control, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tpex_control, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tpex_control, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tpex_control, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tpex_control, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tpex_control@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in control Tpex that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in control Tpex that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Tpex that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in control Tpex that co-express Gzmb + Ly6a + Tcf7

Tstem_RT <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tprimed") & (Treatment == "Buchwald_RT"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tstem_RT, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tstem_RT, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tstem_RT, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tstem_RT, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tstem_RT, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tstem_RT, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tstem_RT@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT Tpex that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT Tpex that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT Tpex that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT Tpex that co-express Gzmb + Ly6a + Tcf7

Teffstem_RT <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Teffstem") & (Treatment == "Buchwald_RT"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Teffstem_RT, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Teffstem_RT, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Teffstem_RT, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Teffstem_RT, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Teffstem_RT, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Teffstem_RT, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Teffstem_RT@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT Teffstem that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT Teffstem that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT Teffstem that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of Teffstem in RT Teffstem that co-express Gzmb + Ly6a + Tcf7


Tpex_RT <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tpex") & (Treatment == "Buchwald_RT"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tpex_RT, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tpex_RT, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tpex_RT, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tpex_RT, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tpex_RT, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tpex_RT, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tpex_RT@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT Tpex that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT Tpex that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT Tpex that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of Tpex in RT Tpex that co-express Gzmb + Ly6a + Tcf7


######
Tstem_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tprimed") & (Treatment == "Buchwald_PD"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tstem_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tstem_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tstem_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tstem_PDL1, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tstem_PDL1, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tstem_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tstem_PDL1@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in PDL1 Tpex that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Tpex that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Tpex that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Tpex that co-express Gzmb + Ly6a + Tcf7

Teffstem_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Teffstem") & (Treatment == "Buchwald_PD"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Teffstem_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Teffstem_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Teffstem_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Teffstem_PDL1, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Teffstem_PDL1, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Teffstem_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Teffstem_RT@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in PDL1 Teffstem that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Teffstem that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Teffstem that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of Teffstem in PDL1 Teffstem that co-express Gzmb + Ly6a + Tcf7


Tpex_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tpex") & (Treatment == "Buchwald_PD"))
Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tpex_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tpex_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tpex_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tpex_PDL1, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tpex_PDL1, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tpex_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tpex_PDL1@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in PDL1 Tpex that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Tpex that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in PDL1 Tpex that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of Tpex in PDL1 Tpex that co-express Gzmb + Ly6a + Tcf7

Tstem_RT_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tprimed") & (Treatment == "Buchwald_RT_PD"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tstem_RT_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tstem_RT_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tstem_RT_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tstem_RT_PDL1, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tstem_RT_PDL1, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tstem_RT_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tstem_RT_PDL1@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT+PDL1 Tstems that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Tstems that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Tstems that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Tstems that co-express Gzmb + Ly6a + Tcf7


Teffstem_RT_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Teffstem") & (Treatment == "Buchwald_RT_PD"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Teffstem_RT_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Teffstem_RT_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Teffstem_RT_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Teffstem_RT_PDL1, vars = Gene2) > Gene2.cutoff 
                                & FetchData(Teffstem_RT_PDL1, vars = Gene1) > Gene1.cutoff
                          & FetchData(Teffstem_RT_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Teffstem_RT_PDL1@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT+PDL1 Teffstems that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that co-express Gzmb + Ly6a + Tcf7


Tpex_RT_PDL1 <- subset(Fig4_v_final_Buchwald, subset = (celltypes4 == "Tpex") & (Treatment == "Buchwald_RT_PD"))

Gene1 <- 'Gzmb' 
Gene2 <- 'Ly6a' 
Gene3 <- 'Tcf7' 
Gene1.cutoff <- 0 
Gene2.cutoff <- 0 
Gene3.cutoff <- 0 

# Time to party
Gene1.cells <- length(which(FetchData(Tpex_RT_PDL1, vars = Gene1) > Gene1.cutoff))
Gene2.cells <- length(which(FetchData(Tpex_RT_PDL1, vars = Gene2) > Gene2.cutoff))
Gene3.cells <- length(which(FetchData(Tpex_RT_PDL1, vars = Gene3) > Gene3.cutoff))

Gene1_Gene2_Gene3.cells <- length(which(FetchData(Tpex_RT_PDL1, vars = Gene2) > Gene2.cutoff 
                                        & FetchData(Tpex_RT_PDL1, vars = Gene1) > Gene1.cutoff
                                        & FetchData(Tpex_RT_PDL1, vars = Gene3) > Gene3.cutoff))

all.cells.incluster <- table(Tpex_RT_PDL1@active.ident)
Gene1.cells/all.cells.incluster*100 # Percentage of cells in RT+PDL1 Teffstems that express Gzmb
Gene2.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that express Ly6a
Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that express Tcf7
Gene1_Gene2_Gene3.cells/all.cells.incluster*100 #Percentage of cells in RT+PDL1 Teffstems that co-express Gzmb + Ly6a + Tcf7

#Figure 4G

PDL1_vs_control<- FindMarkers(object = Fig4_v_final_Buchwald,
                            group.by = 'Treatment', 
                            ident.1 = c('Buchwald_control'),
                            ident.2 = c('Buchwald_PD'),logfc.threshold = .5)

RT_vs_control<- FindMarkers(object = Fig4_v_final_Buchwald,
                                 group.by = 'Treatment', 
                                 ident.1 = c('Buchwald_control'),
                                 ident.2 = c('Buchwald_RT'),logfc.threshold = .5)

RT_PDL1_vs_control<- FindMarkers(object = Fig4_v_final_Buchwald,
                                       group.by = 'Treatment', 
                                       ident.1 = c('Buchwald_control'),
                                       ident.2 = c('Buchwald_RT_PD'),logfc.threshold = .5)


write.csv(RT_PDL1_vs_control,'RT_PDL1_vs_control_test.csv')
Total_DIIFS <- list(RT_vs_control=RT_vs_control,
                    PDL1_vs_control=PDL1_vs_control,
                    RT_PDL1_vs_control=RT_PDL1_vs_control)


deglist <- Total_DIIFS

deglist <- lapply(deglist[1:3], function(x){x <- x[(x$p_val_adj<0.05),]}) 
deglist.up <- lapply(deglist[1:3], function(x){x <- x[(x$avg_log2FC>0.1),]}) 
deglist.down <- lapply(deglist[1:3], function(x){x <- x[(x$avg_log2FC<(-0.1)),]}) 


df <- read_csv("lollipop_Figure_4G.csv")



df <- df %>%
  mutate(paired = rep(1:(n()/2),each=2),
         DEGS=factor(DEGS))
df %>% 
  ggplot(aes(x= DEGS_DIFF, y= reorder(Group,DEGS_DIFF), fill=DEGS)) +
  geom_col(position="dodge")+
  labs(y="Group")

df %>% 
  ggplot(aes(x= DEGS_DIFF, y= (Group), fill=DEGS)) +
  geom_col(position="dodge")+
  labs(y="Group")

df %>% 
  ggplot(aes(x= DEGS_DIFF, y= Group)) +
  geom_line(aes(group = paired))+
  geom_point(aes(color=DEGS), size=4) +
  theme(legend.position="top")


df %>% 
  ggplot(aes(x= DEGS_DIFF, y= reorder(Group,DEGS_DIFF))) +
  geom_line(aes(group = paired))+
  geom_point(aes(color=DEGS), size=4) +
  labs(y="Group")

df %>% 
  group_by(paired) %>%
  ggplot(aes(x= DEGS_DIFF, y= reorder(Group,DEGS_DIFF))) +
  geom_line(aes(group = paired),color="grey")+
  geom_point(aes(color=DEGS), size=7) +
  labs(y="Group")+
  theme_classic(18)+
  theme(legend.position="top") + 
  scale_color_brewer(palette="Set1", direction=-1)+
  scale_x_continuous(  n.breaks = 5)

df %>% 
  ggplot(aes(y= DEGS_DIFF, x= reorder(Group,DEGS_DIFF))) +
  geom_line(aes(group = paired),color="black")+
  geom_point(aes(color=DEGS), size=6) +
  labs(x="Group")+
  theme_classic(18)+
  theme(legend.position="top") + 
  scale_color_brewer(palette="Set1", direction=-1)+
  scale_y_continuous(  n.breaks = 5)



#figure 4H

DotPlot(object = Fig4_v_final_Buchwald, group.by  = 'Treatment', dot.scale = 8,
        features = c(
          'Tox', 'Dapl1', 'Ctla4', 'Dusp1',  'Btla',
          'Lamp1', 'Gzmb','Gzma' , 'Klrk1','Klrc1',
          'Cxcl10','Cxcr3', 'Ly6c2','Icam1',
          'Il18r1', 'Il18rap' , 'Ifngr1','Il7r',
          'Irf7','Isg15','Stat3'
        ), cols = c("RdYlBu") ) + RotatedAxis()+  theme(axis.text.x = element_text(size = 10,face = "bold"))+ theme(axis.text.y.left = element_text(size = 9,face = "bold"))


#Figure 4I

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Trajactory analysis
#'
#' This function is used to construct the possible evolutionary
#' trajectory of the incorporated cells based on their
#' @param object the object
#' 
devtools::load_all("/Users/vdhere/Downloads/monocle")



pkgs <- c(
  "Seurat", "SeuratWrappers", "ggplot2", "batchelor", "circlize",
  "dplyr", "optparse", "reshape2", "data.table", "magrittr",
  "patchwork", "scales", "GSVA", "RColorBrewer", "ggridges",
  "clusterProfiler", "survminer", "survminer", "monocle",
  "psych", "ggrepel", "pheatmap", "escape", "multcomp", "agricolae"
)

lapply(pkgs, function(x) require(package = x, character.only = TRUE, quietly = F, warn.conflicts = FALSE)) # nolint

monomatrix <- as(as.matrix(GetAssayData(Buchwald_RT_PD, slot = "counts")), "sparseMatrix")
feature_ann <- data.frame(gene_id = rownames(monomatrix), gene_short_name = rownames(monomatrix))
rownames(feature_ann) <- rownames(monomatrix)
monofd <- new("AnnotatedDataFrame", data = feature_ann)
sample_ann <- Buchwald_RT_PD@meta.data
monopd <- new("AnnotatedDataFrame", data = sample_ann)

monocds <- newCellDataSet(monomatrix,
                          phenoData = monopd,
                          featureData = monofd,
                          lowerDetectionLimit = 0.1,
                          expressionFamily = negbinomial.size()
)
head(pData(monocds))
head(fData(monocds))
monocds <- estimateSizeFactors(monocds)
monocds <- estimateDispersions(monocds)
monocds <- detectGenes(monocds, min_expr = 0.1)
print(head(fData(monocds)))
expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 50))
monocds <- monocds[expressed_genes, ]



RunMonocle <- function(object) {
  monomatrix <- as(
    as.matrix(GetAssayData(object, slot = "counts")), "sparseMatrix"
  )
  feature_ann <- data.frame(
    gene_id = rownames(monomatrix), gene_short_name = rownames(monomatrix)
  )
  rownames(feature_ann) <- rownames(monomatrix)
  monofd <- new("AnnotatedDataFrame", data = feature_ann)
  sample_ann <- object@meta.data
  monopd <- new("AnnotatedDataFrame", data = sample_ann)
  
  monocds <- newCellDataSet(monomatrix,
                            phenoData = monopd,
                            featureData = monofd,
                            lowerDetectionLimit = 0.1,
                            expressionFamily = negbinomial.size()
  )
  
  head(pData(monocds))
  head(fData(monocds))
  
  monocds <- estimateSizeFactors(monocds)
  monocds <- estimateDispersions(monocds)
  
  monocds <- detectGenes(monocds, min_expr = 0.1)
  print(head(fData(monocds)))
  expressed_genes <- row.names(subset(fData(monocds), num_cells_expressed >= 50)) # nolint
  monocds <- monocds[expressed_genes, ]
  
  disp_table <- dispersionTable(monocds)
  unsup_clustering_genes <- subset(
    disp_table, mean_expression >= 0.05 &
      dispersion_empirical >= 2 * dispersion_fit
  ) #
  monocds <- setOrderingFilter(monocds, unsup_clustering_genes$gene_id)
  
  monocds <- reduceDimension(
    monocds,
    max_components = 2,
    method = "DDRTree"
  )
  monocds <- orderCells(monocds)
  return(monocds)
}
?orderCells()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Hallmarker analysis
#'
#' This function is applied to calculate the hallmarker value for each cell
#'

RunEscape <- function(object,
                      lib = "H") {
  gshallmark <- getGeneSets(library = lib)
  esseurat <- enrichIt(
    obj = object, gene.sets = gshallmark,
    groups = 1000, cores = 30
  )
  write.table(esseurat, "Hall_Marker.txt", quote = F, sep = "\t") # nolint
  object <- Seurat::AddMetaData(object, esseurat)
  return(object)
}

## Trajactory analysis of CAF and NIF
Fig4_v_final_Buchwald <- subset(data, ident = 1:5)
Fig4_v_final_Buchwald@active.ident <- factor(
  paste0("c", Fig4_v_final_Buchwald@active.ident),
  levels = c("c1", "c2", "c4", "c3", "c5")
)
names(Fig4_v_final_Buchwald@active.ident) <- rownames(Fig4_v_final_Buchwald@meta.data)
y<-names(Fig4_v_final_Buchwald@active.ident)
x<-rownames(Fig4_v_final_Buchwald@meta.data)
Fig4_v_final_Buchwald@active.ident


Randomly_Selected_Cells <- append(
  Randomly_Selected_Cells,
  rownames(Fig4_v_final_Buchwald@meta.data[sample(1:nrow(Fig4_v_final_Buchwald@meta.data), 1500), ]) # nolint
)
Randomly_Selected_Cells <- c()
for (i in c('Tpex', 'Tprimed', 'Teffstem', 'Tex','Tcyc')) {
  Fig4_v_final_Buchwald <- subset(Fig4_v_final_Buchwald, ident = i)
  set.seed(1)
  Randomly_Selected_Cells <- append(
    Randomly_Selected_Cells,
    rownames(Fig4_v_final_Buchwald@meta.data[sample(1:nrow(Fig4_v_final_Buchwald@meta.data), 1500), ]) # nolint
  )
}


Fig4_v_final_Buchwald_Randomly <- subset(Fig4_v_final_Buchwald, cells = Randomly_Selected_Cells)
Fig4_v_final_Buchwald_Randomly$group <- ifelse(
  Fig4_v_final_Buchwald_Randomly@active.ident %in% c(1, 2, 4),
  "CAF", "NIF"
)

down<-subset(x = Fig4_v_final_Buchwald, downsample = 2000)

monocds <- RunMonocle(Buchwald_RT_PD)
monocds <- orderCells(monocds, root_state = 2)
monocds
pdf("Fig4_v_final_Buchwald_Trajactory_State.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
  monocds,
  cell_size = 2, color_by = "celltypes4",
  show_tree = FALSE, show_branch_points = FALSE
) +
  scale_color_manual(
    values = c('Tnaive'='#92D0EC',
               'Tcyc'='#6B6ACF',
               "Tprimed"='#E396DE',
               "Tpex" ='#892B73',
               "Tex"='#401593',
               'Teffstem'= '#025959')
  ) +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL,
  )
dev.off()

pdf("Fig4_v_final_Buchwald_Trajactory_Label.pdf", width = 4, height = 4)
plot_cell_trajectory(monocds,
                     cell_size = 1,
                     color_by = "label", show_tree = FALSE, show_branch_points = FALSE
) +
  scale_color_manual(values = colo[c(1, 4, 3, 9)])
dev.off()

## The percentage of clusters in each state
mat <- melt(table(monocds@phenoData@data[, c("State", "label")]))

p <- ggplot(mat, aes(x = State, y = value, fill = label)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(y = "Proportion (%)") +
  scale_fill_manual(values = pal_npg()(4)) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    legend.key.size = unit(0.6, "cm"),
    axis.text.x = element_text(
      color = "black", size = 10
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = c(0, 0.01), labels = percent)
ggsave("Fig4_v_final_Buchwald_Trajactory_Proportion.pdf", p, height = 2.5, width = 5)

beam_res <- BEAM(monocds, branch_point = 1, cores = 40)
beam_res <- beam_res[order(beam_res$qval), ]
beam_res <- beam_res[, c("gene_short_name", "pval", "qval")]

pdf("Fig4_v_final_Buchwald_Trajactory_Heatmap.pdf", height = 8, width = 6)
plot_genes_branched_heatmap(monocds[row.names(subset(
  beam_res,
  qval < 1e-4
)), ],
branch_point = 1,
num_clusters = 4,
cores = 30,
use_gene_short_name = T,
show_rownames = F
)
dev.off()

p <- plot_genes_branched_heatmap(monocds[row.names(subset(
  beam_res,
  qval < 1e-4
)), ],
branch_point = 1,
num_clusters = 4,
cores = 30,
use_gene_short_name = T,
show_rownames = T,
return_heatmap = T
)
write.table(
  p$annotation_row, "Fig4_v_final_Buchwald_Trajactory_Heatmap_Label.txt",
  quote = F, sep = "\t"
)

hallmark <- getGeneSets(library = "H")
data_set <- list(enriched = hallmark[[14]]@geneIds)
EMT_Scores <- gsva(
  as.matrix(Fig4_v_final_Buchwald_Randomly@assays$RNA@data), data_set,
  min.sz = 5, kcdf = "Poisson", method = "ssgsea",
  mx.diff = TRUE, verbose = FALSE, parallel.sz = 20
)

monocds@phenoData@data$EMT <- scale(as.vector(t(EMT_Scores)))
monocds@phenoData@data$EMT <- ifelse(
  monocds@phenoData@data$EMT > 1.5, 1.5, monocds@phenoData@data$EMT
)

pdf("Fig4_v_final_Buchwald_Trajactory_EMT.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
  monocds,
  cell_size = 0.5, color_by = "EMT",
  show_tree = FALSE, show_branch_points = FALSE
) +
  scale_color_gradientn(colors = rev(col1)) +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL
  )
dev.off()

EMT_Value <- data.frame(
  State = monocds@phenoData@data$State,
  EMT = as.vector(t(EMT_Scores)),
  subtype = monocds@phenoData@data$tissue
)

p <- ggplot(EMT_Value, aes(x = State, y = EMT)) +
  geom_jitter(aes(fill = State, color = State),
              width = 0.2, shape = 21, size = 0.5
  ) +
  scale_color_manual(
    values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
  ) +
  scale_fill_manual(
    values = c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
  ) +
  geom_boxplot(
    size = 0.6, fill = "white",
    outlier.fill = NA, outlier.color = NA, outlier.size = 0
  ) +
  theme_classic() +
  labs(title = "EMT scores") +
  theme(
    axis.title.y = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.ticks.length = unit(0.4, "lines"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(hjust = 0.5)
  )
ggsave("Fig4_v_final_Buchwald_Trajactory_EMT_BoxPlot.pdf", p, width = 2.8, height = 3.3)

CREB3L1_Matrix <- data.frame(
  t(monocds@reducedDimS), Fig4_v_final_Buchwald_Randomly@assays$RNA@data["CREB3L1", ]
)
colnames(CREB3L1_Matrix) <- c("Component1", "Component2", "CREB3L1")
CREB3L1_Matrix$group <- monocds@phenoData@data$State
CREB3L1_Matrix <- CREB3L1_Matrix[order(mat$CREB3L1), ]
pdf("creb3l1.pdf", width = 5, height = 3)
ggplot(CREB3L1_Matrix, aes(x = Component1, y = Component2, color = CREB3L1)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = rev(color_for_use)) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL,
  )
dev.off()

CREB3L1_Matrix$group <- paste0("State", CREB3L1_Matrix$group)
pdf("Fig4_v_final_Buchwald_Trajactory_CREB3L1_Box.pdf", width = 2, height = 3)
ggplot(CREB3L1_Matrix, aes(x = group, y = CREB3L1, fill = group)) +
  geom_boxplot(outlier.color = "white", outlier.size = 0) +
  theme_classic() +
  scale_fill_manual(
    values =
    c(pal_aaas()(9)[1], "#fe9929", pal_aaas()(9)[3])
  ) +
  labs(x = "", y = "Expression levels of CREB3L1") +
  theme(legend.position = "none")
dev.off()

CREB3L1_Matrix$tissue <- Fig4_v_final_Buchwald_Randomly$tissue

p <- ggplot(CREB3L1_Matrix, aes(x = tissue, y = CREB3L1, fill = group)) +
  scale_fill_manual(values = c(pal_aaas()(10), rev(pal_aaas()(5)))) +
  geom_boxplot(
    size = 0.6,
    outlier.fill = NA, outlier.color = NA, outlier.size = 0
  ) +
  theme_bw() +
  labs(title = "Expression levels of CREB3L1") +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      size = 11, angle = 45, hjust = 1, color = "black"
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
ggsave("e_state_creb3l1.pdf", p, width = 5.5, height = 3.3)

CREB3L1_Target_Genes <- read.table(
  "CREB3L1_Target_Genes.txt",
  header = T, sep = "\t"
)
CREB3L1_Target_Genes_Enrichment <- gsva(
  as.matrix(Fig4_v_final_Buchwald_Randomly@assays$RNA@data),
  list(CREB3L1_Target_Genes$gene),
  min.sz = 5, kcdf = "Poisson", method = "ssgsea",
  mx.diff = TRUE, verbose = FALSE, parallel.sz = 20
)
monocds@phenoData@data$target <- scale(
  as.vector(t(CREB3L1_Target_Genes_Enrichment))
)
monocds@phenoData@data$target <- ifelse(
  monocds@phenoData@data$target > 1.5, 1.5, monocds@phenoData@data$target
)

pdf("CREB3L1_Target_Genes_Enrichment.pdf", width = 3.5, height = 3.5)
plot_cell_trajectory(
  monocds,
  cell_size = 0.5, color_by = "target",
  show_tree = FALSE, show_branch_points = FALSE
) +
  scale_color_gradientn(colors = rev(col1)) +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL
  )
dev.off()


