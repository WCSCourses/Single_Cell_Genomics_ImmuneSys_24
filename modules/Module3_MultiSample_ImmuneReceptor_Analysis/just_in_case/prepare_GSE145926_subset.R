library(Seurat)
library(tidyverse)
library(SeuratDisk)



#This part is not included in the ggcollab tutorial, as requires more memory than we have there. Nevertheless, I am providing it for the future reference and reproducibility.
# Data is from GSE145926, we need the files in the table below (as they have available TCR data) 

meta_data=tibble(individual=c("C141",
                              "C142",
                              "C144",
                              "C141",
                              "C142",
                              "C144"),
                 type=c(rep("RNA",3),
                        rep("TCR",3)),
                 sample=c("GSM4339769",
                          "GSM4339770",
                          "GSM4339772",
                          "GSM4385990",
                          "GSM4385991",
                          "GSM4385993"))

#########################
######Reading the data in
#########################

#I will refrain from using "C144" as data is much worse than other samples, so taking just two
#reading in both samples, adding QC measures
seurat_list <- list()

for(i in c("C141", "C142")){
  print(i)
  gsm <- meta_data%>%dplyr::filter(individual==i, type=="RNA")%>%.$sample
  print(gsm)
  sc <- Read10X_h5(file.path("data/first_day/GSE145926/", paste0(gsm, "_",i, "_filtered_feature_bc_matrix.h5")))
  sc <- CreateSeuratObject(counts = sc, project = i,   min.cells = 5,
                           min.features = 500)
  
  
  sc[["perc.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  sc[["MALAT1"]] <- PercentageFeatureSet(sc, pattern = "MALAT1")
  sc[["perc.ribo"]] <- PercentageFeatureSet(sc, "^RP[SL]")
  sc[["mt_ribo"]] <- round(sc[["perc.mt"]]/sc[["perc.ribo"]],2)
  sc[["mt_ribo"]][sc[["mt_ribo"]]>100] <- 100
  sc$complexity <-round(sc$nCount_RNA/sc$nFeature_RNA, 2)
  seurat_list[[i]] <-sc
  rm(sc)
}

#########################
###### Merge, QC, clean - very simplistic
#########################

#Joining into one object
sc <-  merge(seurat_list[[1]], y = c(seurat_list[[2]]), add.cell.ids = names(seurat_list)[1:2])

#Investigating what would be good cutoffs for QC filtering
VlnPlot(.,features = c("perc.mt", "MALAT1","perc.ribo", "mt_ribo","nFeature_RNA","nCount_RNA", "complexity"), group.by = "orig.ident", pt.size=0)


nFeature_RNA_cutoff= 1000
perc.mt_cutoff = 10
sc <- subset(sc, subset = nFeature_RNA > nFeature_RNA_cutoff & perc.mt < perc.mt_cutoff )

#19563 features across 6741 samples within 1 assay

#########################
###### Normalisation, Dimensionality Reduction, Clustering
#########################


#variable quality of cells, so making sure cell viability (perc.mt treated as its proxy) is taken into account when normalising
vars.to.regress = "perc.mt"

sc <- SCTransform(sc, method = "glmGamPoi", verbose = FALSE,variable.features.n = 2000, vars.to.regress = vars.to.regress,
                  vst.flavor = "v2", assay="RNA") #by default 3K features


#We want to omit TCR/BCR genes in clustering, as they are not interesting and driving clsutering, especially if we would cluster just T cells. Keeping in the genes defining MAIT/NKT cells
TCRs_pattern="^TRA(J|V)|^TRB(J|V|D)"
BCRs_pattern="^IG(HV|HJ|KV|KJ|LV|LJ)"
TCR_BCR_genes <- setdiff(c(grep(TCRs_pattern,rownames(sc), value = T),
                           grep(BCRs_pattern, rownames(sc), value = T)),
                         c(mait=c("TRAV1-2","TRAJ33", "TRAJ20
                                    ", "TRAJ12"),iNKT=c("TRAV10","TRAJ18" )))

variable_features <- setdiff(sc@assays$SCT@var.features, TCR_BCR_genes)

sc <- RunPCA(sc, verbose = FALSE,  reduction.name="SCT_pca", features = variable_features  )

ElbowPlot(sc , ndims = 50, reduction = "SCT_pca")
sc <- FindNeighbors(sc, dims = 1:30,reduction = "SCT_pca")
sc <- FindClusters(sc, resolution =  0.3 , verbose = FALSE)
sc <- RunUMAP(sc, dims = 1:30,reduction = "SCT_pca", reduction.name="SCT_umap")


#########################
###### Understand what are the cells in the sample
#########################

#For understanding what cells are in these samples: plot clusters, find defining markers
DimPlot(sc, group.by = "orig.ident")
DimPlot(sc)
markers_of_interest <- c("CD4", "CD8A", "CD19",  "CD3E"  ,  "CD3D"  ,  "CD3G" ,"TRBC2","TRDC","TRGD" )
VlnPlot(sc, features = markers_of_interest)


sc <- PrepSCTFindMarkers(sc)
rna.markers <- FindAllMarkers(sc, assay = "SCT", only.pos = TRUE)

rna.markers %>%group_by(cluster)%>%
  arrange(p_val_adj, desc(avg_log2FC))%>%
  filter(pct.1-pct.2>0.2)%>%
  slice_head(n=5)%>%
  print(n=60)

#TCRs are just a small subset of the whole dataset (clusters 2,4,9)

#########################
###### Add the matching TCR data
#########################

#Read in the matching TCR data
tcrs <- list()

for(i in c("C141", "C142")){
  print(i)
  gsm <- meta_data%>%dplyr::filter(individual==i, type=="TCR")%>%.$sample
  print(gsm)
  tcrs[[i]] <- read_csv(file.path("data/first_day/GSE145926/", paste0(gsm, "_",i, "_filtered_contig_annotations.csv.gz")))
}

tcrs <- bind_rows(tcrs, .id = "orig.ident")

tcrs <- tcrs%>%mutate(
  full_barcode= paste(orig.ident, barcode,sep = "_")
)

#data is structured by chain rather than by barcode. It needs to be reshaped into barcode x tcr_features table. I am also deducting multiplet status from multiple TRB chains, more than 2 TRA chains.

tcrs_by_barcode <- tcrs%>%
  select(-c(is_cell, high_confidence,contig_id,raw_consensus_id))%>%
  pivot_wider(id_cols=full_barcode, names_from = chain, values_from=chain:raw_clonotype_id, values_fn = ~ paste(unique(.x), collapse=":"))%>%
  mutate(tcr_multiplet=grepl(pat=":",chain_TRB)|grepl(pat=":.*:",chain_TRA)|chain_Multi=="Multi")%>%
  mutate(cell= case_when(grepl(pat="TRAV1-2", v_gene_TRA) & grepl(pat="TRAJ12|TRAJ20|TRAJ33",j_gene_TRA) ~ "MAIT",
                         grepl(pat="TRAV10", v_gene_TRA) & grepl(pat="TRAJ18",j_gene_TRA) ~ "NKT",
                         !is.na(chain_TRG) ~"GD",
                         !is.na(chain_Multi) ~"multi",
                         3>1 ~ "abT"))%>%
  rename(raw_clonotype_id = raw_clonotype_id_TRA)%>%
  select(-c(starts_with("raw_clonotype_id_")))

# adding the whole tcrs_by_barcode to the meta.data slot
tcrs_by_barcode  <- data.frame(tcrs_by_barcode )
rownames(tcrs_by_barcode ) <- tcrs_by_barcode $full_barcode


sc <- AddMetaData(sc, tcrs_by_barcode %>%select(-c(full_barcode,
                                    chain_Multi, chain_TRG, 
                                    v_gene_Multi, v_gene_TRG ,
                                    starts_with("d_gene"),
                                    j_gene_Multi, j_gene_TRG,
                                    c_gene_Multi, c_gene_TRG,
                                    full_length_Multi, full_length_TRG,
                                    cdr3_Multi, cdr3_TRG,
                                    cdr3_nt_Multi ,cdr3_nt_TRG,
                                    productive_Multi, productive_TRG,
                                    starts_with("reads"),
                                    umis_Multi, umis_TRG)))

sc@meta.data <- sc@meta.data%>%
  mutate(raw_clonotype_id=ifelse(!is.na(raw_clonotype_id), paste(orig.ident, raw_clonotype_id, sep="_"), raw_clonotype_id))

#########################
###### Removing cells established to be multiplets based on TCR data
#########################

sc <- subset(sc, subset = tcr_multiplet, invert=TRUE)


#########################
###### Identyfying clusters with alpha beta T cels
#########################

DimPlot(sc, group.by = "cell")
sc@meta.data%>%group_by(seurat_clusters,cell)%>%
  ggplot(aes(x=seurat_clusters, fill=cell))+geom_bar()


#########################
###### Reducing the dataset - so it includes only T cells and necessary info
#########################

object.size(sc)%>%format(unit='auto')
saveRDS(sc, file=file.path("objects_files_participants/sc.RDS"))
sc <- readRDS("objects_files_participants/sc.RDS")
 
sc_tcells_only <- subset(sc, subset = seurat_clusters%in%c(2,4,9))

object.size(sc_tcells_only)%>%format(unit='auto')
saveRDS(sc_tcells_only, file=file.path("objects_files_participants/sc_tcells_only.RDS"))

sc_tcells_only_diet <-DietSeurat(
  sc_tcells_only,
  layers = c("counts","data"),
  features = c(sc@assays$SCT@var.features,"CD4"),
  assays = c("SCT","RNA"),
  dimreducs = c("SCT_umap"))

object.size(sc_tcells_only_diet)%>%format(unit='auto')
saveRDS(sc_tcells_only_diet, file=file.path("objects_files_participants/sc_tcells_only_diet.RDS"))

sc_tcells_only_diet_clean <- sc_tcells_only_diet
sc_tcells_only_diet_clean@meta.data <- sc_tcells_only_diet_clean@meta.data[,1:12]


#version without TCR info -  this one is for the participants



object.size(sc_tcells_only_diet)%>%format(unit='auto')
#[1] "53.6 Mb"

saveRDS(sc_tcells_only_diet_clean, file=file.path("objects_files_participants/sc_tcells_only_diet_clean.RDS"))



#########################
###### Exercise - Seurat object with T cell data
#########################

sc <- readRDS("objects_files_participants/sc_tcells_only_diet_clean.RDS")

#Have a look on this object
sc

DimPlot(sc)

head(sc@meta.data)

markers_of_interest <- c("CD4", "CD8A",  "CD3E"  ,  "CD3D"  ,  "CD3G" ,"TRBC2","TRDC" )
VlnPlot(sc, features = markers_of_interest)

FeaturePlot(sc, features = markers_of_interest)

#########################
###### Exercise - reading TCRs
#########################
tcrs <- list()

for(i in c("C141", "C142")){
  print(i)
  gsm <- meta_data%>%dplyr::filter(individual==i, type=="TCR")%>%.$sample
  print(gsm)
  tcrs[[i]] <- read_csv(file.path("data/first_day/GSE145926/", paste0(gsm, "_",i, "_filtered_contig_annotations.csv.gz")))
}

tcrs <- bind_rows(tcrs, .id = "orig.ident")

tcrs <- tcrs%>%mutate(
  full_barcode= paste(orig.ident, barcode,sep = "_")
)

#data is structured by chain rather than by barcode. It needs to be reshaped into barcode x tcr_features table. I am also deducting multiplet status from multiple TRB chains, more than 2 TRA chains.

tcrs_by_barcode <- tcrs%>%
  select(-c(is_cell, high_confidence,contig_id,raw_consensus_id))%>%
  pivot_wider(id_cols=full_barcode, names_from = chain, values_from=chain:raw_clonotype_id, values_fn = ~ paste(unique(.x), collapse=":"))%>%
  mutate(tcr_multiplet=grepl(pat=":",chain_TRB)|grepl(pat=":.*:",chain_TRA)|chain_Multi=="Multi")%>%
  mutate(cell= case_when(grepl(pat="TRAV1-2", v_gene_TRA) & grepl(pat="TRAJ12|TRAJ20|TRAJ33",j_gene_TRA) ~ "MAIT",
                         grepl(pat="TRAV10", v_gene_TRA) & grepl(pat="TRAJ18",j_gene_TRA) ~ "NKT",
                         !is.na(chain_TRG) ~"GD",
                         !is.na(chain_Multi) ~"multi",
                         3>1 ~ "abT"))%>%
  rename(raw_clonotype_id = raw_clonotype_id_TRA)%>%
  select(-c(starts_with("raw_clonotype_id_")))

# adding the whole tcrs_by_barcode to the meta.data slot
tcrs_by_barcode  <- data.frame(tcrs_by_barcode )
rownames(tcrs_by_barcode ) <- tcrs_by_barcode $full_barcode


sc <- AddMetaData(sc, tcrs_by_barcode %>%select(-c(full_barcode,
                                                   chain_Multi, chain_TRG, 
                                                   v_gene_Multi, v_gene_TRG ,
                                                   starts_with("d_gene"),
                                                   j_gene_Multi, j_gene_TRG,
                                                   c_gene_Multi, c_gene_TRG,
                                                   full_length_Multi, full_length_TRG,
                                                   cdr3_Multi, cdr3_TRG,
                                                   cdr3_nt_Multi ,cdr3_nt_TRG,
                                                   productive_Multi, productive_TRG,
                                                   starts_with("reads"),
                                                   umis_Multi, umis_TRG)))

sc@meta.data <- sc@meta.data%>%
  mutate(raw_clonotype_id=ifelse(!is.na(raw_clonotype_id), paste(orig.ident, raw_clonotype_id, sep="_"), raw_clonotype_id))



#########################
###### Clonotypes - are there any expansions?
#########################
#count occurences of repeated clonotypes, see the top ones
sc@meta.data$raw_clonotype_id%>%table()%>%sort(decr=T)%>%head()


sc@meta.data%>%select(raw_clonotype_id, orig.ident)%>%
  group_by(orig.ident,raw_clonotype_id )%>%
  summarise(N=n())

#see how many cases of expanded clonotypes 
sc@meta.data%>%select(raw_clonotype_id, orig.ident)%>%
  group_by(orig.ident,raw_clonotype_id )%>%
  summarise(N=n())%>%
  group_by(orig.ident,clonotype_size=N)%>%
  summarise(N=n())


###visualise

clonotypes_to_plot <- sc@meta.data%>%
  filter(!is.na(raw_clonotype_id))%>%
  group_by(orig.ident, raw_clonotype_id)%>%
  summarise(N=n())%>%
  arrange(desc(N))

clonotypes_to_plot%>%
  mutate(raw_clonotype_id=factor(raw_clonotype_id, levels=clonotypes_to_plot$raw_clonotype_id))%>%
  ggplot(aes(x=raw_clonotype_id, y=N))+geom_col()

### Diversity-richness-evenness
vegan::diversity(clonotypes_to_plot$N)

clonotypes_to_plot%>%
  group_by(orig.ident)%>%
  summarise(diversity=vegan::diversity(N),
            richness=length(N),
            evenness=diversity/log(length(N)))


###clusters with expanded clonotypes
sc_tcells_only_diet@meta.data%>%group_by(orig.ident, seurat_clusters, is_expanded=raw_clonotype_id%in%(clonotypes_to_plot%>%filter(N>2)%>%.$raw_clonotype_id))%>%
  summarise(N=n())%>%
  group_by(orig.ident, seurat_clusters)%>%
  mutate(fraq=N/sum(N))

#are they CD8 or CD4s
VlnPlot(sc_tcells_only_diet, features=c("CD4", "CD8A"))


#see expanded clonotypes on the DimPlot
sc_tcells_only_diet@meta.data <-sc_tcells_only_diet@meta.data%>%
  mutate(clono_to_identify = ifelse(raw_clonotype_id%in%(clonotypes_to_plot%>%filter(N>5)%>%.$raw_clonotype_id), raw_clonotype_id,"other"))


DimPlot(sc_tcells_only_diet, cells.highlight = WhichCells(object = sc_tcells_only_diet,expression =clono_to_identify == "other", invert=TRUE))+ ggtitle("All expanded clonotypes")
DimPlot(sc_tcells_only_diet, group.by="seurat_clusters", cells.highlight = WhichCells(object = sc_tcells_only_diet,expression =clono_to_identify == "C141_clonotype1")) + ggtitle("C141_clonotype1")

