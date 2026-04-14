#merFISH data process
human_to_mouse <- function(genes){
  # 全部转小写
  g <- tolower(genes)
  # 首字母大写
  g <- paste0(toupper(substr(g, 1, 1)), substr(g, 2, nchar(g)))
  return(g)
}
write.csv(data.frame(gene=human_to_mouse(unique(deg$gene))),file="seurat_cluster_deg.csv")

library(Seurat)
library(biomaRt)

library(homologene)
human_genes <- rownames(dat)
map_df <- homologene(human_genes, inTax = 9606, outTax = 10090)

# 转换成人-鼠同源基因
map_df <- homologene(human_genes, inTax = 9606, outTax = 10090)

head(map_df)
map_df <- na.omit(map_df)
colnames(map_df) <- c("human", "mouse")

convertHumanToMouse <- function(seurat_obj, assay = "RNA", 
                                inTax = 9606, outTax = 10090, 
                                method = c("first", "sum")) {
  # 加载依赖
  if (!requireNamespace("homologene", quietly = TRUE)) {
    stop("请先安装 homologene: remotes::install_github('oganm/homologene')")
  }
  method <- match.arg(method)
  
  message("提取基因名...")
  human_genes <- rownames(seurat_obj[[assay]]@counts)
  
  message("查询 homologene 映射...")
  map_df <- homologene::homologene(human_genes, inTax = inTax, outTax = outTax)
  map_df <- na.omit(map_df)
  colnames(map_df) <- c("human", "mouse")
  
  if (nrow(map_df) == 0) stop("没有找到任何基因映射！")
  
  # 处理多对一情况
  if (method == "first") {
    map_df <- map_df[!duplicated(map_df$human), ]
    gene_map <- setNames(map_df$mouse, map_df$human)
  } else if (method == "sum") {
    # 多个人类基因对应同一鼠基因 → 合并表达量
    library(Matrix)
    message("合并多对一映射的基因表达矩阵...")
    
    gene_map <- setNames(map_df$mouse, map_df$human)
    counts <- seurat_obj[[assay]]@counts
    common_genes <- intersect(rownames(counts), names(gene_map))
    
    new_counts <- Matrix::aggregate.Matrix(counts[common_genes, ], 
                                           groupings = gene_map[common_genes], 
                                           fun = "sum")
    new_data <- Matrix::aggregate.Matrix(seurat_obj[[assay]]@data[common_genes, ], 
                                         groupings = gene_map[common_genes], 
                                         fun = "sum")
    
    new_scale <- seurat_obj[[assay]]@scale.data
    if (nrow(new_scale) > 0) {
      common_scale <- intersect(rownames(new_scale), names(gene_map))
      new_scale <- Matrix::aggregate.Matrix(new_scale[common_scale, ], 
                                            groupings = gene_map[common_scale], 
                                            fun = "mean")
    }
    
    seurat_obj[[assay]]@counts <- new_counts
    seurat_obj[[assay]]@data <- new_data
    seurat_obj[[assay]]@scale.data <- new_scale
    message("转换完成 ✅")
    return(seurat_obj)
  }
  
  # === method = "first" 情况 ===
  message("应用一对一映射...")
  gene_map <- setNames(map_df$mouse, map_df$human)
  
  counts <- seurat_obj[[assay]]@counts
  data <- seurat_obj[[assay]]@data
  scale <- seurat_obj[[assay]]@scale.data
  
  common_genes <- intersect(rownames(counts), names(gene_map))
  
  new_counts <- counts[common_genes, ]
  rownames(new_counts) <- gene_map[rownames(new_counts)]
  
  new_data <- data[common_genes, ]
  rownames(new_data) <- gene_map[rownames(new_data)]
  
  if (nrow(scale) > 0) {
    common_scale <- intersect(rownames(scale), names(gene_map))
    new_scale <- scale[common_scale, ]
    rownames(new_scale) <- gene_map[rownames(new_scale)]
  } else {
    new_scale <- scale
  }
  
  seurat_obj[[assay]]@counts <- new_counts
  seurat_obj[[assay]]@data <- new_data
  seurat_obj[[assay]]@scale.data <- new_scale
  
  message("转换完成 ✅")
  return(seurat_obj)
}
#usage
# 默认：人类 (9606) → 小鼠 (10090)，取第一个映射
obj <- convertHumanToMouse(dat, assay = "RNA", method = "first")
# 如果想合并多个对应的人类基因 → 一个小鼠基因
mouse_obj <- convertHumanToMouse(seurat_obj, assay = "RNA", method = "sum")
dir.create("ref_h5ad")
setwd("ref_h5ad")
adata<-srt_to_adata(obj)
adata$write_h5ad("subtype.h5ad")
for( s in unique(merFISH_int_slt$orig.ident)){
  tmp<-subset(merFISH_int_slt,orig.ident==s)
  adata<-srt_to_adata(tmp)
  adata$write_h5ad(paste0(s,".h5ad"))
}
obj<-RunDEtest(obj,group_by="SubType",only.pos = T,fc.threshold = 1)
deg<-obj@tools$DEtest_SubType$AllMarkers_wilcox
deg<-deg[deg$p_val_adj<0.05,]
write.csv(data.frame(gene=unique(deg$gene)),file="seurat_cluster_deg.csv")

# run merFISH basic analysis

setwd("/nfs4/chaozhang/proj/Neuron/Aritra/merFISH/script")

# init
source("merFISH_pipeline.R")
source("functions.R")

update_geom_defaults("point", list(shape = 16))


### run pipeline per sample
# 
merFISH_mouse0.1 = merFISH_pipeline("../data/20210911_PFCL1", "mouse0.1") #notfound
merFISH_mouse0.2 = merFISH_pipeline("../data/20211004_PFCL2", "mouse0.2") #no

#
merFISH_mouse1.1 = merFISH_pipeline("../data/20211208_PFCL5", "mouse1.1") #not
merFISH_mouse1.2 = merFISH_pipeline("../data/20211208_PFCL6", "mouse1.2") #not 

#
merFISH_mouse2.1 = merFISH_pipeline("../data/20211024_PFCL1_Results_Pain-d1", "mouse2.1") #no
merFISH_mouse2.2 = merFISH_pipeline("../data/20211101_PFCL1_Results_Pain-d2", "mouse2.2") #no

#
merFISH_mouse3.1 = merFISH_pipeline("../data/20220103_PFCL7", "mouse3.1") #no
merFISH_mouse3.2 = merFISH_pipeline("../data/20220106_PFCL1_Results", "mouse3.2") 
merFISH_mouse3.3 = merFISH_pipeline("../data/20220107_PFCL1_Results", "mouse3.3") 
merFISH_mouse3.4 = merFISH_pipeline("../data/20220128_PFCL1_Results", "mouse3.4") 

#
merFISH_pain1 = merFISH_pipeline("../data/20220218_PFCL1_Results_Pain", "Pain1") 
merFISH_pain2 = merFISH_pipeline("../data/20220210_PFCL1_Results_Pain", "Pain2")

#
merFISH_pain3 = merFISH_pipeline("../data/20220227_PFCL1_Results_Cingulate and NAc", "Pain3") #no

#
merFISH_pain4 = merFISH_pipeline("../data/20220305_RESULTS_Pain and Control_Cingulate and NAc", "Pain4") #no

merFISH_pain5 = merFISH_pipeline("../data/20220416_PFCL1_Results", "Pain5") 

merFISH_pain6 = merFISH_pipeline("../data/20220429_PFCL1_Results", "Pain6") 

merFISH_pain7 = merFISH_pipeline("../data/20220430_PFCL1_Results", "Pain7") 

merFISH_pain8 = merFISH_pipeline("../data/20220513_PFCL1_Results", "Pain8")

merFISH_pain9 = merFISH_pipeline("../data/20220514_PFCL1_Results", "Pain9")

#"20211024_PFCL1_Results_Pain" "20211101_PFCL1_Results_Pain" "20211107_PFCL1_Results"      "20211121_PFCL1_Results"     
#"20211125_PFCL1_Results"   "20220214_PCL1_Results" "20221012_PFCL1_Results"      "20221022_PFCL1_Results" 
merFISH_mouse0.1= merFISH_pipeline("../data/20211024_PFCL1_Results_Pain", "20211024_PFCL1") 
merFISH_mouse0.2= merFISH_pipeline("../data/20211101_PFCL1_Results_Pain", "20211101_PFCL1")
merFISH_mouse1.1= merFISH_pipeline("../data/20211107_PFCL1_Results", "20211107_PFCL1")
merFISH_mouse1.2= merFISH_pipeline("../data/20211121_PFCL1_Results", "20211121_PFCL1")
merFISH_mouse2.1= merFISH_pipeline("../data/20211125_PFCL1_Results", "20211125_PFCL1")
merFISH_mouse2.2= merFISH_pipeline("../data/20220214_PCL1_Results", "20220214_PCL1")
merFISH_mouse3.1= merFISH_pipeline("../data/20221012_PFCL1_Results", "20221012_PFCL1")
merFISH_pain4= merFISH_pipeline("../data/20221022_PFCL1_Results", "20221022_PFCL1")
#以上没有，以下是有的-------------------------------------------------------------------------------------------------------------
# "20220106_PFCL1_Results"  ="mouse3.2"    "20220107_PFCL1_Results"="mouse3.3"     "20220128_PFCL1_Results"="mouse3.4"   
# "20220210_PFCL1_Results_Pain"="Pain2"       "20220218_PFCL1_Results_Pain"="Pain1" "20220416_PFCL1_Results" ="Pain5"    
# "20220429_PFCL1_Results" ="Pain6"     "20220430_PFCL1_Results"="Pain7"      "20220513_PFCL1_Results"="Pain8"      "20220514_PFCL1_Results" ="Pain9"   
merFISH_mouse3.2 = merFISH_pipeline("../data/20220106_PFCL1_Results", "20220106_PFCL1") 
merFISH_mouse3.3 = merFISH_pipeline("../data/20220107_PFCL1_Results", "20220107_PFCL1") 
merFISH_mouse3.4 = merFISH_pipeline("../data/20220128_PFCL1_Results", "20220128_PFCL1") 
merFISH_pain1 = merFISH_pipeline("../data/20220218_PFCL1_Results_Pain", "20220218_PFCL1") 
merFISH_pain2 = merFISH_pipeline("../data/20220210_PFCL1_Results_Pain", "20220210_PFCL1")
merFISH_pain5 = merFISH_pipeline("../data/20220416_PFCL1_Results", "20220416_PFCL1") 
merFISH_pain6 = merFISH_pipeline("../data/20220429_PFCL1_Results", "20220429_PFCL1") 
merFISH_pain7 = merFISH_pipeline("../data/20220430_PFCL1_Results", "20220430_PFCL1") 
merFISH_pain8 = merFISH_pipeline("../data/20220513_PFCL1_Results", "20220513_PFCL1")
merFISH_pain9 = merFISH_pipeline("../data/20220514_PFCL1_Results", "20220514_PFCL1")

### integrated
merFISH_integrated = merge(merFISH_mouse0.1, c(merFISH_mouse0.2,
                                               merFISH_mouse1.1, merFISH_mouse1.2,
                                               merFISH_mouse2.1, merFISH_mouse2.2,
                                               merFISH_mouse3.1, merFISH_mouse3.2, merFISH_mouse3.3, merFISH_mouse3.4,
                                               merFISH_pain1, merFISH_pain2,  merFISH_pain4, merFISH_pain5,
                                               merFISH_pain6, merFISH_pain7, merFISH_pain8, merFISH_pain9)) #merFISH_pain3,

saveRDS(merFISH_integrated, file = "./RDS/merFISH_integrated.seurat0627.rds")
rm(merFISH_mouse0.2,
     merFISH_mouse1.1, merFISH_mouse1.2,
     merFISH_mouse2.1, merFISH_mouse2.2,
     merFISH_mouse3.1, merFISH_mouse3.2, merFISH_mouse3.3, merFISH_mouse3.4,
     merFISH_pain1, merFISH_pain2, merFISH_pain3, merFISH_pain4, merFISH_pain5,
     merFISH_pain6, merFISH_pain7, merFISH_pain8, merFISH_pain9)
gc()

merFISH_integrated = readRDS(file = "./RDS/merFISH_integrated.seurat0627.rds")

# remove low quality 2 samples Pain6 Pain7 
merFISH_integrated = subset(merFISH_integrated, subset = orig.ident %in% c("Pain6", "Pain7"), invert = T)


#10x spatial data process
#https://github.com/satijalab/seurat/discussions/8926
#1 cp tissue_hires_image.png tissue_lowres_image.png
#2 modify scalefactors_json.json file, let tissue_lowres_scalef.value = tissue_hires_scalef.value 
matrix_dir = 'filtered_feature_bc_matrix/' 
counts = Seurat::Read10X(data.dir = matrix_dir)  
data = Seurat::CreateSeuratObject(
  counts = counts , 
  project = 'test', 
  assay = 'Spatial')
data$slice = 1 
data$region = 'test' 
imgpath = "spatial"
img = Seurat::Read10X_Image(image.dir = imgpath)  
Seurat::DefaultAssay(object = img) <- 'Spatial'  
img = img[colnames(x = data)]  
data[['image']] = img  
SpatialFeaturePlot(data, features = "nCount_Spatial")
saveRDS(data,file="sp.rds")


#load data
# 读取表达矩阵
data.dir="GSM4800808/"
setwd(data.dir)
filter_matrix<-Read10X("filtered_feature_bc_matrix")
library(DropletUtils)
write10xCounts("./filtered_feature_bc_matrix.h5", filter_matrix, type = "HDF5",
               genome = "grch38", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))
seurat_obj <- Load10X_Spatial(
  data.dir = data.dir,
  assay = "Spatial"
)

data <- data |>
  NormalizeData() |>
  ScaleData() |>
  FindVariableFeatures() |>
  RunPCA() |>
  RunUMAP(dims=1:20) |>
  FindNeighbors(reduction = "pca", dims = 1:10) |>
  FindClusters(resolution = c(0.2,0.5,1.2))
p<-DimPlot(data,group.by=c("Spatial_snn_res.0.2","Spatial_snn_res.0.5","Spatial_snn_res.1.2"))
p1<-SpatialDimPlot(data,group.by=c("Spatial_snn_res.0.2","Spatial_snn_res.0.5","Spatial_snn_res.1.2"))
pdf("sp_cls.pdf",width=16,height=5)
print(p)
print(p1)
dev.off()
saveRDS(data,file="sp_cls.rds")


library(magick,lib.loc=c("/home/yangmy/R/spatial", "/home/yangmy/R/cellchat","/home/yangmy/R/x86_64-pc-linux-gnu-library/4.3", 
                         "/opt/R/4.3.2/lib/R/library" ))
library(shinyjs,lib.loc=c("/home/yangmy/R/spatial", "/home/yangmy/R/cellchat","/home/yangmy/R/x86_64-pc-linux-gnu-library/4.3", 
                          "/opt/R/4.3.2/lib/R/library" ))
library(semla,lib.loc="/home/yangmy/R/spatial/")

infoTable <- tibble(samples=c("/data2/Project/anding/test_20250910/GSM4800808/filtered_feature_bc_matrix.h5"), 
                    imgs=c("/data2/Project/anding/test_20250910/GSM4800808/spatial/tissue_hires_image.png"), 
                    spotfiles=c("/data2/Project/anding/test_20250910/GSM4800808/spatial/tissue_positions_list.csv"), 
                    json=c("/data2/Project/anding/test_20250910/GSM4800808/spatial/scalefactors_json.json"), # Add required columns
                    sample_id = c("GSM480080")) # Add additional column
se<-ReadVisiumData(infoTable[1,])
saveRDS(se,file="sp_selma.rds")
spatial_data<-GetStaffli(se)
se<-LoadImages(se)
ImagePlot(se)
glist<-c("Slc17a7","Man1a","Syn1","Eno2",
         "Gad1","Gad2",
         "Ly86","Cx3cr1","Ctss",
         "Atp1a2", "Atp13a4",    "F3",
         "Plp1","Mbp","Mobp","Cldn11",
         "Cspg4","Pdgfra")
p<-MapFeatures(se,features=glist,image_use="raw",override_plot_dims=TRUE) & ThemeLegendRight()

se_merged <- MergeSTData(NormalizeData(sp[[1]]), NormalizeData(sp[[2]]))
MapFeatures(se_merged, features = glist[3:5],override_plot_dims=T,pt_size=3) & ThemeLegendRight()
MapFeatures(se_merged, features = glist[6:7],override_plot_dims=T,pt_size=3) & ThemeLegendRight()
MapFeatures(se_merged, features = glist[8:10],override_plot_dims=T,pt_size=3) & ThemeLegendRight()
se_merged1 <- MergeSTData(SCTransform(sp[[1]],assay = "Spatial"), SCTransform(sp[[2]],assay="Spatial"))
MapFeatures(se_merged1, features = glist[3:5],override_plot_dims=T,pt_size=3) & ThemeLegendRight()

se <- se |>
  NormalizeData() |>
  ScaleData() |>
  FindVariableFeatures() |>
  RunPCA() |>
  FindNeighbors(reduction = "pca", dims = 1:10) |>
  FindClusters(resolution = c(0.2,0.5,0.8))
saveRDS(se,file="sp_semla_RD.rds")
p1<-MapLabels(se, column_name = c("Spatial_snn_res.0.2"), ncol = 2,pt_size = 3) & theme(legend.position = "right")
p2<-MapLabels(se, column_name = c("Spatial_snn_res.0.5"), ncol = 2,pt_size = 3) & theme(legend.position = "right")
p3<-MapLabels(se, column_name = c("Spatial_snn_res.0.8"), ncol = 2,pt_size = 3) & theme(legend.position = "right")
pdf("0.spatial_cluster.pdf",width=20,height=10)
print(p1)
print(p2)
print(p3)
dev.off()
pdf("0.spatial_marker.pdf",width=12,height=10)
for( g in glist){
  p<-MapFeatures(se_merged, features = g,override_plot_dims=T,pt_size=3) & ThemeLegendRight()
  print(p)
}
dev.off()

library(dplyr)
library(stringr)
library(cowplot)

wdir<-"/data2/Project/anding/test_20250910/PFC-MERFISH_code/ref_h5ad/output/"
setwd("/data2/Project/anding/test_20250910/PFC-MERFISH_code/downstream")
sp<-readRDS("/data2/Project/anding/test_20250910/PFC-MERFISH_code/merFISH_anno.rds")
test<-readRDS("/data2/Project/anding/test_20250910/subtype.rds")
table(test$sampleinf,test$Group)
Idents(test)<-test$stim
test<-RenameIdents(test,"A"="CON","B"="DEP","D"="KET")
test$Group<-test@active.ident
test<-subset(test,Group %in% c("CON","DEP","KET"))
table(test$sampleinf,test$Group)
#sid<-"20211024_PFCL1"
sid_list<-c("20211121_PFCL1","20211107_PFCL1")
sid<-sid_list[1]
sid<-sid_list[2]
for(sid in unique(sp$orig.ident)){
  message("processing ",sid)
sp2<-subset(sp,orig.ident==sid)
res1<-read.csv(paste0(wdir,sid,"/",sid,"_cells_to_spcoord.txt"),sep="\t")
res1$cxy<-paste0(res1$centroid_1 ,"_",res1$centroid_2)
df1<-sp2@meta.data
df1$cid<-rownames(df1)
df1$cxy<-paste0(df1$centroid_1 ,"_",df1$centroid_2)
rownames(df1)<-df1$cxy
rownames(res1)<-res1$cell
res1$Label1<-df1[res1$cxy,"Label1"]
res1$Label2<-df1[res1$cxy,"Label2"]
test$Label1<-res1[rownames(test@meta.data),"Label1"]
test$Label2<-res1[rownames(test@meta.data),"Label2"]
table(test$SubType,test$Label2)
res1$celltype.stim<-factor(as.character(res1$celltype.stim),
 levels=c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)", "7(Exc5)"  , "4(Exc4)","8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
           "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
      "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)")) 
test$check <- with(test@meta.data, ifelse(
  # Exc + Excitatory
  (str_detect(SubType, "^\\d*\\(Exc") & Label1 == "Excitatory") |
    # Int + Inhibitory
    (str_detect(SubType, "^\\d*\\(Int") & Label1 == "Inhibitory") |
    # 其他 + Non-Neuron
    (!(str_detect(SubType, "^\\d*\\(Exc") | str_detect(SubType, "^\\d*\\(Int")) & Label1 == "Non-Neuron"),
  "pass", "fail"
))
table(test$SubType,test$check)
slt<-subset(test,check=="pass")
saveRDS(slt@meta.data,file=paste0("subtype_",sid,"_emb_meta.rds"))
#slt<-subset(test,Label1=="Non-Neuron",invert=T)
#CellStatPlot(slt,group.by="Group",stat.by="Label2",plot_type="trend")
res1<-res1[rownames(slt@meta.data),]
saveRDS(res1,file=paste0(sid,"_res1.rds"))
p<-ggplot(res1,aes(x=centroid_1,y=centroid_2))+geom_point(aes(color=celltype.stim))+
  theme_scp()+scale_color_manual(values=as.character(SCP::palette_scp(n=24)))
pdf(paste0(sid,"_sptialCoord.pdf"),width=10,height=8)
print(p)
dev.off()
message(sid," p1 done")
stype<-c("4(Exc4)","10(Exc7)","17(Exc9)","18(Exc10)","9(Int2)","13(Mic1)","5(Oli1)","6(Ast1)")
p<-CellStatPlot(slt,group.by="Group",stat.by="Label2",legend.position = "bottom",
                plot_type="trend",split.by="SubType",combine=F)
p<-p[paste0("Group:CON,DEP,KET:",stype)]
p<-cowplot::plot_grid(plotlist=p,nrow=1)
pdf(paste0(sid,"_spLabel2.pdf"),width=16,height=7)
print(p)
dev.off()

res<-c()
for( s in unique(slt$Group)){
  tmp<-subset(slt,Group==s)
  df<-as.data.frame.matrix(prop.table(table(tmp$Label2,tmp$SubType),2)*100)
  df$Group<-s
  df$Labels<-rownames(df)
  if(s == unique(slt$sampleinf)[1]){res<-df}
  else(res<-rbind(res,df))
}
write.xlsx(res,file=paste0(sid,"_subtype_sp_ratio.xlsx"))
message(sid, "all done")
}

test<-readRDS("/data2/Project/anding/test_20250910/subtype.rds")
Idents(test)<-test$stim
test<-RenameIdents(test,"A"="CON","B"="DEP","D"="KET")
test$Group<-test@active.ident
test<-subset(test,Group %in% c("CON","DEP","KET"))
df<-GetAssayData(test)
for(sid in unique(sp$orig.ident)){
  sp2<-subset(sp,orig.ident==sid)
  res1<-read.csv(paste0(wdir,sid,"/",sid,"_cells_to_spcoord.txt"),sep="\t")
  res1$cxy<-paste0(res1$centroid_1 ,"_",res1$centroid_2)
  df1<-sp2@meta.data
  df1$cid<-rownames(df1)
  df1$cxy<-paste0(df1$centroid_1 ,"_",df1$centroid_2)
  rownames(df1)<-df1$cxy
  rownames(res1)<-res1$cell
  res1$Label1<-df1[res1$cxy,"Label1"]
  res1$Label2<-df1[res1$cxy,"Label2"]
  test$Label1<-res1[rownames(test@meta.data),"Label1"]
  test$Label2<-res1[rownames(test@meta.data),"Label2"]
  res1$celltype.stim<-factor(as.character(res1$celltype.stim),
                             levels=c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)", "7(Exc5)"  , "4(Exc4)","8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
                                      "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
                                      "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)")) 
gint2<-toupper(c("Tenm2", 'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1'))
res1<-res1[rownames(res1) %in% colnames(df),]
res1[,gint2]<-df[gint2,rownames(res1)]
p<-list()
for(g in gint2){
  res_sorted <- res1 %>% arrange(.data[[g]]) #res_sorted <- res1 %>% arrange(desc(.data[[g]]))
  p[[g]]<-ggplot(res_sorted, aes(x = centroid_1, y = centroid_2, color = .data[[g]])) +
    geom_point(size=0.5) +
    scale_color_gradient2(low = "snow", high = "red") +
    theme_scp()
}
p<-cowplot::plot_grid(plotlist=p,ncol=5)
pdf(paste0(sid,"_gs1.pdf"),width=30,height=14)
print(p)
dev.off()
ast1<-toupper(c('Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2'))
res1[,ast1]<-df[ast1,rownames(res1)]
p2<-list()
for(g in ast1){
  res_sorted <- res1 %>% arrange(.data[[g]]) #res_sorted <- res1 %>% arrange(desc(.data[[g]]))
  p2[[g]]<-ggplot(res_sorted, aes(x = centroid_1, y = centroid_2, color = .data[[g]])) +
    geom_point(size=0.5) +
    scale_color_gradient2(low = "snow", high = "red") +
    theme_scp()
}
p<-cowplot::plot_grid(plotlist=p2,ncol=5)
pdf(paste0(sid,"_gs2.pdf"),width=30,height=14)
print(p)
dev.off()
gs3<-toupper(c("App","Lrrc4c","Ncam1","Nfasc","Tenm2",
               "Sema6a",
               "Ptn",
               "Nrxn3",
               'Nrxn1',
               "Ncam2"))
res1[,gs3]<-df[gs3,rownames(res1)]
p3<-list()
for(g in gs3){
  res_sorted <- res1 %>% arrange(.data[[g]]) #res_sorted <- res1 %>% arrange(desc(.data[[g]]))
  p3[[g]]<-ggplot(res_sorted, aes(x = centroid_1, y = centroid_2, color = .data[[g]])) +
    geom_point(size=0.5) +
    scale_color_gradient2(low = "snow", high = "red") +
    theme_scp()
}
p<-cowplot::plot_grid(plotlist=p3,ncol=5)
pdf(paste0(sid,"_gs3.pdf"),width=30,height=14)
print(p)
dev.off()
gs4<-toupper(c('Lrrc4c',
               'Nrxn3',
               'Nrxn1',
               'Cadm1',
               'ERBB4',
               'NRG3',
               'NLGN1',
               'CNTNAP2',
               'TENM2',
               'MEF2A'))
res1[,gs4]<-df[gs3,rownames(res1)]
p4<-list()
for(g in gs4){
  res_sorted <- res1 %>% arrange(.data[[g]]) #res_sorted <- res1 %>% arrange(desc(.data[[g]]))
  p4[[g]]<-ggplot(res_sorted, aes(x = centroid_1, y = centroid_2, color = .data[[g]])) +
    geom_point(size=0.5) +
    scale_color_gradient2(low = "snow", high = "red") +
    theme_scp()
}
p<-cowplot::plot_grid(plotlist=p4,ncol=5)
pdf(paste0(sid,"_gs4.pdf"),width=30,height=14)
print(p)
dev.off()
}

dat<-read.table("../ref_h5ad/output/20211024_PFCL1/20211024_PFCL1_scatter_data.csv",sep="\t",header=T)
dim(dat)
head(dat)
dat$celltype.stim<-factor(as.character(dat$celltype.stim),
                           levels=c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)", "7(Exc5)"  , "4(Exc4)","8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
                                    "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
                                    "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)")) 
ggplot(dat,aes(x=centroid_1,y=centroid_2))+
  geom_point(aes(color=celltype.stim))+
  theme_scp(size=0.5)+
  scale_color_manual(values=as.character(SCP::palette_scp(n=24)))

setwd("20250921")
library(dplyr)
library(scatterpie)
#sid<-sid_list[1]
res1<-readRDS(paste0("../downstream/",sid,"_res1.rds"))
# 1. 按 cxy 统计 count 和 ratio
df_summary <- res1 %>%
  group_by(cxy, centroid_1, centroid_2, celltype.stim) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cxy, centroid_1, centroid_2) %>%
  mutate(ratio = count / sum(count))
df_summary$celltype.stim<-factor(as.character(df_summary$celltype.stim),
                                 levels=c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)", "7(Exc5)"  , "4(Exc4)","8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
                                          "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
                                          "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)"))

# 2. 转换为 wide 格式，方便绘制饼图
df_wide <- df_summary %>%
  select(cxy, centroid_1, centroid_2, celltype.stim, ratio) %>%
  tidyr::pivot_wider(names_from = celltype.stim, values_from = ratio, values_fill = 0)

# 3. 绘制饼图散点图
# 你定义好的 levels
ct_levels <- c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)","4(Exc4)", "7(Exc5)"  , "8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
               "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
               "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)")
#11+8+5=24
# 让 df_wide 的列顺序和 levels 对应
cols_for_plot <- intersect(ct_levels, colnames(df_wide))

p<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 60),
    data = df_wide,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
pdf(paste0(sid,"_scatterpie.pdf"),width=50,height=30)
print(p)
dev.off()

ggplot()+  geom_scatterpie(
  aes(x = centroid_1, y = centroid_2, group = cxy, r = 60),
  data = df_wide,
  cols = cols_for_plot,  # 按 levels 顺序传递
  color = NA
) +
  coord_equal()+geom_vline(xintercept =c(-3000,-500,3000) )

df1<-df_wide[df_wide$centroid_1< -3000,]
df2<-df_wide[df_wide$centroid_1> -3000 & df_wide$centroid_1 < -500,]
df3<-df_wide[df_wide$centroid_1> -500 & df_wide$centroid_1<3000,]
df4<-df_wide[df_wide$centroid_1> 5000 ,]
p1<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 30),
    data = df1,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p1  
p2<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 30),
    data = df2,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p2
p3<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 30),
    data = df3,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p3
p4<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 30),
    data = df4,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p4
pdf(paste0(sid,"_scatterpie.pdf"),width=50,height=30)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

#sid2
df1<-df_wide[df_wide$centroid_1< -2000 & df_wide$centroid_2>0,]
df2<-df_wide[df_wide$centroid_1< -2000 & df_wide$centroid_2 <0,]
df3<-df_wide[df_wide$centroid_1> 0 & df_wide$centroid_1<5000,]
df4<-df_wide[df_wide$centroid_1> 5000 ,]
p1<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 20),
    data = df1,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p1  
p2<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 20),
    data = df2,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p2
p3<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 20),
    data = df3,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p3
p4<-ggplot() +
  geom_scatterpie(
    aes(x = centroid_1, y = centroid_2, group = cxy, r = 20),
    data = df4,
    cols = cols_for_plot,  # 按 levels 顺序传递
    color = NA
  ) +
  coord_equal() +
  theme_scp()+scale_fill_manual(values=as.character(SCP::palette_scp(n=24)))
p4
pdf(paste0(sid,"_scatterpie.pdf"),width=50,height=40)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

for(sid in sid_list){
meta<-readRDS(paste0("subtype_",sid,"_emb_meta.rds"))
res<-c()
for( mt in unique(meta$Label1)){
  tmp<-meta[meta$Label1==mt,]
  df<-as.data.frame.matrix(prop.table(table(tmp$Label2,tmp$sampleinf),2)*100)
  df$Labels<-rownames(df)
  if(s == unique(slt$sampleinf)[1]){res<-df}
  else(res<-rbind(res,df))
}
res
write.xlsx(res,file=paste0(sid,"_sampleRatio.xlsx"))
}

for( sid in sid_list){
meta<-read.xlsx(paste0(sid,"_sampleRatio.xlsx"))
plt<-melt(meta)
anno<-unique(data.frame(sample=test$sampleinf,group=test$Group))
rownames(anno)<-anno$sample
plt$Group<-anno[plt$variable,"group"]
p<-ggplot(plt, aes(x = Labels, y = value, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) + 
  geom_point(aes(color=Group),position=position_dodge(width=0.8))+
  ylim(c(0, 50)) +
  stat_compare_means(
    method = "t.test", 
    label = "p.signif", 
    label.y = 48, 
    hide.ns = TRUE, 
    size = 8
  ) +
  theme_scp() +theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
  labs(x = "", y = "Proportion")+
  scale_color_manual(values=as.character(SCP::palette_scp(n=3)))+
  scale_fill_manual(values=as.character(SCP::palette_scp(n=3)))
pdf(paste0(sid,"_sampleRatio.pdf"),width=8,height=4)
print(p)
dev.off()
}
# res2<-compare_means(value~Group,group.by = "Labels",data=plt,method="anova")

cluster2type <- setNames(mtx$celltype2, mtx$cluster)
# 给 Seurat 对象加一列 celltype2
dat$celltype2 <- as.character(cluster2type[as.character(dat$seurat_clusters)])
# 如果想把 celltype2 作为新的分群（Idents）
Idents(dat) <- "celltype2"
CellDimPlot(dat,group.by="celltype2",label=T,label_insitu = T,label_repel = T,legend.position = "none")
saveRDS(dat,file="all.combined_celltype2.rds")


library(homologene)
human_genes <- rownames(dat)
map_df <- homologene(human_genes, inTax = 9606, outTax = 10090)
# 转换成人-鼠同源基因
map_df <- homologene(human_genes, inTax = 9606, outTax = 10090)
head(map_df)
map_df <- na.omit(map_df)
colnames(map_df) <- c("human", "mouse")
seurat_obj<-dat
convertHumanToMouse <- function(seurat_obj, assay = "RNA", 
                                inTax = 9606, outTax = 10090, 
                                method = c("first", "sum")) {
  # 加载依赖
  if (!requireNamespace("homologene", quietly = TRUE)) {
    stop("请先安装 homologene: remotes::install_github('oganm/homologene')")
  }
  method <- match.arg(method)
  
  message("提取基因名...")
  human_genes <- rownames(seurat_obj[[assay]]@counts)
  
  message("查询 homologene 映射...")
  map_df <- homologene::homologene(human_genes, inTax = inTax, outTax = outTax)
  map_df <- na.omit(map_df)
  colnames(map_df) <- c("human", "mouse")
  
  if (nrow(map_df) == 0) stop("没有找到任何基因映射！")
  
  # 处理多对一情况
  if (method == "first") {
    map_df <- map_df[!duplicated(map_df$human), ]
    gene_map <- setNames(map_df$mouse, map_df$human)
  } else if (method == "sum") {
    # 多个人类基因对应同一鼠基因 → 合并表达量
    library(Matrix)
    message("合并多对一映射的基因表达矩阵...")
    
    gene_map <- setNames(map_df$mouse, map_df$human)
    counts <- seurat_obj[[assay]]@counts
    common_genes <- intersect(rownames(counts), names(gene_map))
    
    new_counts <- Matrix::aggregate.Matrix(counts[common_genes, ], 
                                           groupings = gene_map[common_genes], 
                                           fun = "sum")
    new_data <- Matrix::aggregate.Matrix(seurat_obj[[assay]]@data[common_genes, ], 
                                         groupings = gene_map[common_genes], 
                                         fun = "sum")
    
    new_scale <- seurat_obj[[assay]]@scale.data
    if (nrow(new_scale) > 0) {
      common_scale <- intersect(rownames(new_scale), names(gene_map))
      new_scale <- Matrix::aggregate.Matrix(new_scale[common_scale, ], 
                                            groupings = gene_map[common_scale], 
                                            fun = "mean")
    }
    
    seurat_obj[[assay]]@counts <- new_counts
    seurat_obj[[assay]]@data <- new_data
    seurat_obj[[assay]]@scale.data <- new_scale
    message("转换完成 ✅")
    return(seurat_obj)
  }
  
  # === method = "first" 情况 ===
  message("应用一对一映射...")
  gene_map <- setNames(map_df$mouse, map_df$human)
  
  counts <- seurat_obj[[assay]]@counts
  data <- seurat_obj[[assay]]@data
  scale <- seurat_obj[[assay]]@scale.data
  
  common_genes <- intersect(rownames(counts), names(gene_map))
  
  new_counts <- counts[common_genes, ]
  rownames(new_counts) <- gene_map[rownames(new_counts)]
  
  new_data <- data[common_genes, ]
  rownames(new_data) <- gene_map[rownames(new_data)]
  
  if (nrow(scale) > 0) {
    common_scale <- intersect(rownames(scale), names(gene_map))
    new_scale <- scale[common_scale, ]
    rownames(new_scale) <- gene_map[rownames(new_scale)]
  } else {
    new_scale <- scale
  }
  
  seurat_obj[[assay]]@counts <- new_counts
  seurat_obj[[assay]]@data <- new_data
  seurat_obj[[assay]]@scale.data <- new_scale
  
  message("转换完成 ✅")
  return(seurat_obj)
}
#usage
# 默认：人类 (9606) → 小鼠 (10090)，取第一个映射
dat <- convertHumanToMouse(dat, assay = "RNA", method = "first")
# 如果想合并多个对应的人类基因 → 一个小鼠基因
mouse_obj <- convertHumanToMouse(seurat_obj, assay = "RNA", method = "sum")
dir.create("ref_h5ad")
setwd("ref_h5ad")
test<-readRDS("all.combined_celltype2_mouseid.rds")
adata<-srt_to_adata(test)
adata$write_h5ad("subtype.h5ad")
for( s in unique(test$Group)){
  tmp<-subset(test,Group==s)
  adata<-srt_to_adata(tmp)
  adata$write_h5ad(paste0(s,".h5ad"))
}
register(MulticoreParam(workers = length(unique(test$celltype2)), progressbar = TRUE))
test<-RunDEtest(test,group_by="celltype2",only.pos = T,fc.threshold = 1)
deg<-test@tools$DEtest_celltype2$AllMarkers_wilcox
deg<-deg[deg$p_val_adj<0.05,]
write.csv(data.frame(gene=unique(deg$gene)),file="seurat_cluster_deg.csv")
for( s in unique(test$Group)){
  tmp<-subset(test,Group==s)
  tmp<-RunDEtest(tmp,group_by="celltype2",only.pos = T,fc.threshold = 1)
  deg<-tmp@tools$DEtest_celltype2$AllMarkers_wilcox
  deg<-deg[deg$p_val_adj<0.05,]
  write.csv(data.frame(gene=unique(deg$gene)),file=paste0(s,"_seurat_cluster_deg.csv"))
}





#CON_output
wdir<-"/data2/Project/anding/test_20250923/"
test<-readRDS("/data2/Project/anding/test_20250923/all.combined_celltype2_mouseid.rds")
sp<-readRDS("/data2/Project/anding/test_20250910/PFC-MERFISH_code/merFISH_anno.rds")
Idents(test)<-test$stim
df<-GetAssayData(test)
# dir.create("sep_downstream")
# setwd("sep_downstream")
for(sid in unique(sp$orig.ident)){
  sp2<-subset(sp,orig.ident==sid)
  res1<-c()
  for( type in c("CON","DEP","KET")){
    message("reading ",sid," ",type)
  tmp<-read.csv(paste0(wdir,type,"_output/",sid,"/",sid,"_cells_to_spcoord.txt"),sep="\t")
  tmp$cxy<-paste0(tmp$centroid_1 ,"_",tmp$centroid_2)
  df1<-sp2@meta.data
  df1$cid<-rownames(df1)
  df1$cxy<-paste0(df1$centroid_1 ,"_",df1$centroid_2)
  rownames(df1)<-df1$cxy
  rownames(tmp)<-tmp$cell
  tmp$Label1<-df1[tmp$cxy,"Label1"]
  tmp$Label2<-df1[tmp$cxy,"Label2"]
  if(type=="CON"){res1<-tmp}
  else{res1<-rbind(res1,tmp)}
  }
  test$Label1<-res1[rownames(test@meta.data),"Label1"]
  test$Label2<-res1[rownames(test@meta.data),"Label2"]
  res1$celltype.stim<-factor(as.character(res1$celltype.stim),
                             levels=c("0(Exc1)" ,"1(Exc2)" ,"3(Exc3)", "7(Exc5)"  , "4(Exc4)","8(Exc6)", "10(Exc7)", "15(Exc8)" , "17(Exc9)"  , "18(Exc10)",  "20(Exc11)" ,
                                      "2(Int1)", "9(Int2)",  "11(Int3)" ,"12(Int4)" ,   "16(Int5)", "19(Int6)",  "21(Int7)",  "23(Int8)" ,    
                                      "6(Ast1)" ,  "22(Ast2)",  "5(Oli1)"  ,"13(Mic1)"  ,"14(OPC1)")) 
  test$check <- with(test@meta.data, ifelse(
    # Exc + Excitatory
    (str_detect(celltype2 , "^\\d*\\(Exc") & Label1 == "Excitatory") |
      # Int + Inhibitory
      (str_detect(celltype2 , "^\\d*\\(Int") & Label1 == "Inhibitory") |
      # 其他 + Non-Neuron
      (!(str_detect(celltype2 , "^\\d*\\(Exc") | str_detect(celltype2 , "^\\d*\\(Int")) & Label1 == "Non-Neuron"),
    "pass", "fail"
  ))
  table(test$Label1,test$check)
  slt<-subset(test,check=="pass")
  saveRDS(slt@meta.data,file=paste0("subtype_",sid,"_emb_meta.rds"))
  #slt<-subset(test,Label1=="Non-Neuron",invert=T)
  #CellStatPlot(slt,group.by="Group",stat.by="Label2",plot_type="trend")
  res1<-res1[rownames(slt@meta.data),]
  saveRDS(res1,file=paste0(sid,"_res1.rds"))
  p<-ggplot(res1,aes(x=centroid_1,y=centroid_2))+geom_point(aes(color=celltype.stim))+
    theme_scp()+scale_color_manual(values=as.character(SCP::palette_scp(n=24)))
  pdf(paste0(sid,"_sptialCoord.pdf"),width=10,height=8)
  print(p)
  dev.off()
  message(sid," p1 done")
  stype<-c("4(Exc4)","10(Exc7)","17(Exc9)","18(Exc10)","9(Int2)","13(Mic1)","5(Oli1)","6(Ast1)")
  p<-CellStatPlot(slt,group.by="Group",stat.by="Label2",legend.position = "bottom",
                  plot_type="trend",split.by="celltype2",combine=F)
  p<-p[paste0("Group:CON,DEP,KET:",stype)]
  p<-cowplot::plot_grid(plotlist=p,nrow=1)
  pdf(paste0(sid,"_spLabel2.pdf"),width=16,height=7)
  print(p)
  dev.off()
  res<-c()
  slt$celltype2<-factor(as.character(slt$celltype2),levels=unique(slt$celltype2))
  for( s in unique(slt$Group)){
    tmp<-subset(slt,Group==s)
    tmp$celltype2<-factor(as.character(tmp$celltype2),levels=levels(slt$celltype2))
    df<-as.data.frame.matrix(prop.table(table(tmp$Label2,tmp$celltype2),2)*100)
    df$Group<-s
    df$Labels<-rownames(df)
    if(s == unique(slt$Group)[1]){res<-df}
    else(res<-rbind(res,df))
  }
  write.xlsx(res,file=paste0(sid,"_sp_ratio.xlsx"))
  message(sid," all done")
  
  for( s in unique(slt$sampleinf)){
    tmp<-subset(slt,sampleinf==s)
    tmp$celltype2<-factor(as.character(tmp$celltype2),levels=levels(slt$celltype2))
    df<-as.data.frame.matrix(table(tmp$Label2,tmp$celltype2))
    df$Sample<-s
    df$Labels<-rownames(df)
    if(s == unique(slt$sampleinf)[1]){res<-df}
    else(res<-rbind(res,df))
  }
  write.xlsx(res,file=paste0(sid,"_sample_nCells.xlsx"))
  message(sid," done all")
}
for( fh in list.files(path="./",pattern="_emb_meta.rds")){
  slt<-readRDS(fh)
  sid<-paste0(str_split(fh,"_")[[1]][2],"_",str_split(fh,"_")[[1]][3])
  ratio_plot(slt,sample.by="sampleinf",anno.by="celltype2",condition.by="Group",
               strip.col=NULL,save.prefix=sid,width=14,height=18)
}


#20220128和20220218
#1）8种细胞每一种细胞在大脑前额叶的分布，加上大脑皮层6层的大致位置；这样就8张图，目的是证实这8个亚群的存在和分布；
#2）8种细胞在一起的分布彩图，目的是证明这8种细胞互作是有形态基础的，请标一下大脑六层细胞的大致位置；
#3）在DEP和KET处理后的细胞或基因变化；
#4）关键基因的表达分布图


library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(stringr)
library(SCP)       # 你代码里用了 palette_scp
library(ggnewscale)

# 工作路径 & 数据
wdir <- "/data2/Project/anding/test_20250923/"
test <- readRDS("/data2/Project/anding/test_20250923/all.combined_celltype2_mouseid.rds")
sp <- readRDS("/data2/Project/anding/test_20250910/PFC-MERFISH_code/merFISH_anno.rds")

Idents(test) <- test$stim
# 目标样本
target_samples <- c("20220128_PFCL1","20220218_PFCL1")
# 8种细胞类型
stype <- c("4(Exc4)","10(Exc7)","17(Exc9)","18(Exc10)",
           "9(Int2)","13(Mic1)","5(Oli1)","6(Ast1)")

for(sid in target_samples){
  message("Processing ", sid)
  
  sp2 <- subset(sp, orig.ident == sid)
  res1 <- list()
  
  for(type in c("CON","DEP","KET")){
    fn <- file.path(wdir, paste0(type,"_output"), sid, paste0(sid,"_cells_to_spcoord.txt"))
    if(!file.exists(fn)) next
    tmp <- read.table(fn, sep="\t", header=TRUE)
    tmp$cxy <- paste0(tmp$centroid_1, "_", tmp$centroid_2)
    
    df1 <- sp2@meta.data
    df1$cid <- rownames(df1)
    df1$cxy <- paste0(df1$centroid_1, "_", df1$centroid_2)
    rownames(df1) <- df1$cxy
    rownames(tmp) <- tmp$cell
    
    tmp$Label1 <- df1[tmp$cxy,"Label1"]
    tmp$Label2 <- df1[tmp$cxy,"Label2"]
    tmp$Group  <- type
    tmp$cid<-rownames(tmp)
    res1[[type]] <- tmp
  }
  
  res1 <- do.call(rbind, res1)
  rownames(res1)<-res1$cid
  
  # check 匹配
  test$Label1 <- res1[rownames(test@meta.data),"Label1"]
  test$Label2 <- res1[rownames(test@meta.data),"Label2"]
  
  test$check <- with(test@meta.data, ifelse(
    (str_detect(celltype2 , "^\\d*\\(Exc") & Label1 == "Excitatory") |
      (str_detect(celltype2 , "^\\d*\\(Int") & Label1 == "Inhibitory") |
      (!(str_detect(celltype2 , "^\\d*\\(Exc") | str_detect(celltype2 , "^\\d*\\(Int")) & Label1 == "Non-Neuron"),
    "pass","fail"
  ))
  
  slt <- subset(test, check == "pass")
  saveRDS(slt, file=paste0("subtype_", sid, "_meta.rds"))
  saveRDS(res1, file=paste0(sid,"_res1.rds"))
}
#20220218_PFCL1 x.cut=-2000
# 伪造皮层层次标注，实际请根据MERFISH坐标调整
library(ggplot2)

# 假设 res1 中有 celltype 列，stype 是 celltype 向量
# 伪造皮层层次标注，实际请根据MERFISH坐标调整
library(ggplot2)
library(dplyr)
# 指定需要标注的皮层层类型
layer_types <- c("Exc L2/3 IT", "Exc L4/5 IT", "Exc L5 IT", "Exc L5/6 NP")
# 颜色可以自定义，也可以用 palette_scp
layer_colors <- c(
  "Exc L2/3 IT" = "#FFD700",   # 金色
  "Exc L4/5 IT" = "#91ccae",   # 橙色
  "Exc L5 IT"   = "#f6c6d6",   # 深红
  "Exc L5/6 NP" = "#6A5ACD"    # 蓝紫
)
e4<-c("#f5ad65","#91ccae","#795291","#f6c6d6")

dir.create("p1", showWarnings = FALSE)
sid<-"20220218_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data),]
res_split1<-res1[res1$centroid_1 < -2000,]
res_split2<-res1[res1$centroid_1 > -2000,]
res1<-res_split1
# 计算每个皮层层的凸包，用于背景框
  layer_polygons <- res1 %>%
    dplyr::filter(Label2 %in% layer_types) %>%
    dplyr::group_by(Label2) %>%
    dplyr::slice(chull(centroid_1, centroid_2))  # 计算凸包边界点
  
  for (ct in stype) {
    # 高亮当前细胞类型，其余为灰色
    res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
    res1 <- res1[order(res1$plot_color == ct), ]
    
    cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                              "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
    color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
    names(color_values)[2] <- ct
    
    p <- ggplot() +
      # 背景层凸包区域
      #geom_polygon(
      #  data = layer_polygons,
      #  aes(x = centroid_1, y = centroid_2, fill = Label2, group = Label2),
      #  alpha = 0.15,  # 背景透明度
      #  color = NA
      #) +
      scale_fill_manual(values = layer_colors) +
      # 细胞点
      geom_point(
        data = res1,
        aes(x = centroid_1, y = centroid_2, color = plot_color),
        size = 0.2, alpha = 0.8
      )+ 
      scale_color_manual(values = color_values) +
      labs(title = paste0(sid, " - ", ct)) +
      theme_scp() +
      coord_fixed() +
      guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
    p+facet_wrap(~Label2)
    
    pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_left.pdf"), width = 14, height = 14)
    print(p+facet_wrap(~Label2))
    dev.off()
  }
res1<-res_split2
  for (ct in stype) {
    # 高亮当前细胞类型，其余为灰色
    res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
    res1 <- res1[order(res1$plot_color == ct), ]
    
    cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                              "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
    color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
    names(color_values)[2] <- ct
    
    p <- ggplot() +
      # 背景层凸包区域
      #geom_polygon(
      #  data = layer_polygons,
      #  aes(x = centroid_1, y = centroid_2, fill = Label2, group = Label2),
      #  alpha = 0.15,  # 背景透明度
      #  color = NA
      #) +
      scale_fill_manual(values = layer_colors) +
      # 细胞点
      geom_point(
        data = res1,
        aes(x = centroid_1, y = centroid_2, color = plot_color),
        size = 0.2, alpha = 0.8
      )+ 
      scale_color_manual(values = color_values) +
      labs(title = paste0(sid, " - ", ct)) +
      theme_scp() +
      coord_fixed() +
      guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
    p+facet_wrap(~Label2)
    
    pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_right.pdf"), width = 14, height = 14)
    print(p+facet_wrap(~Label2))
    dev.off()
  }  

sid<-"20220128_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 & res1$centroid_2 > -1000,]
res<-res_split[[1]]
plist<-list()
stype<-c("4(Exc4)","10(Exc7)","17(Exc9)","18(Exc10)","9(Int2)","13(Mic1)","5(Oli1)","6(Ast1)")
for (ct in stype[1:4]) {
  res1<-res
  # 高亮当前细胞类型，其余为灰色
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  res1<-res1[res1$Label2 %in% c("Exc L2/3 IT", "Exc L4/5 IT", "Exc L5 IT", "Exc L5/6 NP"),]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  facet_list <- lapply(labels, function(lbl) {
    df <- res1[res1$Label2 == lbl, ]
    ggplot(df, aes(x = centroid_1, y = centroid_2, color = plot_color)) +
      geom_point(size = 0.2, alpha = 0.8) +
      scale_color_manual(values = color_values) +
      labs(fill = "type", x = "", y = sub(".*\\(([^)]+)\\).*", "\\1", ct)) +
      theme_scp() +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
      coord_fixed() +
      guides(color = "none") +
      ggtitle(lbl)  # 给每个 facet 一个标题
  })
  
  # 将这个细胞类型的 facet 列表存入总列表
  plist[[ct]] <- facet_list
}  
for (ct in stype[5]) {
  res1<-res
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  res1<-res1[res1$Label2 %in% c("Int Lamp5", "Int Pvalb", "Int Sst", "Int GABAergic"),]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  facet_list <- lapply(labels, function(lbl) {
    df <- res1[res1$Label2 == lbl, ]
    ggplot(df, aes(x = centroid_1, y = centroid_2, color = plot_color)) +
      geom_point(size = 0.2, alpha = 0.8) +
      scale_color_manual(values = color_values) +
      labs(fill = "type", x = "", y = sub(".*\\(([^)]+)\\).*", "\\1", ct)) +
      theme_scp() +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
      coord_fixed() +
      guides(color = "none") +
      ggtitle(lbl)  # 给每个 facet 一个标题
  })
  
  # 将这个细胞类型的 facet 列表存入总列表
  plist[[ct]] <- facet_list
}  
for (ct in stype[-c(1:5)]) {
  res1<-res
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  res1<-res1[res1$Label2 %in% c( "Non-Neuron Astrocytes", "Non-Neuron Endo","Non-Neuron Oligo", "Non-Neuron Microglia","Non-Neuron OPC"),]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  facet_list <- lapply(labels, function(lbl) {
    df <- res1[res1$Label2 == lbl, ]
    ggplot(df, aes(x = centroid_1, y = centroid_2, color = plot_color)) +
      geom_point(size = 0.2, alpha = 0.8) +
      scale_color_manual(values = color_values) +
      labs(fill = "type", x = "", y = sub(".*\\(([^)]+)\\).*", "\\1", ct)) +
      theme_scp() +
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
      coord_fixed() +
      guides(color = "none") +
      ggtitle(lbl)  # 给每个 facet 一个标题
  })
  
  # 将这个细胞类型的 facet 列表存入总列表
  plist[[ct]] <- facet_list
}  
plot_grid(plist$`9(Int2)`[[1]],plist$`9(Int2)`[[2]],plist$`9(Int2)`[[3]],plist$`9(Int2)`[[3]],
  plist$`4(Exc4)`[[1]],plist$`4(Exc4)`[[3]],plist$`4(Exc4)`[[4]],plist$`4(Exc4)`[[2]],ncol=4,align = "vh",axis="l")



res1<-res_split2
for (ct in stype) {
  # 高亮当前细胞类型，其余为灰色
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  p <- ggplot() +
    scale_fill_manual(values = layer_colors) +
    # 细胞点
    geom_point(
      data = res1,
      aes(x = centroid_1, y = centroid_2, color = plot_color),
      size = 0.2, alpha = 0.8
    )+ 
    scale_color_manual(values = color_values) +
    labs(title = paste0(sid, " - ", ct)) +
    theme_scp() +
    coord_fixed() +
    guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
  p+facet_wrap(~Label2)
  
  pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_region2.pdf"), width = 14, height = 14)
  print(p+facet_wrap(~Label2))
  dev.off()
} 
res1<-res_split3
for (ct in stype) {
  # 高亮当前细胞类型，其余为灰色
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  p <- ggplot() +
    scale_fill_manual(values = layer_colors) +
    # 细胞点
    geom_point(
      data = res1,
      aes(x = centroid_1, y = centroid_2, color = plot_color),
      size = 0.2, alpha = 0.8
    )+ 
    scale_color_manual(values = color_values) +
    labs(title = paste0(sid, " - ", ct)) +
    theme_scp() +
    coord_fixed() +
    guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
  p+facet_wrap(~Label2)
  
  pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_region3.pdf"), width = 14, height = 14)
  print(p+facet_wrap(~Label2))
  dev.off()
} 
res1<-res_split4
for (ct in stype) {
  # 高亮当前细胞类型，其余为灰色
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  p <- ggplot() +
    scale_fill_manual(values = layer_colors) +
    # 细胞点
    geom_point(
      data = res1,
      aes(x = centroid_1, y = centroid_2, color = plot_color),
      size = 0.2, alpha = 0.8
    )+ 
    scale_color_manual(values = color_values) +
    labs(title = paste0(sid, " - ", ct)) +
    theme_scp() +
    coord_fixed() +
    guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
  p+facet_wrap(~Label2)
  
  pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_region4.pdf"), width = 14, height = 14)
  print(p+facet_wrap(~Label2))
  dev.off()
} 
res1<-res_split5
for (ct in stype) {
  # 高亮当前细胞类型，其余为灰色
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                            "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  p <- ggplot() +
    scale_fill_manual(values = layer_colors) +
    # 细胞点
    geom_point(
      data = res1,
      aes(x = centroid_1, y = centroid_2, color = plot_color),
      size = 0.2, alpha = 0.8
    )+ 
    scale_color_manual(values = color_values) +
    labs(title = paste0(sid, " - ", ct)) +
    theme_scp() +
    coord_fixed() +
    guides(color = "none")  # 如果你想保留图例，可以删掉 fill 的 "none"
  p+facet_wrap(~Label2)
  
  pdf(paste0("p1/", sid, "_", ct, "_spatial_single_splitLayer_region5.pdf"), width = 14, height = 14)
  print(p+facet_wrap(~Label2))
  dev.off()
} 

dir.create("p2")
for(sid in target_samples){
  slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
  res1 <- readRDS(paste0(sid, "_res1.rds"))
  res1<-res1[rownames(res1) %in% rownames(slt@meta.data),]
  
  p <- ggplot(subset(res1, celltype.stim %in% stype),
              aes(x=centroid_1, y=centroid_2, color=celltype.stim)) +
    geom_point(size=0.6, alpha=0.7) +
    scale_color_manual(values=cell_colors) +
    theme_scp() + coord_fixed()
  p
  ggsave(paste0("p2/",sid,"_8celltypes_overlay.pdf"), p, width=14, height=8)
}

dir.create("p3")
for(sid in target_samples){
  slt <- readRDS(paste0("subtype_", sid, "_emb_meta.rds"))
  
  p <- CellStatPlot(slt,
                    group.by="Group",
                    stat.by="Label2",
                    legend.position="bottom",
                    plot_type="trend",
                    split.by="celltype2",
                    combine=FALSE)
  p <- p[paste0("Group:CON,DEP,KET:",stype)]
  p <- cowplot::plot_grid(plotlist=p, nrow=1)
  ggsave(paste0(sid,"_DEP_KET_trend.pdf"), p, width=16, height=7)
}
dir.create("p4")
key_genes <- c("Tenm2",'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1',
               'Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2',
               'App', 'Lrrc4c', 'Ncam1','Nfasc','Tenm2', 'Sema6a','Ptn','Nrxn3','Nrxn1','Ncam2',
               'Lrrc4c','Nrxn3','Nrxn1','Cadm1','Erbb4','Nrg3','Nlgn1','Cntnap2','Tenm2','Mef2a')

for( sid in c("20220128_PFCL1","20220218_PFCL1")){
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data),]
#mtx<-GetAssayData(slt,features=g)
cells_use<-rownames(slt@meta.data)
coords <- res1[cells_use, c("centroid_1", "centroid_2")]
coords <- as.matrix(coords)
colnames(coords) <- c("spatial_1", "spatial_2")
rownames(coords) <- cells_use
# 构建 DimReduc 对象
spatial_dr <- CreateDimReducObject(
  embeddings = coords,
  key = "SPATIAL_",
  assay = DefaultAssay(slt)
)
# 写入 slt@reductions
slt[["spatial"]] <- spatial_dr
saveRDS(slt,file=paste0("subtype_",sid,"_spRD.rds"))
key_genes<-key_genes[key_genes %in% rownames(slt)]
for(g in key_genes){
   #p <- FeatureDimPlot(slt,reduction="spatial",features=g,aspect.ratio = 0.5)
   #ggsave(paste0("p4/",sid,"_gene_",g,".pdf"), p, width=14, height=8)
   p1<-FeatureStatPlot(slt,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                    comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
   ggsave(paste0("p4/",sid,"_gene_",g,"_wilcox.pdf"), p1, width=5, height=4)
}
}
int2 <- c("Tenm2",'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1')
Ast1<-c('Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2')	
Oli1<-c('App', 'Lrrc4c', 'Ncam1','Nfasc','Tenm2', 'Sema6a','Ptn','Nrxn3','Nrxn1','Ncam2')	
Mic1<-c('Lrrc4c','Nrxn3','Nrxn1','Cadm1','Erbb4','Nrg3','Nlgn1','Cntnap2','Tenm2','Mef2a')
for( sid in c("20220128_PFCL1","20220218_PFCL1")){
  slt<-readRDS(paste0("subtype_",sid,"_spRD.rds"))
  d1<-subset(slt,celltype2=="9(Int2)")
  int2<-int2[int2 %in% rownames(d1)]
  for( g in int2){
  p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                     comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
  ggsave(paste0("p4/",sid,"_gene_",g,"_Int2Only_wilcox.pdf"), p, width=5, height=4)
  }
  d1<-subset(slt,celltype2=="6(Ast1)")
  Ast1<-Ast1[Ast1 %in% rownames(d1)]
  for(g in Ast1){
  p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                     comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
  ggsave(paste0("p4/",sid,"_gene_",g,"_Ast1Only_wilcox.pdf"), p, width=5, height=4)
  }
  d1<-subset(slt,celltype2=="5(Oli1)")
  Oli1<-Oli1[Oli1 %in% rownames(d1)]
  for( g in Oli1){
  p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                     comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
  ggsave(paste0("p4/",sid,"_gene_",g,"_Oli1Only_wilcox.pdf"), p, width=5, height=4)
  }
  d1<-subset(slt,celltype2=="13(Mic1)")
  Mic1<-Mic1[Mic1 %in% rownames(d1)]
  for( g in Mic1){
  p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                     comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
  ggsave(paste0("p4/",sid,"_gene_",g,"_Mic1Only_wilcox.pdf"), p, width=5, height=4)
  }
}

library(dplyr)
library(tidyr)
library(ggforce)
library(ggplot2)
library(scatterpie)
dir.create("p5")
aggregate_by_region <- function(df, bin_size = 50, label_col = "Label1") {
  df <- df %>%
    mutate(
      bin_x = floor(centroid_1 / bin_size) * bin_size,
      bin_y = floor(centroid_2 / bin_size) * bin_size,
      bin_id = paste0(bin_x, "_", bin_y)
    )
  
  df_bin <- df %>%
    group_by(bin_id, bin_x, bin_y, !!sym(label_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(bin_id) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()
  
  # 转宽格式，用于 scatterpie
  df_wide <- df_bin %>%
    select(bin_id, bin_x, bin_y, !!sym(label_col), freq) %>%
    pivot_wider(names_from = !!sym(label_col), values_from = freq, values_fill = 0)
  
  return(df_wide)
}

# ✅ 使用示例
for( sid in c("20220128_PFCL1","20220218_PFCL1")){
  slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
  res1 <- readRDS(paste0(sid, "_res1.rds"))
  res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
df_region <- aggregate_by_region(res1, bin_size = 100, label_col = "celltype.stim")
df_region<-df_region[,c("bin_id","bin_x","bin_y",stype)]
cols_for_plot <- setdiff(colnames(df_region), c("bin_id","bin_x","bin_y"))
pal <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                  "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), cols_for_plot)

p <- ggplot() +
  geom_scatterpie(
    aes(x = bin_x, y = bin_y, group = bin_id, r = 40),
    data = df_region,
    cols = cols_for_plot,
    color = NA
  ) +
  coord_equal() +
  theme_scp() +
  scale_fill_manual(values = pal) +
  ggtitle("Spatial scatterpie by region bins")
p

pdf(paste0("p5/",sid,"_scatterpie_region_aggregate.pdf"), width = 14, height = 10)
print(p)
dev.off()
}


for( sid in c("20220128_PFCL1","20220218_PFCL1")){
  slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
  res1 <- readRDS(paste0(sid, "_res1.rds"))
  res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
  for(grp in unique(res1$Group)){
   res_grp<-res1[res1$Group==grp,]
  df_region <- aggregate_by_region(res_grp, bin_size = 100, label_col = "celltype.stim")
  df_region<-df_region[,c("bin_id","bin_x","bin_y",stype)]
  cols_for_plot <- setdiff(colnames(df_region), c("bin_id","bin_x","bin_y"))
  pal <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                    "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), cols_for_plot)
  
  p <- ggplot() +
    geom_scatterpie(
      aes(x = bin_x, y = bin_y, group = bin_id, r = 40),
      data = df_region,
      cols = cols_for_plot,
      color = NA
    ) +
    coord_equal() +
    theme_scp() +
    scale_fill_manual(values = pal) +
    ggtitle("Spatial scatterpie by region bins")
  p
  
  pdf(paste0("p5/",sid,"_",grp,"_scatterpie_region_aggregate.pdf"), width = 14, height = 10)
  print(p)
  dev.off()
  }
}
sid<-"20220128_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 & res1$centroid_2 > -1000,]
res_split[[2]]<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 > -1000,]
res_split[[3]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 &res1$centroid_2 < -1000,]
res_split[[4]]<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 < -1000,]
res_split[[5]]<-res1[res1$centroid_1 > 5000,]
for( i in 1:5){
  message ("processing ",sid," region",i)
  res_tmp<-res_split[[i]]
  slt_tmp<-slt[,rownames(res_tmp)]
  ratio_plot(slt_tmp,sample.by="sampleinf",anno.by ="celltype2",condition.by = "Group",strip.col=cell_colors, 
             facet.ncol =4,save.prefix =paste0("p5/",sid,"_region",i,"_celltype2_ratio"),plot_type = "box",return_data = TRUE)
  message("ratio_plot done")
for(grp in unique(res_tmp$Group)){
  res_grp<-res_tmp[res_tmp$Group==grp,]
  df_region <- aggregate_by_region(res_grp, bin_size = 100, label_col = "celltype.stim")
  df_region<-df_region[,c("bin_id","bin_x","bin_y",stype)]
  cols_for_plot <- setdiff(colnames(df_region), c("bin_id","bin_x","bin_y"))
  pal <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                    "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), cols_for_plot)
  
  p <- ggplot() +
    geom_scatterpie(
      aes(x = bin_x, y = bin_y, group = bin_id, r = 40),
      data = df_region,
      cols = cols_for_plot,
      color = NA
    ) +
    coord_equal() +
    theme_scp() +
    scale_fill_manual(values = pal) +
    ggtitle("Spatial scatterpie by region bins")
  p
  write.xlsx(df_region,file=paste0("p5/",sid,"_",grp,"_region",i,"_df.xlsx"))
  pdf(paste0("p5/",sid,"_",grp,"_region",i,"_scatterpie_region_aggregate.pdf"), width = 7, height = 6)
  print(p)
  dev.off()
}
}

sid<-"20220218_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000,]
res_split[[2]]<-res1[res1$centroid_1 > -2000,]
for( i in 1:2){
  message ("processing ",sid," region",i)
  res_tmp<-res_split[[i]]
  slt_tmp<-slt[,rownames(res_tmp)]
  ratio_plot(slt_tmp,sample.by="sampleinf",anno.by ="celltype2",condition.by = "Group",strip.col=cell_colors, 
             facet.ncol =4,save.prefix =paste0("p5/",sid,"_region",i,"_celltype2_ratio"),plot_type = "box",return_data = TRUE)
  message("ratio_plot done")
  for(grp in unique(res_tmp$Group)){
    res_grp<-res_tmp[res_tmp$Group==grp,]
    df_region <- aggregate_by_region(res_grp, bin_size = 100, label_col = "celltype.stim")
    df_region<-df_region[,c("bin_id","bin_x","bin_y",stype)]
    cols_for_plot <- setdiff(colnames(df_region), c("bin_id","bin_x","bin_y"))
    pal <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                      "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), cols_for_plot)
    
    p <- ggplot() +
      geom_scatterpie(
        aes(x = bin_x, y = bin_y, group = bin_id, r = 40),
        data = df_region,
        cols = cols_for_plot,
        color = NA
      ) +
      coord_equal() +
      theme_scp() +
      scale_fill_manual(values = pal) +
      ggtitle("Spatial scatterpie by region bins")
    p
    if(i==1){
      write.xlsx(df_region,file=paste0("p5/",sid,"_",grp,"_regionLeft_df.xlsx"))
    pdf(paste0("p5/",sid,"_",grp,"_regionLeft_scatterpie_region_aggregate.pdf"), width = 7, height = 6)
    print(p)
    dev.off()
    }
    if(i==2){
      write.xlsx(df_region,file=paste0("p5/",sid,"_",grp,"_regionRight_df.xlsx"))
      pdf(paste0("p5/",sid,"_",grp,"_regionRight_scatterpie_region_aggregate.pdf"), width = 7, height = 6)
      print(p)
      dev.off()
    }
  }
}

library(corrplot)
dir.create("p6")
fhlist<-list.files(path="p5",pattern="df.xlsx")
fhlist<-fhlist[!grepl("celltype2",fhlist)]
for(f in fhlist){
  df_region<-read.xlsx(paste0("p5/",f))
# 假设你的相关矩阵叫 mat
mat <- cor(df_region[,-c(1,2,3)])
write.xlsx(as.data.frame.matrix(mat),paste0("p6/",gsub(".xlsx","_cor.xlsx",f)),rowNames=T)
# 只画下三角，去掉对角线
pdf(paste0("p6/",gsub(".xlsx",".pdf",f)),width=5,height=5)
corrplot(mat, method = "color", type = "lower", diag = FALSE,
         tl.col = "black", tl.srt = 45)
dev.off()
}



dir.create("p4_ht")
key_genes <- c("Tenm2",'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1',
               'Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2',
               'App', 'Lrrc4c', 'Ncam1','Nfasc','Tenm2', 'Sema6a','Ptn','Nrxn3','Nrxn1','Ncam2',
               'Lrrc4c','Nrxn3','Nrxn1','Cadm1','Erbb4','Nrg3','Nlgn1','Cntnap2','Tenm2','Mef2a')
for( sid in c("20220128_PFCL1","20220218_PFCL1")){
  slt<-readRDS(paste0("subtype_",sid,"_spRD.rds"))
  key_genes<-key_genes[key_genes %in% rownames(slt)]
  p1<-GroupHeatmapy(slt,group.by="Group",features=unique(key_genes),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_ht.pdf"), p1$plot, width=12, height=3)
}
int2 <- c("Tenm2",'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1')
Ast1<-c('Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2')	
Oli1<-c('App', 'Lrrc4c', 'Ncam1','Nfasc','Tenm2', 'Sema6a','Ptn','Nrxn3','Nrxn1','Ncam2')	
Mic1<-c('Lrrc4c','Nrxn3','Nrxn1','Cadm1','Erbb4','Nrg3','Nlgn1','Cntnap2','Tenm2','Mef2a')
sid<-"20220218_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000,]
res_split[[2]]<-res1[res1$centroid_1 > -2000,]
for( i in 1:2){
  res_tmp<-res_split[[i]]
  slt_tmp<-slt[,rownames(res_tmp)]
  p1<-GroupHeatmapy(slt_tmp,group.by="Group",features=unique(key_genes),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,".pdf"), p1$plot, width=12, height=3)
  for( g in unique(key_genes)){
    p<-FeatureStatPlot(slt_tmp,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt_tmp,celltype2=="9(Int2)")
  int2<-int2[int2 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(int2),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Int2Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(int2)){
  p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                     comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
  ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Int2Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="6(Ast1)")
  Ast1<-Ast1[Ast1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Ast1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Ast1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Ast1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Ast1Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="5(Oli1)")
  Oli1<-Oli1[Oli1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Oli1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Oli1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Oli1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Oli1Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="13(Mic1)")
  Mic1<-Mic1[Mic1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Mic1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Mic1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Mic1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Mic1Only_wilcox.pdf"), p, width=5, height=4)
  }
}

sid<-"20220128_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 & res1$centroid_2 > -1000,]
res_split[[2]]<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 > -1000,]
res_split[[3]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 &res1$centroid_2 < -1000,]
res_split[[4]]<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 < -1000,]
res_split[[5]]<-res1[res1$centroid_1 > 5000,]
for( i in 1:5){
  res_tmp<-res_split[[i]]
  slt_tmp<-slt[,rownames(res_tmp)]
  p1<-GroupHeatmapy(slt_tmp,group.by="Group",features=unique(key_genes),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,".pdf"), p1$plot, width=12, height=3)
  for( g in unique(key_genes)){
    p<-FeatureStatPlot(slt_tmp,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt_tmp,celltype2=="9(Int2)")
  int2<-int2[int2 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(int2),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Int2Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(int2)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Int2Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="6(Ast1)")
  Ast1<-Ast1[Ast1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Ast1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Ast1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Ast1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Ast1Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="5(Oli1)")
  Oli1<-Oli1[Oli1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Oli1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Oli1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Oli1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Oli1Only_wilcox.pdf"), p, width=5, height=4)
  }
  
  d1<-subset(slt,celltype2=="13(Mic1)")
  Mic1<-Mic1[Mic1 %in% rownames(d1)]
  p1<-GroupHeatmapy(d1,group.by="Group",features=unique(Mic1),flip=T,row_title = "",column_title = "",
                    nlabel=0,show_row_names = T,show_column_names = T,cluster_columns = T,exp_legend_title = "zscore",
                    add_dot=T,add_bg = T)
  p1$plot
  ggsave(paste0("p4_ht/",sid,"_region",i,"_Mic1Only.pdf"), p1$plot, width=5, height=3)
  for( g in unique(Mic1)){
    p<-FeatureStatPlot(d1,stat.by=g,group.by="Group",add_box = T,sig_label = "p.format",
                       comparisons = list(c("CON","DEP"),c("CON","KET"),c("DEP","KET")))
    ggsave(paste0("p4_ht/",sid,"_gene_",g,"_region",i,"_Mic1Only_wilcox.pdf"), p, width=5, height=4)
  }
}


#-------------------------------
sid<-"20220128_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data),]
res_split1<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 & res1$centroid_2 > -1000,]
res_split2<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 > -1000,]
res_split3<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 &res1$centroid_2 < -1000,]
res_split4<-res1[res1$centroid_1 > -2000 & res1$centroid_1 <  5000 & res1$centroid_2 < -1000,]
res_split5<-res1[res1$centroid_1 > 5000,]
res1<-res_split1
res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
res1 <- res1[order(res1$plot_color == ct), ]
cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                          "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
res1<-res1[res1$celltype.stim %in% stype,]
ggplot(res1, aes(x = centroid_1, y = centroid_2, color = celltype.stim)) +
  geom_jitter(size = 2, alpha = 0.6, stroke = 0) +
  scale_color_manual(values = cell_colors) +facet_wrap(~Group)+
  coord_fixed() +
  theme_scp()


df <- res1 %>%
  mutate(bin_x = centroid_1,
         bin_y = centroid_2) %>%
  group_by(bin_x, bin_y, celltype.stim) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(bin_x, bin_y) %>%
  mutate(total = sum(n),
         prop = n / total) %>%
  ungroup() %>%
  mutate(bin_id = paste0(bin_x, "_", bin_y)) %>%
  select(bin_id, bin_x, bin_y, celltype.stim, prop)
df_wide <- df %>%
  pivot_wider(
    names_from = celltype.stim,
    values_from = prop,
    values_fill = 0
  ) %>%
  mutate(across(-c(bin_id, bin_x, bin_y), ~ round(.x, 6))) %>%
  arrange(bin_y, bin_x)  # 可选：按坐标排序

ggplot(df_wide, aes(x = bin_x, y = bin_y)) +
  geom_point(shape=19,color="lightgrey")+
  geom_point(aes(color=`4(Exc4)`),data=df_wide[df_wide$`4(Exc4)` !=0,])+
  #scale_color_gradient2(low,mid=,high=,midpoint = 0.5) +
  scale_color_gradientn(colors=colorRamps::blue2green2red(40))+
  coord_fixed() +
  theme_scp() +
  theme(panel.grid = element_blank()) +
  labs(color = "4(Exc4) proportion", title = "Exc4 spatial distribution")

library(ggplot2)
library(patchwork)
library(colorRamps)
library(ggplot2)
library(patchwork)
library(colorRamps)
library(rlang)  # 用于 !!sym()
celltypes <- stype
plot_list <- lapply(celltypes, function(ct) {
  ggplot(df_wide, aes(x = bin_x, y = bin_y)) +
    geom_point(shape = 19,size=0.8, color = "lightgrey") +  # 底灰点
    geom_point(aes(color = .data[[ct]]),size=0.8, data = df_wide[df_wide[[ct]] != 0, ]) +
    scale_color_gradientn(colors = colorRamps::blue2green2red(40)) +
    coord_fixed() +
    theme_scp() +
    theme(panel.grid = element_blank(), legend.position="top",legend.direction = "vertical",legend.justification = "left",legend.text.align = 0) +
    labs(color = paste0(ct, " proportion"), title = ct)
})
final_plot <- wrap_plots(plot_list, nrow = 2, ncol = 4)
final_plot
pdf("region1_test.pdf",width=14,height=12)
print(final_plot)
dev.off()
#------------------------------------------------------------------


dir.create("p4_128Region1_geneExp")
key_genes <- unique(c("Tenm2",'Tenm4','Flrt3','App','Cadm1','Cntn1','Efna5','Flrt2','Gad1','Slc6a1',
               'Negr1','Nrxn1','Tenm3','Lama2','Lrrc4c','Ncam1','Nrxn2','Nrxn3','Cadm1','Cdh2',
               'App', 'Lrrc4c', 'Ncam1','Nfasc','Tenm2', 'Sema6a','Ptn','Nrxn3','Nrxn1','Ncam2',
               'Lrrc4c','Nrxn3','Nrxn1','Cadm1','Erbb4','Nrg3','Nlgn1','Cntnap2','Tenm2','Mef2a'))
for( sid in c("20220128_PFCL1")){
  slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
  cells_use<-rownames(slt@meta.data)
  cid<-cells_use[cells_use %in% rownames(res_split1)]
  slt<-slt[,cid]
  coords <- res1[cid, c("centroid_1", "centroid_2")]
  coords <- as.matrix(coords)
  colnames(coords) <- c("spatial_1", "spatial_2")
  rownames(coords) <- cid
  # 构建 DimReduc 对象
  spatial_dr <- CreateDimReducObject(
    embeddings = coords,
    key = "SPATIAL_",
    assay = DefaultAssay(slt)
  )
  # 写入 slt@reductions
  slt[["spatial"]] <- spatial_dr
  saveRDS(slt,file=paste0("p4_region1_",sid,"_spRD.rds"))
  key_genes<-key_genes[key_genes %in% rownames(slt)]
    p <- FeatureDimPlot(slt,reduction="spatial",features=sort(key_genes),ncol=7,theme_use = "theme_blank",show_stat = F)
    ggsave(paste0("p4_128Region1_geneExp/28gene_spRD.pdf"), p, width=20, height=12)
}

key_genes<-sort(c("Nrxn2",'Negr1','App','Lrrc4c','Cdh2','Tenm2','Nrxn3','Cntn1','Gad1','Slc6a1',
                  'Ptn','Ncam1','Nrxn1','Nfasc','Cadm1'))
slt<-readRDS("p4_128Region1_geneExp/p4_region1_20220128_PFCL1_spRD.rds")
p<-FeatureDimPlot(slt,reduction="spatial",features=key_genes,ncol=5,theme_use = "theme_blank",show_stat = F) & 
  theme(
    plot.title = element_text(face = "italic"),           # 标题斜体
    strip.text = element_text(face = "italic")            # facet标签斜体
  )
ggsave(paste0("p4_128Region1_geneExp/15gene_spRD.pdf"), p, width=16, height=12)

slt<-subset(slt,celltype2 %in% stype)
Idents(slt)<-slt$celltype2
slt<-RenameIdents(slt,"9(Int2)"=expression("Int2^(Rarb+)"),
                  "13(Mic1)"="Mic1(Runx1+)",
                  "6(Ast1)"="Ast1(Gpc5+)",
                 "5(Oli1)"="Oli1(St18+)",
                  "4(Exc4)"="Exc4(Ndst4+)",
                  "10(Exc7)"="Exc7(Tafa1+)",
                  "17(Exc9)"="Exc9(Tmem163+)",
                 "18(Exc10)"= "Exc10(Abi3bp+)")
slt$type<-slt@active.ident
#cell_colors
cols<-c("#6262ff", "#1F78B4", "#3ba997", "#33A02C", "#FDBF6F", "#FF7F00", "#FB9A99", "#E31A1C")
names(cols)<-c("Exc4(Ndst4+)", "Exc7(Tafa1+)","Exc9(Tmem163+)", "Exc10(Abi3bp+)","Int2(Rarb+)", "Mic1(Runx1+)","Oli1(St18+)","Ast1(Gpc5+)")

p<-CellDimPlot(slt,reduction="spatial",theme_use="theme_blank",split.by="type",group.by="type",bg_color = "lightgrey",ncol=4,legend.position = "none",
            palcolor = list(c("#FDBF6F"    ,  "#FF7F00"  ,    "#E31A1C"    ,  "#FB9A99",      "#6262ff"   ,   "#1F78B4" ,     "#3ba997"    ,  "#33A02C")))
p
pdf("p4_128Region1_geneExp/p1.pdf",width=14,height=8)
print(p)
dev.off()

library(ggplot2)
library(cowplot)

library(ggplot2)
library(cowplot)

setwd("/data2/Project/anding/test_20250923/")
sid<-"20220128_PFCL1"
slt <- readRDS(paste0("subtype_",sid, "_check_meta.rds"))
slt$celltype2<-factor(as.character(slt$celltype2),levels=stype)
res1 <- readRDS(paste0(sid, "_res1.rds"))
res1<-res1[rownames(res1) %in% rownames(slt@meta.data) & res1$celltype.stim %in% stype,]
res_split<-list()
res_split[[1]]<-res1[res1$centroid_1 < -2000 & res1$centroid_1 < 5000 & res1$centroid_2 > -1000,]
res<-res_split[[1]]

stype<-c("4(Exc4)","10(Exc7)","17(Exc9)","18(Exc10)","9(Int2)","13(Mic1)","5(Oli1)","6(Ast1)")
cell_colors <- setNames(c("#6262ff", "#1F78B4", "#3ba997", "#33A02C" ,
                          "#FDBF6F", "#FF7F00", "#FB9A99" ,"#E31A1C"), stype)
plist<-list()

for (ct in stype) {
  res1 <- res
  res1$plot_color <- ifelse(res1$celltype == ct, ct, "Other")
  res1 <- res1[order(res1$plot_color == ct), ]
  # 根据 ct 设置 Label2
  if (ct %in% stype[1:4]) {
    labels <- c("Exc L2/3 IT", "Exc L4/5 IT", "Exc L5 IT", "Exc L5/6 NP")
    res1 <- res1[res1$Label2 %in% labels, ]
  } else if (ct %in% stype[5]) {
    labels <- c("Int Lamp5", "Int Pvalb", "Int Sst", "Int GABAergic")
    res1 <- res1[res1$Label2 %in% labels, ]
  } else {
    labels <- c("Non-Neuron Astrocytes", "Non-Neuron Endo","Non-Neuron Oligo",
                "Non-Neuron Microglia","Non-Neuron OPC")
    res1 <- res1[res1$Label2 %in% labels, ]
  }
  color_values <- c("Other" = "lightgrey", ct = cell_colors[ct])
  names(color_values)[2] <- ct
  # 每个 Label2 生成单独图
  facet_list <- lapply(labels, function(lbl) {
    df <- res1[res1$Label2 == lbl, ]
    ggplot(df, aes(x = centroid_1, y = centroid_2, color = plot_color)) +
      geom_point(size = 0.2, alpha = 0.8) +
      scale_color_manual(values = color_values) +
      labs(fill = "type", x = "",y="") +
      theme_void() + scale_x_continuous(breaks = c(-7000, -5000))+
      coord_fixed() +
      guides(color = "none") +
      ggtitle(lbl)+theme(
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
  })
  
  plist[[ct]] <- facet_list
}
theme_y<-theme(axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
empty_plot <- ggplot() + theme_void()

legend_plot <- ggplot(data.frame(stype, x=1, y=1), aes(x, y, color=stype)) +
  geom_point(size=4) +
  scale_color_manual(values = cell_colors, name = "Type") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size=12, face="bold"),
    legend.text = element_text(size=10),
    legend.key.size = unit(0.6, "cm")
  )

legend_grob <- cowplot::get_legend(legend_plot)
theme_small_margin <- theme(  plot.margin = margin(-2, -5, -2,-5)) # 上右下左，单位 pt
theme_small_margin1 <- theme(  plot.margin = margin(-2, -5, -2,0)) # 上右下左，单位 pt

combine_plot<-plot_grid(plist$`9(Int2)`[[1]]+labs(y="Int2")+theme_y+theme_small_margin1,
                        plist$`9(Int2)`[[2]]+theme_small_margin,
                        plist$`9(Int2)`[[3]]+theme_small_margin,
                        plist$`9(Int2)`[[4]]+theme_small_margin,empty_plot,
                        plist$`13(Mic1)`[[1]]+labs(y="Mic1")+theme_y+theme_small_margin1,
                        plist$`13(Mic1)`[[2]]+theme_small_margin,
                        plist$`13(Mic1)`[[3]]+theme_small_margin,
                        plist$`13(Mic1)`[[4]]+theme_small_margin,
                        plist$`13(Mic1)`[[5]]+theme_small_margin,
                        plist$`6(Ast1)`[[1]]+labs(y="Ast1")+theme_y+theme_small_margin1,
                        plist$`6(Ast1)`[[2]]+theme_small_margin,
                        plist$`6(Ast1)`[[3]]+theme_small_margin,
                        plist$`6(Ast1)`[[4]]+theme_small_margin,
                        plist$`6(Ast1)`[[5]]+theme_small_margin,
                        plist$`5(Oli1)`[[1]]+labs(y="Oli1")+theme_y+theme_small_margin1,
                        plist$`5(Oli1)`[[2]]+theme_small_margin,
                        plist$`5(Oli1)`[[3]]+theme_small_margin,
                        plist$`5(Oli1)`[[4]]+theme_small_margin,
                        plist$`5(Oli1)`[[5]]+theme_small_margin,
                        plist$`4(Exc4)`[[1]]+labs(y="Exc4")+theme_y+theme_small_margin1,
                        plist$`4(Exc4)`[[3]]+theme_small_margin,
                        plist$`4(Exc4)`[[4]]+theme_small_margin,
                        plist$`4(Exc4)`[[2]]+theme_small_margin,empty_plot,
                        plist$`10(Exc7)`[[1]]+labs(y="Exc7")+theme_y+theme_small_margin1,
                        plist$`10(Exc7)`[[3]]+theme_small_margin,
                        plist$`10(Exc7)`[[4]]+theme_small_margin,
                        plist$`10(Exc7)`[[2]]+theme_small_margin,empty_plot,
                        plist$`17(Exc9)`[[1]]+labs(y="Exc9")+theme_y+theme_small_margin1,
                        plist$`17(Exc9)`[[3]]+theme_small_margin,
                        plist$`17(Exc9)`[[4]]+theme_small_margin,
                        plist$`17(Exc9)`[[2]]+theme_small_margin,empty_plot,
                        plist$`18(Exc10)`[[1]]+labs(y="Exc10")+theme_y+theme_small_margin1,
                        plist$`18(Exc10)`[[3]]+theme_small_margin,
                        plist$`18(Exc10)`[[4]]+theme_small_margin,
                        plist$`18(Exc10)`[[2]]+theme_small_margin,empty_plot,ncol=5,
                        align = "vh",axis="l",  rel_heights = rep(1, 8*5),
                        rel_widths  = rep(1, 5),   
                        hjust = -2,                
                        vjust = -2 )
final_plot <- cowplot::plot_grid(
  combine_plot, legend_grob,
  ncol = 2,
  rel_widths = c(5,1)
)

pdf("merge1.pdf",width=12,height=16)
print(final_plot)
dev.off()