library(Seurat)
library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(ggrastr)
library(magick)
library(stringr)
library(ComplexHeatmap)
library(ggridges)
library(ggpubr)
library(ggbeeswarm)
library(AUCell)

###########################################################
###function
#parse geneset from MSigDB
parse_geneset_term<-function(geneSet="~/work/MSigDB/human/v6.2/h.all.v6.2.symbols.gmt"){
	tmp<-readLines(geneSet)
	gs<-list()
	for(i in 1:length(tmp)){
		tmp2<-strsplit(tmp[[i]],split="\t")[[1]]
		gs[[tmp2[1]]]<-tmp2[-1:-2]
	}
	gs
}

#get coordination information of image
get_spot_coord<-function(json_file="alignment.json"){
	require(rjson)
	fromJSON(file=json_file)->alignment
	coord<-do.call(rbind,lapply(alignment$oligo,function(x){
		if(length(x)==7) x$tissue<-FALSE
		data.frame(x[c("x","y","row","col","dia","imageX","imageY","tissue")])
	}))
	rownames(coord)<-paste(coord$row,coord$col,sep="_")
	fiducial<-do.call(rbind,lapply(alignment$fiducial,function(x){
		if(length(x)==7){x$fidName<-""}
		data.frame(x[c("x","y","row","col","dia","imageX","imageY","fidName")])
	}))
	rownames(fiducial)<-paste(fiducial$row,fiducial$col,sep="_")
	list(coord=coord,fiducial=fiducial)
}

#parse image and save a low resolution copy
get_image<-function(image_path,resolution=600){
	require(magick)
	low_res<-paste0(image_path,".lowres.png")
	image_read(image_path)->img
	orig.info<-image_info(img)
	if(!file.exists(low_res)){
		image_write(image_scale(img,resolution),low_res,format="png",depth=16)
	}
	reduced.info<-image_info(image_read(low_res))
	scalef<-min(reduced.info$height/orig.info$width,reduced.info$width/orig.info$width)
	list(low_res=low_res,scalef=scalef)
}

#add image path and coordination information into seurat object
#saved in misc slot
add_image_seurat<-function(so,img_path,coord_path){
	coord_fid<-get_spot_coord(coord_path)
	img<-get_image(img_path)
	
	so$RowCol<-paste(so$col-1,so$row-1,sep="_")
	coord<-coord_fid$coord
	rownames(coord)<-paste(coord$row,coord$col,sep="_")
	so$imageX<-coord[so$RowCol,"imageX"]
	so$imageY<-coord[so$RowCol,"imageY"]
	so$tissue<-coord[so$RowCol,"tissue"]
	
	Misc(so,slot="img")<-img$low_res
	Misc(so,slot="scalef")<-img$scalef
	Misc(so,slot="fiducial")<-coord_fid$fiducial
	so
}

#create seurat object from a path with cellranger output structure
createSO<-function(path="MTF"){
	so = CreateSeuratObject(counts = Read10X(paste0(path,"/counts_unfiltered/cellranger/")), assay = "Spatial", project=path)
	#genes
	read.table(paste0(path,"/counts_unfiltered/cellranger/genes.tsv"),header=F,sep="\t")->symbol
	unique(symbol[str_detect(symbol$V1,"^ENSG"),"V2"])->hg38
	setdiff(symbol$V2,hg38)->mm10
	so$nCount_hg38<-Matrix::colSums(so[hg38,]@assays$Spatial@counts) 
	so$nCount_mm10<-Matrix::colSums(so[mm10,]@assays$Spatial@counts)
	so$perc_nCount<-so$nCount_mm10/(so$nCount_Spatial)
	so$nFeature_hg38<-Matrix::colSums(so[hg38,]@assays$Spatial@counts>=1)
	so$nFeature_mm10<-Matrix::colSums(so[mm10,]@assays$Spatial@counts>=1)
	so$perc_nFeature<-so$nFeature_mm10/(so$nFeature_Spatial)
	#mitochondria and ribosome RNA
	so$perc_mt <- PercentageFeatureSet(so, pattern = "^(mt|MT)-")
	so$perc_rs <- PercentageFeatureSet(so, pattern = "^R(p|P)(s|l|S|L)\\d+")
	so$perc_hb <- PercentageFeatureSet(so, pattern = "^Hb.*-")
	#
	read.table("/opt/tools/spaceranger-1.3.1/lib/python/cellranger/barcodes/visium-v1_coordinates.txt",row.names=1)->coord
	rownames(coord)<-paste(rownames(coord),1,sep="-")
	colnames(coord)<-c("row","col")
	so$row<-coord[colnames(so),"row"]
	so$col<-coord[colnames(so),"col"]
	so
}
#proto
GeomSpatial <- ggproto("GeomSpatial",Geom,
	setup_data = function(self, data, params) {
	  data <- ggproto_parent(Geom, self)$setup_data(data, params)
	  data
	},
	draw_group = function(data, panel_scales, coord) {
		ggplot2:::ggname("spatial",grid::rasterGrob(data$img[[1]],x=data$x,y=data$y,
			default.units = "npc"
		))
	},
	required_aes = c("img","x","y")
)
#geom class for spatial image
geom_spatial <-  function(mapping = NULL,
	data = NULL,
	stat = "identity",
	position = "identity",
	na.rm = FALSE,
	show.legend = NA,
	inherit.aes = FALSE,
	...) {
	
	layer(
		geom = GeomSpatial, 
		mapping = mapping,
		data = data,
		stat = stat,
		position = position,
		show.legend = show.legend,
		inherit.aes = inherit.aes,
		params = list(na.rm = na.rm, ...)
	)
}
#parse low resolution image to get height and width
parse_image<-function(path=image_paths){
	img<-readbitmap::read.bitmap(path)
	img.d<-dim(img)
	tibble(img=list(img),height=img.d[1],width=img.d[2])
}

#plot feature in seurat object, with image and fiducial at bottom layer
plot_feature<-function(object,feature="nCount_Spatial",slot="data",
	image_only=FALSE,add_image=FALSE,add_fiducial=FALSE,
	alpha=c(0.05, 0.95),limits=c(0.1,0.9),pt.size=2){
	require(ggpubr)
	if(image_only){
		add_image<-TRUE
		add_fiducial<-TRUE
	}
	
	data.df<-object[[]]
	if(! feature %in% colnames(data.df)){
		fdata<-FetchData(object,feature,slot=slot)
		colnames(fdata)<-feature
		data.df<-cbind(data.df,fdata)
	}
	
	if(add_image) {
		scalef<-Misc(object,slot="scalef")
		img<-parse_image(Misc(object,"img"))
		g<-ggplot(data.df,aes(x=imageX*scalef,y=imageY*scalef))+
			geom_spatial(data=img, aes(img=img), x=0.5, y=0.5)
	}
	else{
		g<-ggplot(data.df,aes(x=row,y=col))
	}
	if(add_fiducial) {
		fid<-as.data.frame(Misc(object,"fiducial"))
		g<-g+geom_point(data=fid,color="black",shape=21,stroke=0,show.legend=F,size=pt.size)
	}
	
	if(!image_only){
		g<-g+geom_point(aes_(fill=as.name(feature),alpha=as.name(feature)),
			shape=21,stroke=NA,size=pt.size)
	}
	else{g<-g+ggtitle("Image only")}
	g<-g+ggpubr::theme_pubr(border=T)
	if(add_image) {
		g<-g+xlim(0,img$width)+ylim(img$height,0)+
			coord_fixed(expand=F)+Seurat::NoAxes()
	}
	else{
		g<-g+coord_cartesian(ylim=c(77,0),xlim=c(0,127))
	}
	g<-g+labs(x=NULL,y=NULL)+
		scale_fill_gradientn(limits = quantile(data.df[,feature],limits,na.rm=T),
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)+
		scale_alpha(range = alpha,guide=NULL)
	list(ggplot=g,data=data.df)
}
#QC plot
plot_QC<-function(object,outfile="QC.pdf"){
	pdf(outfile)
	print(plot_feature(object,feature="nCount_Spatial")$ggplot)
	print(plot_feature(object,feature="perc_nCount",alpha=c(0.1,1))$ggplot)
	print(plot_feature(object,feature="nFeature_Spatial")$ggplot)
	print(plot_feature(object,feature="perc_nFeature",alpha=c(0.1,1))$ggplot)
	print(plot_feature(object,feature="perc_mt",alpha=c(0.1,1),limits=c(0.05, 0.95))$ggplot)
	print(plot_feature(object,feature="perc_rs",alpha=c(0.1,1),limits=c(0.05, 0.95))$ggplot)
	print(plot_feature(object,feature="perc_hb",alpha=c(0.1,1),limits=c(0.05, 0.95))$ggplot)
	dev.off()
}
##end of function
################################################

########################
#Main#
########################

#create seurat object
MTF<-createSO("MTF")
#MTF<-add_image_seurat(MTF,img_path="image/Control_Composite (RGB)_frame.jpeg",coord_path="image/MTF.alignment.json")
PAR<-createSO("Parental")
plot_QC(MTF,"QC.MTF.pdf")
plot_QC(PAR,"QC.Parental.pdf")

pdf("Elbowplot.pdf",width=4,height=3)
ElbowPlot(MTF,ndims=50)
ElbowPlot(PAR,ndims=50)
dev.off()

###pre-process 
MTF <- MTF %>% 
	SCTransform(assay = "Spatial", verbose = FALSE) %>%
	RunPCA(assay = "SCT", verbose = FALSE) %>%
	FindNeighbors(reduction = "pca", dims = 1:30) %>%
	FindClusters(verbose = FALSE) %>%
	RunUMAP(reduction = "pca", dims = 1:30) %>%
	RunTSNE(reduction = "pca", dims = 1:30)
PAR <- PAR %>% 
	SCTransform(assay = "Spatial", verbose = FALSE) %>%
	RunPCA(assay = "SCT", verbose = FALSE) %>%
	FindNeighbors(reduction = "pca", dims = 1:30) %>%
	FindClusters(verbose = FALSE) %>%
	RunUMAP(reduction = "pca", dims = 1:30) %>%
	RunTSNE(reduction = "pca", dims = 1:30)
	

###
pdf("Dimreduction.MTF.pdf")
DimPlot(MTF, reduction = "tsne", label = TRUE)
DimPlot(MTF, reduction = "umap", label = TRUE)
MTF[[]]->test2
label<-as.data.frame(t(sapply(split(test2[,c("row","col")],test2$seurat_clusters),function(x){apply(x,2,median)})))
label$label<-rownames(label)
ggplot(subset(test2,nCount_Spatial>=200),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.25,size=2)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(78,0))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete()+
	guides(fill=guide_legend(nrow=2))
color<-scales:::hue_pal()(nlevels(test2$seurat_clusters))
lapply(levels(test2$seurat_clusters),function(x){
	test2$seurat_clusters<-test2$seurat_clusters==x
	g<-ggplot(subset(test2,nCount_Spatial>=200),aes(x=row,y=col))+
		geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.25,size=2)+
		ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(78,0))+labs(x=NULL,y=NULL,alpha=NULL,title=paste0("cluster",x))+
		scale_fill_manual(values=c("white",color[as.numeric(x)+1]),guide=NULL)
	print(g)
	NULL
})
dev.off()

#####
pdf("Dimreduction.Parental.pdf")
DimPlot(PAR, reduction = "tsne", label = TRUE)
DimPlot(PAR, reduction = "umap", label = TRUE)

PAR[[]]->test3
label<-as.data.frame(t(sapply(split(test3[,c("row","col")],test3$seurat_clusters),function(x){apply(x,2,median)})))
label$label<-rownames(label)
ggplot(subset(test3,nCount_Spatial>=200),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.25,size=2)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(78,0))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete()+
	guides(fill=guide_legend(nrow=2))
color<-scales:::hue_pal()(nlevels(test3$seurat_clusters))
lapply(levels(test3$seurat_clusters),function(x){
	test3$seurat_clusters<-test3$seurat_clusters==x
	g<-ggplot(subset(test3,nCount_Spatial>=200),aes(x=row,y=col))+
		geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.25,size=2)+
		ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(78,0))+labs(x=NULL,y=NULL,alpha=NULL,title=paste0("cluster",x))+
		scale_fill_manual(values=c("white",color[as.numeric(x)+1]),guide=NULL)
	print(g)
	NULL
})
dev.off()

#####
#keep only spot with high quality
MTF.ss <- subset(MTF, idents = c(3,8))
PAR.ss <- subset(PAR, idents = c(3,5:8,10))
merged <- merge(MTF.ss, PAR.ss)

DefaultAssay(merged) <- "SCT"
VariableFeatures(merged) <- c(VariableFeatures(MTF.ss), VariableFeatures(PAR.ss))
merged <- RunPCA(merged,verbose = FALSE) 
pdf("Elbowplot.pdf",width=4,height=3)
ElbowPlot(merged,ndims=50)
dev.off()

merged <- merged %>% 
	FindNeighbors(reduction = "pca",dims = 1:30) %>%
	FindClusters(verbose = FALSE) %>%
	RunUMAP(reduction = "pca", dims = 1:30) %>%
	RunTSNE(reduction = "pca", dims = 1:30)
merged$Sample<-merged$orig.ident

#cluster plot
pdf("Dimreduction.merged.pdf",width=3.5,height=3)
DimPlot(merged, reduction = "tsne", label = TRUE)
DimPlot(merged, reduction = "tsne", label = TRUE, group.by="Sample")
DimPlot(merged, reduction = "tsne", label = TRUE, split.by="Sample")

merged[[]]->test4
label<-as.data.frame(t(sapply(split(test4[,c("row","col")],test4$seurat_clusters),function(x){apply(x,2,median)})))
label$label<-rownames(label)
ggplot(subset(test4,Sample=="MTF"),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F,size=2)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete(guide=NULL)#+guides(fill=guide_legend(nrow=2))
ggplot(subset(test4,Sample=="Parental"),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F,size=2)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete(guide=NULL)#+guides(fill=guide_legend(nrow=2))
dev.off()

###
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmp <- merged.markers %>% group_by(cluster)
write.table(tmp,file="marker genes.txt",sep="\t",quote=F)

#####
require(ComplexHeatmap)
sort(merged$seurat_clusters)->clusters
rlv<-c(1,5,6,8,0,2,9,4,7,3)
factor(clusters,levels=rlv)->clusters
color<-scales:::hue_pal()(nlevels(clusters))

merged.markers<-merged.markers[order(factor(merged.markers$cluster,levels=rlv)),]

merged@assays$SCT@scale.data[,]->dd2
read.table("MTF/counts_unfiltered/cellranger/genes.tsv",header=F,sep="\t")->symbol
unique(symbol[str_detect(symbol$V1,"^ENSG"),"V2"])->hg38
setdiff(symbol$V2,hg38)->mm10

hm<-intersect(intersect(merged.markers$gene,rownames(dd2)),hg38)
mm<-intersect(intersect(merged.markers$gene,rownames(dd2)),mm10)
pdf("Heatmap.marker.pdf",width=5,height=6)
Heatmap(dd2[c(hm,mm),names(clusters)],
	show_column_names=F,cluster_columns=F,cluster_rows=F,name="Expression",
	col=circlize::colorRamp2(seq(-2,2,length=50),colors=Seurat:::PurpleAndYellow(50)),
	show_row_names=T,row_names_gp=gpar(fontsize=1),row_names_side ="left",
	column_split = clusters,column_gap = unit(0.5, "mm"),#column_title = NULL,
	row_split= rep(c("hg38","mm10"),c(length(hm),length(mm))),
	column_title_gp = gpar(fontsize=12),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()

#####
merged$ratio<-colSums(GetAssayData(merged[rownames(merged) %in% mm10,],assay="SCT",slot="data"))/colSums(GetAssayData(merged,assay="SCT",slot="data"))
merged$hg38<- 1-merged$ratio
pdf("Mouse_ratio.pdf")
plot_feature(merged[,merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="MTF")
plot_feature(merged[,merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Parental")
plot_feature(merged[,merged$seurat_clusters==0 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster0")
plot_feature(merged[,merged$seurat_clusters==1 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster1")
plot_feature(merged[,merged$seurat_clusters==2 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster2")
plot_feature(merged[,merged$seurat_clusters==3 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster3")
plot_feature(merged[,merged$seurat_clusters==4 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster4")
plot_feature(merged[,merged$seurat_clusters==5 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster5")
plot_feature(merged[,merged$seurat_clusters==6 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster6")
plot_feature(merged[,merged$seurat_clusters==7 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster7")
plot_feature(merged[,merged$seurat_clusters==8 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster8")
plot_feature(merged[,merged$seurat_clusters==9 & merged$Sample=="Parental"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster9")
dev.off()

####
#plot expression in both samples side by side
plot_feature_gene<-function(object=merged,gene="KLK3",limits=NULL){
	so1<-object[,object$Sample=="MTF" & object$seurat_clusters %in% c(4,7)]
	so2<-object[,object$Sample=="Parental" & object$seurat_clusters %in% c(2,5,8)]
	rg<-range(c(FetchData(so1,gene)[,1],FetchData(so2,gene)[,1]))
	limits<-limits %||% round(rg,2)
	p1<-plot_feature(so1,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="MTF")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	p2<-plot_feature(so2,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="Parental")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	gridExtra::grid.arrange(p1,p2,ncol=2)
}

genes<-c("KLK3", "KLK2", "FOLH1", "TMPRSS2", "NKX3-1", "FOXA1", "DDX17", "DDX5", "NCOR1", "SPP1", "PARK7",
	"FGF23", "FOXP1", "AP1M1", "CADM1", "CNNM2",
	"CTNND1", "CTNNB1",
	"HLA-B","HLA-DRA","CYBA","HLA-A","TUBB","CALR","CORO1A","ACTG1","HLA-E",
	"CTSD","CTSB", "CTSC", "CD63",
	"VIM", "RHOA", "CDC42", "MSN", "RHOB", 
	"CXCR4",
	"MYC",
	"RPS17","RPL32","RPS19","RPL22","RPL11","RPL14","RPS11","RPL29","RPL18","RPS27A","RPL28"
)
pdf("Map.genes.pdf",height=4,width=7)
lapply(genes,function(x){
	plot_feature_gene(merged,x)
	NULL
})
dev.off()

##################
#3.2, 3.2, 1.1, 1.2, 8.1, 8.2
merged<-StashIdent(merged,save.name="old")
as.character(merged$seurat_clusters)->clusters
clusters[which(merged$seurat_clusters==3&merged$Sample=="MTF")]<-"3.2"
clusters[which(merged$seurat_clusters==3&merged$Sample=="Parental")]<-"3.1"
clusters[which(merged$seurat_clusters==1&merged$Sample=="MTF")]<-"1.2"
clusters[which(merged$seurat_clusters==1&merged$Sample=="Parental")]<-"1.1"
clusters[which(merged$seurat_clusters==8&merged$Sample=="MTF")]<-"8.2"
clusters[which(merged$seurat_clusters==8&merged$Sample=="Parental")]<-"8.1"

rlv<-c(6,8.1,1.1,5,3.1,9,2,0,1.2,8.2,3.2,4,7)
factor(clusters,levels=rlv)->clusters
Idents(merged)<-clusters
merged$seurat_clusters<-clusters
names(clusters)<-colnames(merged)


######
#pie plot
p1<-table(subset(merged[[]],Sample=="MTF")$old)
p2<-table(subset(merged[[]],Sample=="Parental")$old)

p1<-data.frame(round(p1/sum(p1)*100,2))
p2<-data.frame(round(p2/sum(p2)*100,2))
p1$label<-paste0(p1$Var1,":",p1$Freq,"%")
p2$label<-paste0(p2$Var1,":",p2$Freq,"%")

color<-scales:::hue_pal()(10)
color[4]<-"#A020F0"

pdf("Pie.cluster.pdf",height=3,width=3)
ggpie(p1, "Freq", label = "label",
   lab.pos = "out",
   fill = "Var1", color = NA,
   palette = color)+guides(fill=guide_legend(title=NULL))+
   coord_polar(theta = "y", start = 0,clip="off")+
   theme(axis.text.x=element_blank())
ggpie(p2, "Freq", label = "label",
   lab.pos = "out", 
   fill = "Var1", color = NA,
   palette = color)+guides(fill=guide_legend(title=NULL))+
   coord_polar(theta = "y", start = 0,clip="off")+
   theme(axis.text.x=element_blank())
dev.off()

######
#ratio plot
merged$ratio<-colSums(GetAssayData(merged[rownames(merged) %in% mm10,],assay="SCT",slot="data"))/colSums(GetAssayData(merged,assay="SCT",slot="data"))
ratio_df<-merged[[]][,c("perc_nCount","ratio","seurat_clusters")]
ratio_df$seurat_clusters<-factor(ratio_df$seurat_clusters,rlv)
ratio_df$hg38<- 1-ratio_df$ratio
ratio_df2<-reshape2::melt(ratio_df[,c(3,2,4)])
#quasirandom, pseudorandom, smiley or frowney

g<-ggplot(ratio_df2,aes(seurat_clusters,value*100,color=variable))+
	#geom_violin()+
	geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
	geom_beeswarm(aes(group=variable),priority="random",method="frowney",
		cex=0.3,dodge.width=0.5,size=0.1,show.legend=T,position=position_dodge2()
	)
tmp<-ggplot_build(g)
yrange<-tmp$layout$panel_scales_y[[1]]$range$range
yrange[2]<-yrange[2]+5
yrange[1]<-yrange[1]-2
xmin=seq(0.5,by=1,length=13)
xmax=seq(1.5,by=1,length=13)
fill=rep_len(c("white","blue"),13)
alpha=rep_len(c(0,0.05),13)

pdf("Ratio.beeswarm.pdf",width=4.8,height=2)
ggplot(ratio_df2,aes(seurat_clusters,value*100,color=variable))+
	annotate("rect",xmin=xmin,xmax=xmax,ymin=yrange[1],ymax=yrange[2],fill=fill, alpha=alpha)+
	scale_x_discrete()+
	geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
	geom_beeswarm(aes(group=variable),priority="random",method="frowney",
		cex=0.3,dodge.width=0.5,size=0.1,show.legend=T,position=position_dodge2())+
	theme_pubr(border=T)+theme(legend.position="right")+
	scale_color_manual(values=c("blue","red"),labels=c("mm10","hg38"))+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="Percentage (%)",color=NULL) + 
	coord_cartesian(expand=F)
dev.off()
########

###############################
###############################
#single cell level pathway enrichment
###################
custom1<-parse_geneset_term("custom.gs.human.txt") #human
custom2<-parse_geneset_term("custom.gs.mouse.txt") #mouse
names(custom2)<-paste0("M-",names(custom2))
parse_geneset_term("h.all.v6.1.symbols.gmt")->Hallmark
parse_geneset_term("h.all.v6.1.symbols_mouse.gmt")->hallmark.mm10
names(hallmark.mm10)<-paste0("M-",names(Hallmark))

merged[rownames(merged) %in% hg38,]->merged.hg38
AUCell_buildRankings(merged.hg38@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
AUCell_calcAUC(c(Hallmark,custom1),cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
merged.hg38[["Hallmark"]]<-CreateAssayObject(data=cells_AUC@assays@data$AUC)

merged[rownames(merged) %in% mm10,]->merged.mm10
AUCell_buildRankings(merged.mm10@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings2
AUCell_calcAUC(c(hallmark.mm10,custom2),cells_rankings2,aucMaxRank=500,nCores=32)->cells_AUC2 
merged.mm10[["Hallmark"]]<-CreateAssayObject(data=cells_AUC2@assays@data$AUC)

#
pw<-c("NK cell","MDSCs","M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)
#
m1<-do.call(cbind,lapply(split(ratio_df[,-27],ratio_df$seurat_clusters),function(x)apply(x,2,mean)))
m2<-do.call(cbind,lapply(split(ratio_df2[,-27],ratio_df2$seurat_clusters),function(x)apply(x,2,mean)))
rownames(m1)<-rownames(m2)<-pw
#lvs<-names(sort(apply(m1[,8:10],1,sum)))
lvs<-rev(pw[c(1,2,4,3,5,7,9,8,25,10,16,17,18,15,20,12,11,14,6,21,19,23,24,26,22,13)])
g1<-ggplot(reshape2::melt(m1))+geom_point(aes(x=factor(Var2,levels=rlv),y=factor(Var1,levels=lvs),size=value),shape=19,color="red")+ #A50026
	labs(x=NULL,y=NULL)+
	scale_x_discrete(position="top")+
	scale_size_area(max_size=3,breaks=c(0.05,0.1,0.15),limits=c(0,0.15),oob=scales::squish)+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.x= element_text(size=8),axis.text.y= element_text(size=8),axis.ticks=element_blank())
g2<-ggplot(reshape2::melt(m2))+geom_point(aes(x=factor(Var2,levels=rlv),y=factor(Var1,levels=lvs),size=value),shape=19,color="blue")+ #313695
	labs(x=NULL,y=NULL)+
	scale_y_discrete(labels=NULL,breaks=NULL)+
	scale_x_discrete(position="top")+
	scale_size_area(max_size=3,breaks=c(0.04,0.08,0.12),limits=c(0,0.12),oob=scales::squish)+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.x= element_text(size=8),axis.text.y= element_text(size=8),axis.ticks=element_blank())
pdf("Pathway.circle.pdf",height=4,width=7)
g1 %+% g2
dev.off()
#t.test(subset(ratio_df2,seurat_clusters %in% c(6,8.1,1.1,5,3.1,9,2,0))[,"hallmark_B cell"],
#	subset(ratio_df2,!seurat_clusters %in% c(6,8.1,1.1,5,3.1,9,2,0))[,"hallmark_B cell"])
#
plot_feature_gene_v2<-function(object=merged,gene="KLK3",limits=NULL){
	so.list1<-object[,object$Sample=="MTF" & object$seurat_clusters %in% c(4,7)]
	so.list2<-object[,object$Sample=="Parental" & object$seurat_clusters %in% c(1.1,5,6)]
	rg<-range(c(FetchData(so.list1,gene)[,1],FetchData(so.list2,gene)[,1]))
	limits<-limits %||% round(rg,2)
	p1<-plot_feature(so.list1,gene,alpha=c(1,1),pt.size=1)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="MTF")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	p2<-plot_feature(so.list2,gene,alpha=c(1,1),pt.size=1)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="Parental")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	gridExtra::grid.arrange(p1,p2,ncol=2)
}

#map plot for individual gene or pathway
pw1<-c("M1 marker","M2 marker","General Macrophage","MDSCs","NK cell",
	"VIM","MYC","TUBB"
)
pw2<-c("M1 marker","M2 marker","General Macrophage","MDSCs","NK cell",
	"Vim","Col1a1","Ctsb"
)
pdf("Map.tmp.human.all.pdf",height=4,width=7)
lapply(pw1,function(x){
	plot_feature_gene_v2(merged.hg38,x)
	NULL
})
dev.off()

pdf("Map.tmp.mouse.all.pdf",height=4,width=7)
lapply(pw2,function(x){
	plot_feature_gene_v2(merged.mm10,x)
	NULL
})
dev.off()

#####
as.data.frame(t(merged.hg38@assays$Hallmark@data[,names(clusters)]))->tmp1
tmp1$Cluster<-clusters
as.data.frame(t(merged.mm10@assays$Hallmark@data[,names(clusters)]))->tmp2
tmp2$Cluster<-clusters

pdf("M1_M2.plot.pdf",width=3.8,height=2.3)
lapply(c("M1 marker","M2 marker","Myc Targets V1","Androgen Response"),function(x){
	g<-ggplot(tmp1,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="AUC score",color=NULL,title=paste0(x," in Human")) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
lapply(c("M1 marker","M2 marker"),function(x){
	g<-ggplot(tmp2,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="area",size=0.2,show.legend=F,color="blue")+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="AUC score",color=NULL,title=paste0(x," in Mouse")) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
lapply(c("NK cell","MDSCs"),function(x){
	g<-ggplot(tmp2,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="area",size=0.2,show.legend=F,color="blue")+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="AUC score",color=NULL,title=paste0(x," in Mouse")) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()
