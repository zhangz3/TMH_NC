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

#Load10X_Spatial("test/outs/","filtered_feature_bc_matrix.h5")->tmp

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

plot_spatial<-function(data=bcs_merge,
	fill="sum_umi",size=0.7,stroke=0,rast=TRUE,
	add_image=T,img_data=NULL,width=img_data$width,height=img_data$height,
	legend.title="Total UMI"){
	g<-ggplot(data=data,aes_string(x="imagecol",y="imagerow",fill=fill))
	if(add_image){
		g<-g+geom_spatial(data=img_data, aes(img=img), x=0.5, y=0.5)
	}
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
	if(rast){
		require(ggrastr)
		g<-g+geom_point_rast(shape = 21, colour = "black", size = size, stroke = stroke)
	}
	else{
		g<-g+geom_point(shape = 21, colour = "black", size = size, stroke = stroke)
	}
	g<-g+coord_fixed(expand=FALSE)+
		scale_fill_gradientn(colours = myPalette(100))+
		xlim(0,width)+
		ylim(height,0)+
		xlab("") +
		ylab("") +
		labs(fill = legend.title)+
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(), 
			axis.line = element_blank(),
			axis.text = element_blank(),
			axis.ticks = element_blank()
		)
	g
}

parse_image<-function(path=image_paths){
	img<-readbitmap::read.bitmap(path)
	img.d<-dim(img)
	tibble(img=list(img),height=img.d[1],width=img.d[2])
}
parse_bcs<-function(path=tissue_paths,scalef=scales$tissue_hires_scalef){
	bcs <- read.csv(path,col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
	bcs$imagerow <- bcs$imagerow * scalef
	bcs$imagecol <- bcs$imagecol * scalef
	bcs$tissue <- as.factor(bcs$tissue)
	bcs
}

createSO<-function(path="MTF"){
	so = CreateSeuratObject(counts = Read10X(paste0(path,"/counts_unfiltered/cellranger/")), assay = "Spatial")
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
	#
	read.table("/opt/tools/spaceranger-1.3.1/lib/python/cellranger/barcodes/visium-v1_coordinates.txt",row.names=1)->coord
	rownames(coord)<-paste(rownames(coord),1,sep="-")
	colnames(coord)<-c("row","col")
	so$row<-coord[colnames(so),"row"]
	so$col<-coord[colnames(so),"col"]
	so
}
plot_feature<-function(object,feature="nCount_Spatial",slot="data",
	alpha=c(0.05, 0.95),limits=c(0.1,0.9),pt.size=2){
	data.df<-object[[]]
	if(! feature %in% colnames(data.df)){
		fdata<-FetchData(object,feature,slot=slot)
		colnames(fdata)<-feature
		data.df<-cbind(data.df,fdata)
	}
	g<-ggplot(data.df,aes(x=row,y=col))+
		geom_point(aes_(fill=as.name(feature),alpha=as.name(feature)),shape=21,stroke=0.2,size=pt.size)+
		ggpubr::theme_pubr(border=T)+
		coord_cartesian(ylim=c(77,0))+labs(x=NULL,y=NULL)+
		#coord_fixed(ylim=c(79,-2),xlim=c(-2,129),expand=F)+labs(x=NULL,y=NULL)+
		scale_fill_gradientn(limits = quantile(data.df[,feature],limits),
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)+
		scale_alpha(range = alpha,guide=NULL)
	list(ggplot=g,data=data.df)
}
plot_QC<-function(object,outfile="QC.pdf"){
	pdf(outfile)
	print(plot_feature(object,feature="nCount_Spatial"))
	print(plot_feature(object,feature="perc_nCount",alpha=c(0.1,1)))
	print(plot_feature(object,feature="nFeature_Spatial"))
	print(plot_feature(object,feature="perc_nFeature",alpha=c(0.1,1)))
	print(plot_feature(object,feature="perc_mt",alpha=c(0.1,1),limits=c(0.05, 0.95)))
	print(plot_feature(object,feature="perc_rs",alpha=c(0.1,1),limits=c(0.05, 0.95)))
	dev.off()
}

#structure(.Data = object, class = "scalefactors")
sf<-scalefactors(spot, fiducial, hires, lowres)
img<-readbitmap::read.bitmap(image_paths)
#
                   # tissue row col imagerow imagecol
# AAACAAGTATCTCCCA-1      1  50 102     7409     8455
# AAACACCAATAACTGC-1      1  59  19     8487     2740
coord<-read.csv(tissue_paths,row.names=1, col.names=c("","tissue","row","col","imagerow","imagecol"), header = FALSE)
Read10X_Image<-function (image.dir, filter.matrix = TRUE, ...) 
{
	image <- readPNG(source = file.path(image.dir, "tissue_lowres_image.png"))
	scale.factors <- fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
	tissue.positions <- read.csv(file = file.path(image.dir, 
		"tissue_positions_list.csv"), col.names = c("barcodes", 
		"tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, 
		as.is = TRUE, row.names = 1)
	if (filter.matrix) {
		tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 1), , drop = FALSE]
	}
	unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
		scale.factors$tissue_lowres_scalef
	spot.radius <- unnormalized.radius/max(dim(x = image))
	return(new(Class = "VisiumV1", 
		image = image, 
		scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
			fiducial = scale.factors$fiducial_diameter_fullres, 
			hires = scale.factors$tissue_hires_scalef, 
			lowres=scale.factors$tissue_lowres_scalef), 
		coordinates = tissue.positions, 
		spot.radius = spot.radius))
}

rjson::fromJSON(file="alignment.json")->alignment
coord<-do.call(rbind,lapply(alignment$oligo,function(x){
	if(length(x)==7) x$tissue<-FALSE
	data.frame(x[c("x","y","row","col","dia","imageX","imageY","tissue")])
}))
fiducial<-do.call(rbind,lapply(alignment$fiducial,function(x){
	if(length(x)==7){x$fidName<-""}
	data.frame(x[c("x","y","row","col","dia","imageX","imageY","fidName")])
}))

library(magick)
image_read("/opt/tools/spaceranger-1.3.1/external/spaceranger_tiny_inputs/image/tinyimage.jpg")->img2
image_write(image_scale(img2,"600"),"test.png",format="png",depth=16)
orig.info<-image_info(img2)
reduced.info<-image_info(image_read("test.png"))
scalef<-min(reduced.info$height/orig.info$width,reduced.info$width/orig.info$width)

parse_image("test.png")->img
pdf("test3.pdf")
ggplot(fiducial,aes(x=imageX*scalef,y=imageY*scalef))+
	geom_spatial(data=img, aes(img=img), x=0.5, y=0.5)+
	geom_point(aes(color=fidName),shape=16,stroke=0,show.legend=F)+
	geom_point(data=subset(coord,tissue),fill="red",shape=21,stroke=0)+
	xlim(0,img$width)+
	ylim(img$height,0)+
	coord_fixed(expand=F)+
	theme_pubr(border=T)+NoAxes()
dev.off()


##
MTF<-createSO("MTF")
PAR<-createSO("Parental")
plot_QC(MTF,"QC.MTF.pdf")
plot_QC(PAR,"QC.Parental.pdf")
DescTools::Gini(MTF$nCount_Spatial)
DescTools::Gini(PAR$nCount_Spatial)
DropletUtils::emptyDrops(PAR@assays$Spatial@counts)->PAR.empty
DropletUtils::emptyDrops(MTF@assays$Spatial@counts)->MTF.empty

VlnPlot(MTF, features = c("perc_mt", "perc_rs"),pt.size=0,ncol=2)
FeatureScatter(MTF, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")

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
MTF[[]]->test2
pdf("Dimreduction.MTF.pdf")
DimPlot(MTF, reduction = "tsne", label = TRUE)
DimPlot(MTF, reduction = "umap", label = TRUE)

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
PAR[[]]->test3
pdf("Dimreduction.Parental.pdf")
DimPlot(PAR, reduction = "tsne", label = TRUE)
DimPlot(PAR, reduction = "umap", label = TRUE)

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

#######################
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmp<- so.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
write.table(tmp,file="marker genes.txt",sep="\t",quote=F)

top10 <- so.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
pdf("Heatmap.marker.pdf",width=5,height=5)
DoHeatmap(so, features = top10$gene,size=3.5,angle=0,hjust=0.5)+theme(axis.text.y=element_text(size=5))
dev.off()

#######
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
Sample<-rep("MTF",ncol(merged))
Sample[stringr::str_detect(colnames(merged),"_2")]<-"PAR"
merged$Sample<-Sample

merged[[]]->test4
pdf("Dimreduction.merged.pdf",width=3.5,height=3)
DimPlot(merged, reduction = "tsne", label = TRUE)
DimPlot(merged, reduction = "tsne", label = TRUE, group.by="Sample")
DimPlot(merged, reduction = "tsne", label = TRUE, split.by="Sample")

label<-as.data.frame(t(sapply(split(test4[,c("row","col")],test4$seurat_clusters),function(x){apply(x,2,median)})))
label$label<-rownames(label)
ggplot(subset(test4,Sample=="MTF"),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F,size=2)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete(guide=NULL)#+guides(fill=guide_legend(nrow=2))
ggplot(subset(test4,Sample=="PAR"),aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1)+
	geom_label(aes(label=label,fill=label),data=label,show.legend=F,size=2)+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL)+
	scale_fill_discrete(guide=NULL)#+guides(fill=guide_legend(nrow=2))
dev.off()

ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==1.1),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster1.1")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)

merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmp <- merged.markers %>% group_by(cluster)
write.table(tmp,file="marker genes.txt",sep="\t",quote=F)

top50 <- merged.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
pdf("Heatmap.marker.pdf",width=5,height=5)
DoHeatmap(merged, features = top50$gene,size=3.5,angle=0,hjust=0.5)+theme(axis.text.y=element_text(size=2))
dev.off()

#####
require(ComplexHeatmap)

sort(merged$seurat_clusters)->clusters
#rlv<-c(1,5,6,8,3,0,2,9,4,7)
rlv<-c(1,5,6,8,0,2,9,4,7,3)
#rlv<-0:9
factor(clusters,levels=rlv)->clusters
color<-scales:::hue_pal()(nlevels(clusters))

top50 <- merged.markers
#top50 <- merged.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
top50<-top50[order(factor(top50$cluster,levels=rlv)),]

merged@assays$SCT@scale.data[,]->dd2
hm<-intersect(intersect(top50$gene,rownames(dd2)),hg38)
mm<-intersect(intersect(top50$gene,rownames(dd2)),mm10)
pdf("test3.pdf",width=5,height=6)
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

####
Heatmap(dd2[c(hm),names(clusters)],
	show_column_names=F,cluster_columns=F,cluster_rows=F,name="hg38",
	col=circlize::colorRamp2(seq(-2,2,length=50),colors=Seurat:::PurpleAndYellow(50)),
	show_row_names=T,row_names_gp=gpar(fontsize=1),row_names_side ="left",
	column_split = clusters,column_gap = unit(0.5, "mm"),#column_title = NULL,
	#row_split= rep(c("hg38","mm10"),c(length(hm),length(mm))),
	column_title_gp = gpar(fontsize=12),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)->h1
Heatmap(dd2[c(mm),names(clusters)],
	show_column_names=F,cluster_columns=F,cluster_rows=F,name="mm10",
	col=circlize::colorRamp2(seq(-2,2,length=50),colors=CustomPalette("magenta","blue","black",50)),
	show_row_names=T,row_names_gp=gpar(fontsize=1),row_names_side ="left",
	column_split = clusters,column_gap = unit(0.5, "mm"),#column_title = NULL,
	#row_split= rep(c("hg38","mm10"),c(length(hm),length(mm))),
	use_raster=T, raster_device="png",raster_quality=8,
	column_title_gp = gpar(fontsize=12)
)->h2
pdf("test4.pdf",width=5,height=6)
draw(h1 %v% h2,gap=unit(1,"mm"))
dev.off()

merged@assays$SCT@data[,]->dd3
hm<-intersect(intersect(top50$gene,rownames(dd3)),hg38)
mm<-intersect(intersect(top50$gene,rownames(dd3)),mm10)
Heatmap(t(scale(t(as.matrix(dd3[c(hm,mm),names(clusters)])))),
	show_column_names=F,cluster_columns=F,cluster_rows=F,name="Expression",
	col=circlize::colorRamp2(seq(-2,2,length=50),colors=Seurat:::PurpleAndYellow(50)),
	show_row_names=T,row_names_gp=gpar(fontsize=1),row_names_side ="left",
	column_split = clusters,column_gap = unit(0.5, "mm"),#column_title = NULL,
	row_split= rep(c("hg38","mm10"),c(length(hm),length(mm))),
	column_title_gp = gpar(fontsize=12),
	top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
#####
MTF2<-SCTransform(MTF,assay = "Spatial", verbose = FALSE,return.only.var.genes=F)
PAR2<-SCTransform(PAR,assay = "Spatial", verbose = FALSE,return.only.var.genes=F)
MTF@assays$SCT@scale.data<-MTF2@assays$SCT@scale.data
PAR@assays$SCT@scale.data<-PAR2@assays$SCT@scale.data

MTF.ss <- subset(MTF, idents = c(3,8))
PAR.ss <- subset(PAR, idents = c(3,5:8,10))
merged <- merge(MTF.ss, PAR.ss)


#####
pdf("Mouse_ratio.pdf")
plot_feature(merged[,merged$Sample=="MTF"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="MTF")
plot_feature(merged[,merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="PAR")
plot_feature(merged[,merged$seurat_clusters==0 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster0")
plot_feature(merged[,merged$seurat_clusters==1 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster1")
plot_feature(merged[,merged$seurat_clusters==2 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster2")
plot_feature(merged[,merged$seurat_clusters==3 & merged$Sample=="MTF"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster3")
plot_feature(merged[,merged$seurat_clusters==4 & merged$Sample=="MTF"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster4")
plot_feature(merged[,merged$seurat_clusters==5 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster5")
plot_feature(merged[,merged$seurat_clusters==6 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster6")
plot_feature(merged[,merged$seurat_clusters==7 & merged$Sample=="MTF"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster7")
plot_feature(merged[,merged$seurat_clusters==8 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster8")
plot_feature(merged[,merged$seurat_clusters==9 & merged$Sample=="PAR"],"perc_nCount",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster9")
dev.off()

merged$ratio<-colSums(GetAssayData(merged[rownames(merged) %in% mm10,],assay="SCT",slot="data"))/colSums(GetAssayData(merged,assay="SCT",slot="data"))
merged$hg38<- 1-merged$ratio
pdf("Mouse_ratio2.pdf")
plot_feature(merged[,merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="MTF")
plot_feature(merged[,merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="PAR")
plot_feature(merged[,merged$seurat_clusters==0 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster0")
plot_feature(merged[,merged$seurat_clusters==1 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster1")
plot_feature(merged[,merged$seurat_clusters==2 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster2")
plot_feature(merged[,merged$seurat_clusters==3 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster3")
plot_feature(merged[,merged$seurat_clusters==4 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster4")
plot_feature(merged[,merged$seurat_clusters==5 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster5")
plot_feature(merged[,merged$seurat_clusters==6 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster6")
plot_feature(merged[,merged$seurat_clusters==7 & merged$Sample=="MTF"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster7")
plot_feature(merged[,merged$seurat_clusters==8 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster8")
plot_feature(merged[,merged$seurat_clusters==9 & merged$Sample=="PAR"],"ratio",alpha=c(1,1))$ggplot+xlim(c(0,127))+labs(fill="Perc",title="Cluster9")
dev.off()
pdf("Mouse_ratio.v3.pdf",width=3.5,height=4)
plot_feature(merged[,merged$Sample=="MTF" & merged$seurat_clusters %in% c(4,7)],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=CustomPalette("white","darkblue",k=20),oob = scales::squish)
plot_feature(merged[,merged$Sample=="MTF" & merged$seurat_clusters %in% c(4,7)],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=CustomPalette("white","darkred",k=20),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR" & merged$seurat_clusters %in% c(2,5,8)],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=CustomPalette("white","darkblue",k=20),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR" & merged$seurat_clusters %in% c(2,5,8)],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=CustomPalette("white","darkred",k=20),oob = scales::squish)
dev.off()

pdf("Mouse_ratio.v4.pdf",width=3.5,height=4)
plot_feature(merged[,merged$Sample=="MTF" & merged$seurat_clusters %in% c(4,7)],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="MTF" & merged$seurat_clusters %in% c(4,7)],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR" & merged$seurat_clusters %in% c(2,5,8)],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR" & merged$seurat_clusters %in% c(2,5,8)],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
dev.off()

pdf("Mouse_ratio.v5.pdf",width=3.5,height=4)
plot_feature(merged[,merged$Sample=="MTF"],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="MTF"],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR"],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR"],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
dev.off()


plot_feature_gene<-function(object=merged,gene="KLK3",limits=NULL){
	so1<-object[,object$Sample=="MTF" & object$seurat_clusters %in% c(4,7)]
	so2<-object[,object$Sample=="PAR" & object$seurat_clusters %in% c(2,5,8)]
	rg<-range(c(FetchData(so1,gene)[,1],FetchData(so2,gene)[,1]))
	limits<-limits %||% round(rg,2)
	p1<-plot_feature(so1,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="MTF")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	p2<-plot_feature(so2,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="PAR")+
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

genes2<-unique(read.table("genes2.txt",header=F)$V1)
genes2<-mm10[which(str_to_upper(mm10) %in% str_to_upper(genes2))]
pdf("Map.genes.v2.pdf",height=4,width=7)
lapply(genes2,function(x){
	plot_feature_gene(merged,x)
	NULL
})
dev.off()

#
ratio_df<-merged[[]][,c("perc_nCount","ratio","seurat_clusters")]
ratio_df$seurat_clusters<-factor(ratio_df$seurat_clusters,rlv)
color<-scales:::hue_pal()(nlevels(ratio_df$seurat_clusters))
pdf("Ratio.beeswarm.pdf",width=4,height=2)
ggplot(ratio_df,aes(seurat_clusters,perc_nCount*100))+
	#geom_violin()+
	geom_hline(yintercept=c(25,50,75),linetype=2,col="grey")+
	geom_beeswarm(aes(col=seurat_clusters),priority="random",
		cex=0.3,dodge.width=0,size=0.1,show.legend=F)+
	theme_pubr()+
	scale_color_manual(values=color[rlv+1])+
	labs(x=NULL,y="Percentage (%)")
ggplot(ratio_df,aes(seurat_clusters,ratio*100))+
	#geom_violin()+
	geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
	geom_beeswarm(aes(col=seurat_clusters),priority="random",
		cex=0.3,dodge.width=0,size=0.1,show.legend=F)+
	theme_pubr()+
	scale_color_manual(values=color[rlv+1])+
	labs(x=NULL,y="Percentage (%)")
dev.off()
#
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
xmin=seq(0.5,by=1,length=10)
xmax=seq(1.5,by=1,length=10)
fill=rep_len(c("white","blue"),10)
alpha=rep_len(c(0,0.05),10)

pdf("Ratio.beeswarm.v2.pdf",width=4.8,height=2)
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

#Ridge plot
library(ggridges)
pdf("Ridge.perc.pdf",width=4.5,height=3)
ggplot(merged[[]],aes(y=seurat_clusters,x=perc_nCount))+
	geom_density_ridges(aes(fill=seurat_clusters),show.legend=F,rel_min_height = 0.01,scale=2)+
	scale_x_continuous(expand = c(0, 0)) +
	scale_y_discrete(expand = c(0, 0)) +
	#coord_cartesian(clip = "off") + 
	theme_ridges(grid = T, center_axis_labels = TRUE)+
	#scale_x_continuous(trans=scales::log10_trans(),n.breaks=5)+
	labs(x="Mouse_perc",y=NULL)
dev.off()

####
genes<-c("PTPRC","EPCAM","MUC1","FOLH1","KLK3","TMPRSS2","NKX3-1","AR","VIM","FN1","HLA-A","HLA-C")
FetchData(merged,genes)->tmp
cbind(merged[[]],tmp)->genes.df
pdf("Violin2.genes.pdf",width=4,height=3)
lapply(genes,function(x){
	g<-ggplot(genes.df,aes_(x=~seurat_clusters,y=as.name(x),fill=~seurat_clusters))+
		geom_violin(scale="width")+ggpubr::theme_pubr(base_size=16)+NoLegend()+
		labs(x="identity",y="Expression level",title=x)
	print(g)
	NULL
})
dev.off()

pdf("Violin.genes.pdf",width=4,height=3)
lapply(genes,function(x){
	g<-VlnPlot(merged,x,stack=F,pt.size=0)+NoLegend()+
		theme(axis.text.x=element_text(angle=0,hjust=0.5))
	print(g)
	NULL
})
dev.off()

#
read.table("MTF/counts_unfiltered/cellranger/genes.tsv",header=F,sep="\t")->symbol
unique(symbol[str_detect(symbol$V1,"^ENSG"),"V2"])->hg38
setdiff(symbol$V2,hg38)->mm10
cluster.genes<-split(merged.markers$gene,merged.markers$cluster)
cluster.genes.mm10<-lapply(cluster.genes,function(x){intersect(mm10,x)})
cluster.genes.hg38<-lapply(cluster.genes,function(x){intersect(hg38,x)})

require(enrichR)
require(xlsx)

dbs<-sort(enrichR::listEnrichrDbs()$libraryName)[c(67,104,129,146,160)]
enriched.hg38 <- lapply(cluster.genes.hg38,function(x){
	enrichr(x, dbs)
})
lapply(names(enriched.hg38),function(x){
	lapply(names(enriched.hg38[[x]]),function(y){
		if(nrow(enriched.hg38[[x]][[y]] > 0)){
			write.xlsx(enriched.hg38[[x]][[y]],file=paste0("hg38.cluster",x,".xlsx"),sheetName=y,append=T)
		}
		NULL
	})
	NULL
})

enriched.mm10 <- lapply(cluster.genes.mm10,function(x){
	enrichr(x, dbs)
})
lapply(names(enriched.mm10),function(x){
	lapply(names(enriched.mm10[[x]]),function(y){
		if(nrow(enriched.mm10[[x]][[y]] > 0)){
			write.xlsx(enriched.mm10[[x]][[y]],file=paste0("mm10.cluster",x,".xlsx"),sheetName=y,append=T)
		}
		NULL
	})
	NULL
})

##
load("/data/genome/enrichr/.RData")
#parse_geneset_term("~/work/MSigDB/mouse/c2.cgp.v6.1.symbols_mouse.gmt")->m.cgp
#GO, KEGG_human, KEGG_mouse, Panther, "Reactome", "Hallmark"
#enrichr
run_GO<-function(genes,gs,n=20000){
	k<-length(genes)
	sapply(names(gs),function(go){
		m<-length(gs[[go]])
		overlap<-intersect(gs[[go]],genes)
		x<-length(overlap)
		if(x==0){
			c(1,0,0,m,"")
		}
		else{
			#ftest<-fisher.test(matrix(c(x, m-x, k-x, n-m-(k-x)),nrow=2,byrow=F), alternative='greater')
			#c(ftest$p.value,ftest$estimate,x,m)
			c(phyper(x-1,m,n-m,k,lower.tail=F),x/(k-x)/((m-x)/(n-k+x)),x,m,paste(overlap,collapse=";"))
		}
	})->tmp
	as.data.frame(t(tmp))->tmp
	tmp[, 1:4] <- sapply(tmp[, 1:4], as.numeric)
	tmp<-tmp[order(tmp[,1]),]
	subset(tmp,V3!=0)->tmp
	p.adjust(tmp$V1,method="BH")->tmp$V6
	colnames(tmp)<-c("pval","odds.ratio","overlap","gs.size","genes","padj")
	tmp
}
run_GO(cluster.genes.hg38$`3`,GO)->tmp

###############################
###############################
#single cell level pathway enrichment
library(AUCell)
generate_random<-function(gs,genes=hg38,n=1000){
	size<-sapply(gs,length)
	lapply(names(size),function(x){
		gs.list<-replicate(n,sample(genes,size[x]),simplify=F)
		names(gs.list)<-paste("random",1:n)
		gs.list
	})->tmp
	names(tmp)<-names(size)
	tmp
}

generate_random(KEGG_human,genes=hg38,n=1000)->rd1000
AUCell_buildRankings(merged.hg38@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
lapply(rd1000,function(x){
	AUCell_calcAUC(x,cells_rankings,aucMaxRank=500,nCores=32,verbose = FALSE)->cells_AUC
	cells_AUC@assays@data$AUC
})->rd1000.res

#hg38
merged[rownames(merged) %in% hg38,]->merged.hg38
AUCell_buildRankings(merged.hg38@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
AUCell_calcAUC(KEGG_human,cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, assign=TRUE,nCores=32) 
AUC<-cells_AUC@assays@data$AUC
merged.hg38[["AUC"]]<-CreateAssayObject(data=AUC)

# sapply(rd1000.res,function(x){apply(x,2,mean)})->tmp1
# sapply(rd1000.res,function(x){apply(x,2,sd)})->tmp2
# (AUC-t(tmp1))/t(tmp2)->AUC.norm
# merged.hg38[["AUCn"]]<-CreateAssayObject(data=AUC.norm)


plot_feature(merged.hg38,"Phagosome")$ggplot

FeaturePlot(merged.hg38,reduction="tsne",features="Phagosome")+
	scale_color_gradientn(colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)

sort(merged.hg38$seurat_clusters)->clusters
rlv<-c(1,5,6,8,0,2,9,4,7,3)
factor(clusters,levels=rlv)->clusters
color<-scales:::hue_pal()(nlevels(clusters))
#t(scale(t(merged@assays$AUC@data[,names(clusters)])))

merged.hg38@assays$AUC@data[,names(clusters)]->auc.mat
mat<-(sweep(auc.mat,1,apply(auc.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
#o = seriate(mat, method = "BEA_TSP")
pdf("test.pdf",height=7,width=5)
Heatmap(mat,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.2),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::PurpleAndYellow(50)),
	cluster_columns=F,
	clustering_distance_rows="pearson",
	#row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.5, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()

#
AUCell_calcAUC(Hallmark,cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
merged.hg38[["Hallmark"]]<-CreateAssayObject(data=cells_AUC@assays@data$AUC)

merged.hg38@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))

do.call(rbind,lapply(levels(clusters),function(y){
	lapply(1:50,function(x){
		t.test(mat2[x,]~(clusters==y))->tmp
		c(tmp$p.val,tmp$est[2]-tmp$est[1],as.numeric(y))
	})->tmp
	do.call(rbind,tmp)
}))->tmp
as.data.frame(tmp)->tmp
colnames(tmp)<-c("pval","diff","cluster")
tmp$gs<-rownames(mat2)
tmp1<-dcast(tmp,gs~cluster,value.var="diff")
rownames(tmp1)<-tmp1[,1];tmp1<-tmp1[,-1]
tmp2<-dcast(tmp,gs~cluster,value.var="pval")
rownames(tmp2)<-tmp2[,1];tmp2<-tmp2[,-1]


o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.pdf",height=4,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.4),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,
	clustering_distance_rows="pearson",
	row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()


merged.hg38@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
mat2<-mat2[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49),]
o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.v2.pdf",height=3,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,
	clustering_distance_rows="pearson",
	row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()

#
pw<-c("Interferon Gamma Response","Inflammatory Response","TNF-alpha Signaling via NF-kB")
pdf("Map.pathway.human.pdf",height=4,width=7)
lapply(pw,function(x){
	plot_feature_gene(merged.hg38,x)
	NULL
})
dev.off()

pdf("Map.pathway.mouse.pdf",height=4,width=7)
lapply(pw,function(x){
	plot_feature_gene(merged.mm10,x)
	NULL
})
dev.off()

plot_feature_gene2<-function(object=merged,gene="KLK3",limits=NULL){
	so1<-object[,object$Sample=="MTF"]
	so2<-object[,object$Sample=="PAR"]
	rg<-range(c(FetchData(so1,gene)[,1],FetchData(so2,gene)[,1]))
	limits<-limits %||% round(rg,2)
	p1<-plot_feature(so1,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="MTF")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	p2<-plot_feature(so2,gene,alpha=c(1,1),pt.size=0.8)$ggplot+
		xlim(c(0,127))+labs(fill=gene,title="PAR")+
		scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
	gridExtra::grid.arrange(p1,p2,ncol=2)
}

pw<-c("Interferon Gamma Response","Inflammatory Response","TNF-alpha Signaling via NF-kB")
pdf("Map.pathway.human.all.pdf",height=4,width=7)
lapply(pw,function(x){
	plot_feature_gene2(merged.hg38,x)
	NULL
})
dev.off()

pdf("Map.pathway.mouse.all.pdf",height=4,width=7)
lapply(pw,function(x){
	plot_feature_gene2(merged.mm10,x)
	NULL
})
dev.off()

#pathway by beeswarm
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)

coeff<-sapply(1:3,function(x)max(ratio_df2[,x])/max(ratio_df[,x]))
ratio_df2[,1:3]<-sweep(ratio_df2[,1:3],2,coeff,"/")

ratio_df<-reshape2::melt(ratio_df)
ratio_df$group<-"hg38"
ratio_df2<-reshape2::melt(ratio_df2)
ratio_df2$group<-"mm10"
rbind(ratio_df,ratio_df2)->ratio_df3
ratio_df3$group<-factor(ratio_df3$group,levels=c("mm10","hg38"))
#quasirandom, pseudorandom, smiley or frowney

pdf("Pathway.beeswarm.pdf",width=5,height=2.3)
lapply(1:3,function(x){
	df<-subset(ratio_df3,variable==paste0("hallmark_",pw[x]))
	g<-ggplot(df,aes(seurat_clusters,value,color=group))+
		geom_beeswarm(aes(group=group),priority="random",method="frowney",
			cex=1,dodge.width=0,size=0.1,show.legend=T,position=position_dodge2()
		)
	tmp<-ggplot_build(g)
	yrange<-tmp$layout$panel_scales_y[[1]]$range$range
	yrange[2]<-yrange[2]*1.05
	yrange[1]<-yrange[1]*0.95
	xmin=seq(0.5,by=1,length=10)
	xmax=seq(1.5,by=1,length=10)
	fill=rep_len(c("white","blue"),10)
	alpha=rep_len(c(0,0.05),10)
	g<-ggplot(df,aes(seurat_clusters,value,color=factor(group)))+
		annotate("rect",xmin=xmin,xmax=xmax,ymin=yrange[1],ymax=yrange[2],fill=fill,alpha=alpha)+
		scale_x_discrete()+
		scale_y_continuous(name="mm10",sec.axis=sec_axis(~.*coeff[x], name="hg38"))+
		#geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
		geom_beeswarm(aes(group=group),priority="random",method="frowney",
			cex=0.4,dodge.width=0.5,size=0.1,show.legend=T,position=position_dodge2()
		)+
		theme_pubr(border=T)+
		theme(legend.position="right",axis.text.y.left= element_text(color="blue"),axis.text.y.right= element_text(color="red"),
			axis.ticks.y.left= element_line(color="blue"),axis.ticks.y.right= element_line(color="red"))+
		scale_color_manual(breaks=c("mm10","hg38"),values=c("blue","red"),labels=c("mm10","hg38"))+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		guides(y = guide_axis(position="right"), y.sec = guide_axis(position="left"))+
		labs(x=NULL,y="AUC score",color=NULL,title=pw[x]) + 
		coord_cartesian(expand=F)
	print(g)
	NULL
})
dev.off()

#mm10
parse_geneset_term("~/work/MSigDB/mouse/h.all.v6.1.symbols_mouse.gmt")->hallmark.mm10
names(hallmark.mm10)<-names(Hallmark)
merged[rownames(merged) %in% mm10,]->merged.mm10
AUCell_buildRankings(merged.mm10@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings2
AUCell_calcAUC(hallmark.mm10,cells_rankings2,aucMaxRank=500,nCores=32)->cells_AUC2 
merged.mm10[["Hallmark"]]<-CreateAssayObject(data=cells_AUC2@assays@data$AUC)

merged.mm10@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.mouse.pdf",height=4,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.4),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,
	clustering_distance_rows="pearson",
	row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()


merged.mm10@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
mat2<-mat2[c(40,36,1,31,7,18,19,28,26,30,23,10,34,15,25,48,37,22,2,35,43),]
o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.mouse.v2.pdf",height=3,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,
	clustering_distance_rows="pearson",
	row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		#labels_gp = gpar(col = "white", fontsize = 10)
		),
		height=unit(3,"mm")
	)
)
dev.off()

###############################
#Trendsceek
so <- FindSpatiallyVariableFeatures(so, assay = "SCT", features = VariableFeatures(so),
    selection.method = "markvariogram")

top.features <- head(SpatiallyVariableFeatures(so, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(so, features = top.features, ncol = 3, alpha = c(0.1, 1))
#cell number distribution
pdf("Cell_number.pdf",width=2.3,height=3)
ggplot(so@meta.data,aes(fill=group,x=seurat_clusters))+geom_bar(width=0.5)+theme_pubr(border=T)+ylab("Cell count")
dev.off()



######################
#3.1 and 3.2
merged<-StashIdent(merged,save.name="old")
as.character(merged$seurat_clusters)->clusters
clusters[which(merged$seurat_clusters==3&merged$Sample=="MTF")]<-"3.1"
clusters[which(merged$seurat_clusters==3&merged$Sample=="PAR")]<-"3.2"

rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
factor(clusters,levels=rlv)->clusters
Idents(merged)<-clusters
merged$seurat_clusters<-clusters
names(clusters)<-colnames(merged)

compare.c3.1 <- FindMarkers(merged, only.pos = TRUE,assay="Spatial",ident.1 = "3.1",ident.2=NULL,min.pct = 0.25, logfc.threshold = 0.25)
compare.c3.2 <- FindMarkers(merged, only.pos = TRUE,assay="Spatial",ident.1 = "3.2",ident.2=NULL,min.pct = 0.25, logfc.threshold = 0.25)
tmp <- merged.markers %>% group_by(cluster)

###
custom1<-parse_geneset_term("custom.gs3.txt") #human
custom2<-parse_geneset_term("custom.gs4.txt") #mouse

merged[rownames(merged) %in% hg38,]->merged.hg38
AUCell_buildRankings(merged.hg38@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
AUCell_calcAUC(c(Hallmark,custom1),cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
merged.hg38[["Hallmark"]]<-CreateAssayObject(data=cells_AUC@assays@data$AUC)

parse_geneset_term("~/work/MSigDB/mouse/h.all.v6.1.symbols_mouse.gmt")->hallmark.mm10
names(hallmark.mm10)<-names(Hallmark)
merged[rownames(merged) %in% mm10,]->merged.mm10
AUCell_buildRankings(merged.mm10@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings2
AUCell_calcAUC(c(hallmark.mm10,custom2),cells_rankings2,aucMaxRank=500,nCores=32)->cells_AUC2 
merged.mm10[["Hallmark"]]<-CreateAssayObject(data=cells_AUC2@assays@data$AUC)


#pathway by beeswarm
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)

coeff<-sapply(1:24,function(x)max(ratio_df2[,x])/max(ratio_df[,x]))
ratio_df2[,1:24]<-sweep(ratio_df2[,1:24],2,coeff,"/")

ratio_df<-reshape2::melt(ratio_df)
ratio_df$group<-"hg38"
ratio_df2<-reshape2::melt(ratio_df2)
ratio_df2$group<-"mm10"
rbind(ratio_df,ratio_df2)->ratio_df3
ratio_df3$group<-factor(ratio_df3$group,levels=c("mm10","hg38"))
#quasirandom, pseudorandom, smiley or frowney
pdf("Pathway.beeswarm.pdf",width=4.5,height=2.3)
lapply(1:24,function(x){
	df<-subset(ratio_df3,value>0&seurat_clusters %in% c(3.1,3.2)&variable==paste0("hallmark_",pw[x]))
	g<-ggplot(df,aes(seurat_clusters,value,color=factor(group)))+
		scale_x_discrete()+
		scale_y_continuous(name="mm10",sec.axis=sec_axis(~.*coeff[x], name="hg38"))+
		#geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
		geom_violin(scale="area",position=position_dodge(width=0.8))+
		geom_boxplot(show.legend=F,outlier.shape=NA,position=position_dodge(width=0.8),width=0.3)+
		geom_beeswarm(aes(group=group),priority="random",method="frowney",
			cex=1,dodge.width=0.8,size=0.2,show.legend=F,alpha=0.3
		)+
		theme_pubr(border=T)+
		theme(legend.position="right",axis.text.y.left= element_text(color="blue"),axis.text.y.right= element_text(color="red"),
			axis.ticks.y.left= element_line(color="blue"),axis.ticks.y.right= element_line(color="red"))+
		scale_color_manual(breaks=c("mm10","hg38"),values=c("blue","red"),labels=c("mm10","hg38"))+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		guides(y = guide_axis(position="right"), y.sec = guide_axis(position="left"))+
		labs(x=NULL,y="AUC score",color=NULL,title=pw[x]) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()

#
rlv<-c(1,5,6,8,0,2,9,4,7,3)
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
FetchData(merged.hg38,pw)->ratio_df
ratio_df$old<-factor(merged.hg38$old,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$old<-factor(merged.mm10$old,rlv)
#
m1<-do.call(cbind,lapply(split(ratio_df[,-25],ratio_df$old),function(x)apply(x,2,mean)))
m2<-do.call(cbind,lapply(split(ratio_df2[,-25],ratio_df2$old),function(x)apply(x,2,mean)))
rownames(m1)<-rownames(m2)<-pw
lvs<-names(sort(apply(m1[,8:10],1,sum)))
g1<-ggplot(reshape2::melt(m1))+geom_point(aes(x=factor(Var2,levels=rlv),y=factor(Var1,levels=lvs),size=value),shape=19,color="#A50026")+
	labs(x=NULL,y=NULL)+
	scale_x_discrete(position="top")+
	scale_size_area(max_size=4,breaks=c(0.05,0.1,0.15),limits=c(0,0.15),oob=scales::squish)+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.y= element_text(size=8),axis.ticks=element_blank())
g2<-ggplot(reshape2::melt(m2))+geom_point(aes(x=factor(Var2,levels=rlv),y=factor(Var1,levels=lvs),size=value),shape=19,color="#313695")+
	labs(x=NULL,y=NULL)+
	scale_y_discrete(labels=NULL,breaks=NULL)+
	scale_x_discrete(position="top")+
	scale_size_area(max_size=4,breaks=c(0.04,0.08,0.12),limits=c(0,0.12),oob=scales::squish)+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.y= element_text(size=8),axis.ticks=element_blank())
pdf("Pathway.circle.pdf",height=4,width=7)
g1 %+% g2
dev.off()


rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)

ratio_df<-subset(ratio_df,seurat_clusters %in% c(3.1,3.2))
ratio_df2<-subset(ratio_df2,seurat_clusters %in% c(3.1,3.2))

m3<-do.call(cbind,lapply(split(ratio_df[,-25],ratio_df$seurat_clusters),function(x)apply(x,2,mean)))[,10:11]
m4<-do.call(cbind,lapply(split(ratio_df2[,-25],ratio_df2$seurat_clusters),function(x)apply(x,2,mean)))[,10:11]
rownames(m3)<-rownames(m4)<-pw
g3<-ggplot(reshape2::melt(m3))+geom_point(aes(x=factor(Var2,levels=c(3.1,3.2)),y=factor(Var1,levels=lvs),size=value),shape=19,color="#A50026")+
	labs(x=NULL,y=NULL)+scale_x_discrete(position="top")+
	scale_size(range=c(0,4),breaks=c(0.05,0.1,0.15),limits=c(0,0.15))+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.y= element_text(size=8),axis.ticks=element_blank())
g4<-ggplot(reshape2::melt(m4))+geom_point(aes(x=factor(Var2,levels=c(3.1,3.2)),y=factor(Var1,levels=lvs),size=value),shape=19,color="#313695")+
	labs(x=NULL,y=NULL)+scale_y_discrete(labels=NULL,breaks=NULL)+scale_x_discrete(position="top")+
	scale_size(range=c(0,4),breaks=c(0.04,0.08,0.12),limits=c(0,0.12))+
	theme_pubr(border=T)+
	theme(legend.position="right",axis.text.y= element_text(size=8),axis.ticks=element_blank())
pdf("Cluster3.pathway.circle.pdf",height=4,width=5)
g3 %+% g4
dev.off()


#
coeff<-sapply(1:24,function(x)max(ratio_df2[,x])/max(ratio_df[,x]))
ratio_df2[,1:24]<-sweep(ratio_df2[,1:24],2,coeff,"/")

ratio_df<-reshape2::melt(ratio_df)
ratio_df$group<-"hg38"
ratio_df2<-reshape2::melt(ratio_df2)
ratio_df2$group<-"mm10"
rbind(ratio_df,ratio_df2)->ratio_df3
ratio_df3$group<-factor(ratio_df3$group,levels=c("mm10","hg38"))
ratio_df3$xaxis<-paste(ratio_df3$group,ratio_df3$old,sep=":")
ratio_df3$xaxis<-factor(ratio_df3$xaxis,paste(rep(c("mm10","hg38"),each=10),rlv,sep=":"))
pdf("Pathway.beeswarm.v2.pdf",width=5.5,height=2.3)
lapply(1:24,function(x){
	df<-subset(ratio_df3,value>0&variable==paste0("hallmark_",pw[x]))
	g<-ggplot(df,aes(factor(xaxis),value,color=group))+
		scale_x_discrete(labels=rep(rlv,2))+
		scale_y_continuous(name="mm10",sec.axis=sec_axis(~.*coeff[x], name="hg38"))+
		#geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
		geom_violin(scale="area",size=0.2)+
		#geom_boxplot(show.legend=F,outlier.shape=NA,width=0.2)+
		geom_beeswarm(priority="random",method="frowney",
			cex=0.4,size=0.1,show.legend=F,alpha=0.2
		)+
		theme_pubr(border=T)+
		theme(legend.position="right",axis.text.y.left= element_text(color="blue"),axis.text.y.right= element_text(color="red"),
			axis.ticks.y.left= element_line(color="blue"),axis.ticks.y.right= element_line(color="red"))+
		scale_color_manual(breaks=c("mm10","hg38"),values=c("blue","red"),labels=c("mm10","hg38"))+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		guides(y = guide_axis(position="right"), y.sec = guide_axis(position="left"))+
		labs(x=NULL,y="AUC score",color=NULL,title=pw[x]) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()

####################
MTF[[]]->test2
PAR[[]]->test3
merged[[]]->test4
test2$seurat_clusters<-0
test3$seurat_clusters<-0
test2[str_replace(rownames(subset(test4,seurat_clusters==3.1)),"_1",""),"seurat_clusters"]<-3.1
test3[str_replace(rownames(subset(test4,seurat_clusters==3.2)),"_2",""),"seurat_clusters"]<-3.2
test2$seurat_clusters<-test2$seurat_clusters==3.1
test3$seurat_clusters<-test3$seurat_clusters==3.2
pdf("cluster3.pdf",width=3.5,height=3)
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster3.1")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster3.2")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
dev.off()

#######
merged$ratio<-colSums(GetAssayData(merged[rownames(merged) %in% mm10,],assay="SCT",slot="data"))/colSums(GetAssayData(merged,assay="SCT",slot="data"))
merged$hg38<- 1-merged$ratio
pdf("Mouse_ratio.v6.pdf",width=3.5,height=4)
plot_feature(merged[,merged$Sample=="MTF"&merged$seurat_clusters==3.1],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="MTF"&merged$seurat_clusters==3.1],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR"&merged$seurat_clusters==3.2],"ratio",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="mm10")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
plot_feature(merged[,merged$Sample=="PAR"&merged$seurat_clusters==3.2],"hg38",alpha=c(1,1),pt.size=0.8)$ggplot+
	xlim(c(0,127))+labs(fill="Perc",title="hg38")+
	scale_fill_gradientn(limits=c(0,1),
		colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
dev.off()

#######
merged$ratio<-colSums(GetAssayData(merged[rownames(merged) %in% mm10,],assay="SCT",slot="data"))/colSums(GetAssayData(merged,assay="SCT",slot="data"))
ratio_df<-merged[[]][,c("perc_nCount","ratio","seurat_clusters")]
ratio_df$seurat_clusters<-factor(ratio_df$seurat_clusters,rlv)
color<-scales:::hue_pal()(nlevels(ratio_df$seurat_clusters))
ratio_df$hg38<- 1-ratio_df$ratio
ratio_df2<-reshape2::melt(ratio_df[,c(3,2,4)])
#quasirandom, pseudorandom, smiley or frowney

g<-ggplot(subset(ratio_df2,seurat_clusters%in%c(3.1,3.2)),aes(seurat_clusters,value*100,color=variable))+
	#geom_violin()+
	geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
	geom_beeswarm(aes(group=variable),priority="random",method="frowney",
		cex=0.3,dodge.width=0.5,size=0.1,show.legend=T,position=position_dodge2()
	)
tmp<-ggplot_build(g)
yrange<-tmp$layout$panel_scales_y[[1]]$range$range
yrange[2]<-yrange[2]+5
yrange[1]<-yrange[1]-2
xmin=seq(0.5,by=1,length=10)
xmax=seq(1.5,by=1,length=10)
fill=rep_len(c("white","blue"),10)
alpha=rep_len(c(0,0.05),10)

pdf("cluster3.Ratio.pdf",width=2.5,height=2)
ggplot(subset(ratio_df2,seurat_clusters%in%c(3.1,3.2)),aes(seurat_clusters,value*100,color=variable))+
	scale_x_discrete()+
	geom_hline(yintercept=c(25,50,75,100),linetype=2,col="grey")+
	geom_beeswarm(aes(group=variable),priority="random",method="frowney",
		cex=0.6,dodge.width=0.6,size=0.1,show.legend=T,position=position_dodge2())+
	theme_pubr(border=T)+theme(legend.position="right")+
	scale_color_manual(values=c("blue","red"),labels=c("mm10","hg38"))+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="Percentage (%)",color=NULL) + 
	coord_cartesian(expand=F)
dev.off()

##########
rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
pw<-c("M1 marker","M2 marker")
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)
colnames(ratio_df)<-colnames(ratio_df2)<-c("M1","M2","seurat_clusters")

pdf("cluster3.M1_M2.pdf",height=3,width=3.5)
ggplot(subset(ratio_df,seurat_clusters %in% c(3.1)),aes(M1,M2))+
	geom_jitter(shape=21,color="red")+theme_pubr(border=T)+
	labs(title="Human 3.1")+
	coord_cartesian(xlim=c(0,0.2),ylim=c(0,0.25))
ggplot(subset(ratio_df,seurat_clusters %in% c(3.2)),aes(M1,M2))+
	geom_jitter(shape=21,color="red")+theme_pubr(border=T)+
	labs(title="Human 3.2")+
	coord_cartesian(xlim=c(0,0.2),ylim=c(0,0.25))
ggplot(subset(ratio_df2,seurat_clusters %in% c(3.1)),aes(M1,M2))+
	geom_jitter(shape=21,color="blue")+theme_pubr(border=T)+
	labs(title="Mouse 3.1")+
	coord_cartesian(xlim=c(0,0.2),ylim=c(0,0.25))
ggplot(subset(ratio_df2,seurat_clusters %in% c(3.2)),aes(M1,M2))+
	geom_jitter(shape=21,color="blue")+theme_pubr(border=T)+
	labs(title="Mouse 3.2")+
	coord_cartesian(xlim=c(0,0.2),ylim=c(0,0.25))
dev.off()

#########################

pw<-names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)]
FetchData(merged.hg38,pw)->ratio_df
ratio_df$old<-factor(merged.hg38$old,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$old<-factor(merged.mm10$old,rlv)

rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
pw<-names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)]
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)
colnames(ratio_df)<-c(paste0("H_",pw),"seurat_clusters")
colnames(ratio_df2)<-c(paste0("M_",pw),"seurat_clusters")
ratio_df3<-cbind(ratio_df,ratio_df2)

pdf("cluster3.IFNg.pdf",height=3,width=3.5)
lapply(pw,function(x){
	g<-ggplot(subset(ratio_df3,seurat_clusters %in% c(3.1,3.2)),aes_string(as.name(paste0("H_",x)),as.name(paste0("M_",x))))+
		geom_jitter(aes(color=seurat_clusters),shape=21,show.legend=F)+theme_pubr(border=T)+
		labs(title="3.1 and 3.2")+
		scale_color_manual(values=c("green","purple"))#+
		#coord_cartesian(xlim=c(0,0.05),ylim=c(0,0.06))
	print(g)
	NULL
})
dev.off()

ratio_df4<-subset(ratio_df3,seurat_clusters %in% c(3.1,3.2))
c1<-data.frame(row.names=c("I","II","III","IV"))
c2<-data.frame(row.names=c("I","II","III","IV"))
pdf("cluster3.IFNg.v2.pdf",height=3,width=3.5)
lapply(pw,function(x){
	xin<-median(ratio_df4[,paste0("H_",x)])
	yin<-median(ratio_df4[,paste0("M_",x)])
	tmp1<-subset(ratio_df4,seurat_clusters==3.1)
	tmp2<-subset(ratio_df4,seurat_clusters==3.2)
	c1[,x]<<-table(tmp1[,paste0("H_",x)]>xin,tmp1[,paste0("M_",x)]>yin)[c(4,2,3,1)]
	c2[,x]<<-table(tmp2[,paste0("H_",x)]>xin,tmp2[,paste0("M_",x)]>yin)[c(4,2,3,1)]
	g<-ggplot(ratio_df4,aes_string(as.name(paste0("H_",x)),as.name(paste0("M_",x))))+
		geom_jitter(aes(color=seurat_clusters),shape=21,show.legend=F)+theme_pubr(border=T)+
		geom_vline(xintercept=xin,linetype=2)+geom_hline(yintercept=yin,linetype=2)+
		labs(title="3.1 and 3.2")+
		scale_color_manual(values=c("green","purple"))#+
		#coord_cartesian(xlim=c(0,0.05),ylim=c(0,0.06))
	print(g)
	NULL
})
dev.off()

c1$group<-rownames(c1)
c2$group<-rownames(c2)
c1<-reshape2::melt(c1)
c2<-reshape2::melt(c2)
c1$value<-c1$value/226*100
c2$value<-c2$value/69*100

subset(c2,group=="III")->tmp2
c2$variable<-factor(c2$variable,levels=tmp2[order(tmp2$value,decreasing=T),"variable"])
c1$variable<-factor(c1$variable,levels=tmp2[order(tmp2$value,decreasing=T),"variable"])

pdf("cluster3.barplot.pdf",height=4,width=4.5)
ggplot(c1)+geom_bar(aes(x=variable,y=value,fill=group),stat="identity")+
	scale_fill_manual(values=brewer.pal(n=11,name="RdYlBu")[c(1,2,11,10)])+
	theme_pubr()+
	theme(axis.text.x=element_text(angle=45,hjust=1,size=8), plot.margin=unit(c(1,1,1,15),"mm"))+
	labs(x=NULL,y="Percentage (%)",title="3.1")
ggplot(c2)+geom_bar(aes(x=variable,y=value,fill=group),stat="identity")+
	scale_fill_manual(values=brewer.pal(n=11,name="RdYlBu")[c(1,2,11,10)])+
	theme_pubr()+
	theme(axis.text.x=element_text(angle=45,hjust=1,size=8), plot.margin=unit(c(1,1,1,15),"mm"))+
	labs(x=NULL,y="Percentage (%)",title="3.2")
dev.off()

#################
rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
FetchData(merged.hg38,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged.hg38$seurat_clusters,rlv)
FetchData(merged.mm10,pw)->ratio_df2
ratio_df2$seurat_clusters<-factor(merged.mm10$seurat_clusters,rlv)
coeff<-sapply(1:24,function(x)max(ratio_df2[,x])/max(ratio_df[,x]))
ratio_df2[,1:24]<-sweep(ratio_df2[,1:24],2,coeff,"/")

ratio_df<-reshape2::melt(ratio_df)
ratio_df$group<-"hg38"
ratio_df2<-reshape2::melt(ratio_df2)
ratio_df2$group<-"mm10"
rbind(ratio_df,ratio_df2)->ratio_df3
ratio_df3$group<-factor(ratio_df3$group,levels=c("hg38","mm10"))
ratio_df3$xaxis<-paste(ratio_df3$group,ratio_df3$seurat_clusters,sep=":")
ratio_df3$xaxis<-factor(ratio_df3$xaxis,paste(rep(c("hg38","mm10"),each=11),rlv,sep=":"))

pdf("cluster3.pathways.pdf",width=3.5,height=2)
lapply(1:24,function(x){
	df<-subset(ratio_df3,seurat_clusters%in%c(3.1,3.2) & value>0 & variable==paste0("hallmark_",pw[x]))
	g<-ggplot(df,aes(factor(xaxis),value,color=group))+
		scale_x_discrete(labels=rep(c(3.1,3.2),2))+
		scale_y_continuous(name="hg38",sec.axis=sec_axis(~.*coeff[x], name="mm10"))+
		geom_violin(scale="area",size=0.2)+
		#geom_boxplot(show.legend=F,outlier.shape=NA,width=0.2)+
		geom_beeswarm(priority="random",method="frowney",
			cex=1,dodge.width=0,size=0.1,show.legend=F,alpha=0.2
		)+
		#stat_compare_means(aes(label = signif(..p.signif..,3)),comparisons=list(c(1,2),c(3,4)),paired=F,method="t.test",tip.length=0.01,vjust = 2)+
		theme_pubr(border=T)+
		theme(legend.position="right",axis.text.y.left= element_text(color="red"),axis.text.y.right= element_text(color="blue"),
			axis.ticks.y.left= element_line(color="red"),axis.ticks.y.right= element_line(color="blue"))+
		scale_color_manual(breaks=c("mm10","hg38"),values=c("blue","red"),labels=c("mm10","hg38"))+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		guides(y = guide_axis(position="left"), y.sec = guide_axis(position="right"))+
		labs(x=NULL,y="AUC score",color=NULL,title=pw[x]) + 
		coord_cartesian(expand=T)
	#g$layers[[3]]$aes_params$textsize <- 2
	print(g)
	NULL
})
dev.off()

######
#pie
p1<-table(subset(merged[[]],Sample=="MTF")$old)
p2<-table(subset(merged[[]],Sample=="PAR")$old)

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



###################
custom1<-parse_geneset_term("custom.gs3.txt") #human
custom2<-parse_geneset_term("custom.gs4.txt") #mouse
names(custom2)<-paste0("M-",names(custom2))
parse_geneset_term("~/work/MSigDB/mouse/h.all.v6.1.symbols_mouse.gmt")->hallmark.mm10
names(hallmark.mm10)<-paste0("M-",names(Hallmark))

AUCell_buildRankings(merged@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
AUCell_calcAUC(c(Hallmark,hallmark.mm10,custom1,custom2),cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
merged[["Hallmark"]]<-CreateAssayObject(data=cells_AUC@assays@data$AUC)

rlv<-c(1,5,6,8,0,2,9,4,7,3.1,3.2)
pw<-names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)]
FetchData(merged,pw)->ratio_df
ratio_df$seurat_clusters<-factor(merged$seurat_clusters,rlv)
FetchData(merged,paste0("M-",pw))->ratio_df2
ratio_df2$seurat_clusters<-factor(merged$seurat_clusters,rlv)
colnames(ratio_df)<-c(paste0("H_",pw),"seurat_clusters")
colnames(ratio_df2)<-c(paste0("M_",pw),"seurat_clusters")
ratio_df3<-cbind(ratio_df,ratio_df2)

pdf("cluster3.IFNg.pdf",height=3,width=3.5)
lapply(pw,function(x){
	g<-ggplot(subset(ratio_df3,seurat_clusters %in% c(3.1,3.2)),aes_string(as.name(paste0("H_",x)),as.name(paste0("M_",x))))+
		geom_jitter(aes(color=seurat_clusters),shape=21,show.legend=F)+theme_pubr(border=T)+
		labs(title="3.1 and 3.2")+
		scale_color_manual(values=c("green","purple"))#+
		#coord_cartesian(xlim=c(0,0.05),ylim=c(0,0.06))
	print(g)
	NULL
})
dev.off()

##################
#3.2, 3.2, 1.1, 1.2, 8.1, 8.2
merged<-StashIdent(merged,save.name="old")
as.character(merged$seurat_clusters)->clusters
clusters[which(merged$seurat_clusters==3&merged$Sample=="MTF")]<-"3.2"
clusters[which(merged$seurat_clusters==3&merged$Sample=="PAR")]<-"3.1"
clusters[which(merged$seurat_clusters==1&merged$Sample=="MTF")]<-"1.2"
clusters[which(merged$seurat_clusters==1&merged$Sample=="PAR")]<-"1.1"
clusters[which(merged$seurat_clusters==8&merged$Sample=="MTF")]<-"8.2"
clusters[which(merged$seurat_clusters==8&merged$Sample=="PAR")]<-"8.1"

rlv<-c(1.1,5,6,8.1,0,2,9,3.1,3.2,8.2,1.2,4,7)
factor(clusters,levels=rlv)->clusters
Idents(merged)<-clusters
merged$seurat_clusters<-clusters
names(clusters)<-colnames(merged)

MTF[[]]->test2
PAR[[]]->test3
merged[[]]->test4
test2$seurat_clusters<- -1
test3$seurat_clusters<- -1
test2[str_replace(rownames(subset(test4,seurat_clusters==1.2)),"_1",""),"seurat_clusters"]<-1.2
test2[str_replace(rownames(subset(test4,seurat_clusters==3.2)),"_1",""),"seurat_clusters"]<-3.2
test2[str_replace(rownames(subset(test4,seurat_clusters==4)),"_1",""),"seurat_clusters"]<-4
test2[str_replace(rownames(subset(test4,seurat_clusters==7)),"_1",""),"seurat_clusters"]<-7
test2[str_replace(rownames(subset(test4,seurat_clusters==8.2)),"_1",""),"seurat_clusters"]<-8.2

test3[str_replace(rownames(subset(test4,seurat_clusters==1.1)),"_2",""),"seurat_clusters"]<-1.1
test3[str_replace(rownames(subset(test4,seurat_clusters==2)),"_2",""),"seurat_clusters"]<-2
test3[str_replace(rownames(subset(test4,seurat_clusters==3.1)),"_2",""),"seurat_clusters"]<-3.1
test3[str_replace(rownames(subset(test4,seurat_clusters==5)),"_2",""),"seurat_clusters"]<-5
test3[str_replace(rownames(subset(test4,seurat_clusters==6)),"_2",""),"seurat_clusters"]<-6
test3[str_replace(rownames(subset(test4,seurat_clusters==8.1)),"_2",""),"seurat_clusters"]<-8.1
test3[str_replace(rownames(subset(test4,seurat_clusters==9)),"_2",""),"seurat_clusters"]<-9
test3[str_replace(rownames(subset(test4,seurat_clusters==0)),"_2",""),"seurat_clusters"]<-0

#test2$seurat_clusters<-test2$seurat_clusters==3.1
#test3$seurat_clusters<-test3$seurat_clusters==3.2
pdf("cluster3.v2.pdf",width=3.5,height=3)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==1.1),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster1.1")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==2),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster2")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==3.1),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster3.1")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==5),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster5")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==6),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster6")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==8.1),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster8.1")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==9),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster9")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test3,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==0),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster0")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
	
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==1.2),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster1.2")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==3.2),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster3.2")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==4),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster4")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==7),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster7")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
ggplot(test2,aes(x=row,y=col))+
	geom_point(aes(fill=seurat_clusters==8.2),shape=21,stroke=0.1,size=1,color="grey80")+
	ggpubr::theme_pubr(border=T)+coord_cartesian(ylim=c(77,0),xlim=c(0,127))+labs(x=NULL,y=NULL,alpha=NULL,title="cluster8.2")+
	scale_fill_manual(values=c("white","purple"),guide=NULL)
dev.off()


#

rlv<-c(1.1,5,6,8.1,0,2,9,3.1,3.2,8.2,1.2,4,7)
rlv<-c(6,8.1,1.1,5,3.1,9,2,0,1.2,8.2,3.2,4,7)
factor(clusters,levels=rlv)->clusters
Idents(merged)<-clusters
merged$seurat_clusters<-clusters
names(clusters)<-colnames(merged)

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

pdf("Ratio.beeswarm.v3.pdf",width=4.8,height=2)
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

#
custom1<-parse_geneset_term("custom.gs3.txt") #human
custom2<-parse_geneset_term("custom.gs4.txt") #mouse

merged[rownames(merged) %in% hg38,]->merged.hg38
AUCell_buildRankings(merged.hg38@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings
AUCell_calcAUC(c(Hallmark,custom1),cells_rankings,aucMaxRank=500,nCores=32)->cells_AUC
merged.hg38[["Hallmark"]]<-CreateAssayObject(data=cells_AUC@assays@data$AUC)

parse_geneset_term("~/work/MSigDB/mouse/h.all.v6.1.symbols_mouse.gmt")->hallmark.mm10
names(hallmark.mm10)<-names(Hallmark)
merged[rownames(merged) %in% mm10,]->merged.mm10
AUCell_buildRankings(merged.mm10@assays$SCT@data,plotStats=F,nCores=32)->cells_rankings2
AUCell_calcAUC(c(hallmark.mm10,custom2),cells_rankings2,aucMaxRank=500,nCores=32)->cells_AUC2 
merged.mm10[["Hallmark"]]<-CreateAssayObject(data=cells_AUC2@assays@data$AUC)


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

########################
DefaultAssay(merged.hg38)<-"Hallmark"
merged.hg38[lvs[1:8],]->hg38.ss1
merged.hg38[lvs[9:24],]->hg38.ss2
merged.hg38[pw,]->hg38.all

hg38.all <- hg38.all %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs)
hg38.ss1 <- hg38.ss1 %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs[1:8])
hg38.ss2 <- hg38.ss2 %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs[9:24])

DefaultAssay(merged.mm10)<-"Hallmark"
merged.mm10[lvs[1:8],]->mm10.ss1
merged.mm10[lvs[9:24],]->mm10.ss2
merged.mm10[pw,]->mm10.all

mm10.all <- mm10.all %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs)
mm10.ss1 <- mm10.ss1 %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs[1:8])
mm10.ss2 <- mm10.ss2 %>% 
	ScaleData() %>%
	RunPCA(assay = "Hallmark", verbose = FALSE,feature=lvs[9:24])
	
pdf("Distance.pathways.pdf",width=5,height=4)
DimPlot(hg38.all, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: All pathways")
DimPlot(hg38.ss2, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: Category I")
DimPlot(hg38.ss1, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: Category II")

DimPlot(mm10.all, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: All pathways")
DimPlot(mm10.ss2, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: Category I")
DimPlot(mm10.ss1, reduction = "pca", label = TRUE,group.by="seurat_clusters")+labs(title="Human: Category II")
dev.off()

#
color<-scales:::hue_pal()(10)
merged.hg38@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
mat2<-mat2[pw[c(14,15,16,13,18,10,9,12,4,19,17,21,22,24,20,11,2,1,3,5,7,6,23,8)],]
#o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.human.v3.pdf",height=3,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,cluster_rows=F,
	clustering_distance_rows="pearson",
	#row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		labels_gp = gpar(fontsize = 1)
		),
		height=unit(3,"mm")
	)
)
Heatmap(Hallmark.mat[pw[c(14,15,16,13,18,10,9,12,4,19,17,21,22,24,20,11,2,1,3,5,7,6,23,8)],],name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,0.15,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,cluster_rows=F,
	clustering_distance_rows="pearson",
	#row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		labels_gp = gpar(fontsize = 1)
		),
		height=unit(3,"mm")
	)
)
dev.off()

merged.mm10@assays$Hallmark@data[,names(clusters)]->Hallmark.mat
mat2<-(sweep(Hallmark.mat,1,apply(Hallmark.mat,1,function(x)quantile(x[x!=0],0.75)),"/"))
pw<-c("M1 marker","M2 marker",names(Hallmark)[c(40,10,19,18,23,31,2,35,33,30,27,28,26,37,48,25,6,15,43,7,1,49)])
mat2<-mat2[pw[c(14,15,16,13,18,10,9,12,4,19,17,21,22,24,20,11,2,1,3,5,7,6,23,8)],]
#o = seriate(dist(mat2), method = "OLO_complete")
pdf("SC.MSigDB.Hallmark.mouse.v3.pdf",height=3,width=5)
Heatmap(mat2,name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,2,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,cluster_rows=F,
	clustering_distance_rows="pearson",
	#row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		labels_gp = gpar(fontsize = 1)
		),
		height=unit(3,"mm")
	)
)
Heatmap(Hallmark.mat[pw[c(14,15,16,13,18,10,9,12,4,19,17,21,22,24,20,11,2,1,3,5,7,6,23,8)],],name="AUCscore",
	show_column_names=F,row_names_gp=gpar(cex=0.5),
	col=circlize::colorRamp2(seq(0,0.15,length=50),colors=Seurat:::CustomPalette(low = "green", high = "red", mid = "black",k=50)),
	cluster_columns=F,cluster_rows=F,
	clustering_distance_rows="pearson",
	#row_order=get_order(o,1),
	column_split = clusters,column_gap = unit(0.35, "mm"),
	use_raster=T, raster_device="png",raster_quality=8,
	top_annotation=HeatmapAnnotation(foo = anno_block(gp = gpar(fill = color[rlv+1],col=NA),
		#labels = levels(clusters), 
		labels_gp = gpar(fontsize = 1)
		),
		height=unit(3,"mm")
	)
)
dev.off()


###
as.data.frame(t(merged.hg38@assays$Hallmark@data[,names(clusters)]))->tmp1
tmp1$Cluster<-clusters
as.data.frame(t(merged.mm10@assays$Hallmark@data[,names(clusters)]))->tmp2
tmp2$Cluster<-clusters
pdf("Myc.plot.pdf",width=3.5,height=2.3)
ggplot(tmp1,aes_string("Cluster",as.name("Myc Targets V1"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="Myc Targets in Human") + 
	coord_cartesian(expand=T)
ggplot(tmp1,aes_string("Cluster",as.name("Androgen Response"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="Androgen Response in Human") + 
	coord_cartesian(expand=T)
ggplot(tmp1,aes_string("Cluster",as.name("Fatty Acid Metabolism"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="Fatty Acid Metabolism in Human") + 
	coord_cartesian(expand=T)

dev.off()

pdf("M1_M2.plot.pdf",width=3.5,height=2.3)
ggplot(tmp1[tmp1[,"M1 marker"]>0,],aes_string("Cluster",as.name("M1 marker"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="M1 marker in Human") + 
	coord_cartesian(expand=T)
ggplot(tmp1[tmp1[,"M2 marker"]>0,],aes_string("Cluster",as.name("M2 marker"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="M2 marker in Human") + 
	coord_cartesian(expand=T)
ggplot(tmp2[tmp2[,"M1 marker"]>0,],aes_string("Cluster",as.name("M1 marker"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="M1 marker in Mouse") + 
	coord_cartesian(expand=T)
ggplot(tmp2[tmp2[,"M2 marker"]>0,],aes_string("Cluster",as.name("M2 marker"),color="Cluster"))+
	geom_violin(scale="area",size=0.2,show.legend=F,color="red")+
	geom_beeswarm(priority="random",method="frowney",
		cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
	)+
	theme_pubr(border=T)+
	#scale_color_manual(values="red")+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x=NULL,y="AUC score",color=NULL,title="M2 marker in Mouse") + 
	coord_cartesian(expand=T)
dev.off()


#####
genes<-c("TJP1", "TJP2", "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "OCLN", 
	"VIM", "COL1A1", "FN1", "TGFB1", "SNAI1", "MMP2", "MMP9", "TWIST2", "ZEB1", "ITGAV", "VCL", "PXN")
FetchData(merged.hg38,genes,cells=names(clusters))->tmp1
tmp1$Cluster<-clusters

pdf("EMT.Marker_genes.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp1,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="red",trim=T)+
		# geom_beeswarm(priority="random",method="frowney",
			# cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
		# )+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()
pdf("EMT.Marker_genes.v2.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp1,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="red",trim=T)+
		geom_beeswarm(priority="density",method="frowney",
			cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
		)+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()


###

#####
genes<-c("VIM","CDH1","Epithelial Mesenchymal Transition")
FetchData(merged.hg38,genes,cells=names(clusters))->tmp1
colnames(tmp1)<-genes
tmp1$Cluster<-clusters

pdf("beeswarm.EMT.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp1,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="red",trim=T)+
		# geom_beeswarm(priority="random",method="frowney",
			# cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
		# )+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()
pdf("beeswarm.EMT.v2.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp1,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="red",trim=T)+
		geom_beeswarm(priority="density",method="frowney",
			cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
		)+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()

#
genes<-c("Vim","Cdh1","Epithelial Mesenchymal Transition")
FetchData(merged.mm10,genes,cells=names(clusters))->tmp2
colnames(tmp2)<-genes
tmp2$Cluster<-clusters

pdf("beeswarm.mouse.EMT.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp2,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="blue",trim=T)+
		# geom_beeswarm(priority="random",method="frowney",
			# cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="red"
		# )+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()
pdf("beeswarm.mouse.EMT.v2.pdf",width=3.5,height=2.3)
lapply(genes,function(x){
	g<-ggplot(tmp2,aes_string("Cluster",as.name(x),color="Cluster"))+
		geom_violin(scale="width",size=0.2,show.legend=F,color="blue",trim=T)+
		geom_beeswarm(priority="density",method="frowney",
			cex=0.4,size=0.1,show.legend=F,alpha=0.2,color="blue"
		)+
		theme_pubr(border=T)+
		#scale_color_manual(values="red")+
		guides(color=guide_legend(override.aes = list(size = 5)))+
		labs(x=NULL,y="Expression level",color=NULL,title=x) + 
		coord_cartesian(expand=T)
	print(g)
	NULL
})
dev.off()


###
genes<-c("MYC","M1 marker","M2 marker")
FetchData(merged.hg38,"MYC",cells=names(clusters),slot="scale.data")->tmp1
FetchData(merged.hg38,c("M1 marker","M2 marker"),cells=names(clusters))->tmp2
cbind(tmp1,tmp2)->tmp1
colnames(tmp1)<-genes
tmp1$Cluster<-clusters
tmp1$Sample<-merged.hg38$Sample

pdf("Cor.MYC.M1_M2.pdf",width=2.6,height=2.5)
ggplot(subset(tmp1,Sample=="PAR"),aes_string("MYC",as.name("M1 marker")))+
	geom_point(shape=21,size=0.7)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="MYC",y="M1 marker",color=NULL,title="Par\ncor=-0.01,pv=0.67")
ggplot(subset(tmp1,Sample=="PAR"),aes_string("MYC",as.name("M2 marker")))+
	geom_point(shape=21,size=0.7)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="MYC",y="M2 marker",color=NULL,title="Par\ncor=-0.0074,pv=0.75")
ggplot(subset(tmp1,Sample=="MTF"),aes_string("MYC",as.name("M1 marker")))+
	geom_point(shape=21,size=0.7)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="MYC",y="M1 marker",color=NULL,title="TMH\ncor=0.20,pv=1.03e-7")
ggplot(subset(tmp1,Sample=="MTF"),aes_string("MYC",as.name("M2 marker")))+
	geom_point(shape=21,size=0.7)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="MYC",y="M2 marker",color=NULL,title="TMH\ncor=0.16,pv=4.94e-5")
dev.off()


###
SCTransform(merged,assay = "Spatial", return.only.var.genes=F,verbose = FALSE)->merged2
merged2[rownames(merged2) %in% hg38,]->merged.hg38.2
FetchData(merged.hg38.2,"CXCR4",slot="scale.data")

genes<-c("CD44", "CXCR4", "VIM", "FN1", "KRT19", "EPCAM","M1 marker","M2 marker")
FetchData(merged.hg38,c("CD44", "CXCR4", "VIM"),cells=names(clusters),slot="scale.data")->tmp1
FetchData(merged.hg38,c("FN1", "KRT19", "EPCAM"),cells=names(clusters))->tmp3

FetchData(merged.hg38,c("M1 marker","M2 marker"),cells=names(clusters))->tmp2
cbind(cbind(tmp1,tmp3),tmp2)->tmp1
colnames(tmp1)<-genes
tmp1$Cluster<-clusters
tmp1$Sample<-merged.hg38$Sample

pdf("Cor.test.pdf",width=3,height=3.5)
lapply(c("CD44", "CXCR4", "VIM", "FN1", "KRT19", "EPCAM"),function(x){
limits<-quantile(tmp1[,x],c(0.05,0.95))
g1<-ggplot(subset(tmp1,Sample=="MTF"),aes_string(as.name("M1 marker"),as.name("M2 marker")))+
	geom_point(aes_string(fill=x),shape=21,size=1)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="M1 marker",y="M2 marker",color=NULL,title="MTF")+
	scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
g2<-ggplot(subset(tmp1,Sample=="PAR"),aes_string(as.name("M1 marker"),as.name("M2 marker")))+
	geom_point(aes_string(fill=x),shape=21,size=1)+
	theme_pubr(border=T)+
	guides(color=guide_legend(override.aes = list(size = 5)))+
	labs(x="M1 marker",y="M2 marker",color=NULL,title="PAR")+
	scale_fill_gradientn(limits=limits,
			colours=rev(RColorBrewer::brewer.pal(11,"Spectral")),oob = scales::squish)
print(g1)
print(g2)
NULL
}
)
dev.off()

