#DEseq2 for GBR common garden RNAseq - Orpheus and Wilke 

#installation of DESeq2
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
source("http://bioconductor.org/biocLite.R") 
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library(DESeq2)
library(arrayQualityMetrics)
library(VennDiagram)
library(genefilter)
library(gplots)
library(vegan)
library(rgl)
library(ape)
library(WGCNA)

# exporting data for WGCNA 

degs.select=which(wcwd$pvalue<0.022)
degs.select2=which(wco$pvalue<0.021)
degs.select3=which(wdo$pvalue<0.021)
degs.sel=c(degs.select, degs.select2, degs.select3)
degs.unique=unique(degs.sel)
length(degs.unique) # 6138
degs=vsd[degs.unique,]
save(vsd,degs,samples,file="wgcna_input.RData")

o0=c(rep(1,3), rep(0, 54))
o10=c(rep(0,3), rep(1,3), rep(0,51))
o1=c(rep(0,6), rep(1,3), rep(0,48))
o2=c(rep(0,9), rep(1,3), rep(0, 45))
o3=c(rep(0,12), rep(1,3), rep(0,42))
o4=c(rep(0,15), rep(1,3), rep(0,39))
o6=c(rep(0,18), rep(1,3), rep(0,36))
o7=c(rep(0,21), rep(1,3), rep(0,33))
o8=c(rep(0,24), rep(1,2), rep(0,31))
o9=c(rep(0,26), rep(1,3), rep(0,28))
oM1=c(rep(0,29), rep(1,3), rep(0,25))
w22=c(rep(0,32), rep(1,2), rep(0,23))
w23=c(rep(0,34), rep(1,3), rep(0,20))
w24=c(rep(0,37), rep(1,3), rep(0,17))
w25=c(rep(0,40), rep(1,2), rep(0,15))
w26=c(rep(0,42), rep(1,3), rep(0,12))
w27=c(rep(0,45), rep(1,3), rep(0,9))
w28=c(rep(0,48), rep(1,3), rep(0,6))
w30=c(rep(0,51), rep(1,3), rep(0,3))
w31=c(rep(0,54), rep(1,3))

Wc=c(rep(0,32), rep(1,13), rep(0,12))
Wd=c(rep(0,45), rep(1,12))
O=c(rep(1,32), rep(0,25))

group=c(rep("Orph", 32), rep("WilkieC", 13), rep("WilkieD", 12))
ind=c(rep("o0", 3), rep("o10",3), rep("o1",3), rep("o2",3), rep("o3",3), rep("o4",3), rep("o6",3), rep("o7",3), rep("o8",2), rep("o9",3), rep("oM1",3), rep("w22",2), rep("w23",3), rep("w24",3), rep("w25",2), rep("w26",3), rep("w27",3), rep("w28",3), rep("w30",3), rep("w31",3))

samples=cbind(group, ind, Wc, Wd, O, o0, o10, o1, o2, o3, o4, o6, o7, o8, o9, oM1, w22, w23, w24, w25, w26, w27, w28, w30, w31)
samples=data.frame(samples)

#############################################################################

lnames=load(file = "wgcna_input.RData")
lnames
datt=t(degs)
samples

library(WGCNA)
options(stringsAsFactors=FALSE)


# Choose a set of soft-thresholding powers

powers = c(c(1:10), seq(from = 12, to=26, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datt, powerVector = powers, verbose = 5,networkType="signed")
# Plot the results:
#sizeGrWindow(9, 5)
#pdf("soft_threshold_signed.pdf",height=4, width=8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()


# making modules
library(flashClust)
s.th=12 # re-specify according to previous section
adjacency = adjacency(datt, power = s.th,type="signed");
TOM = TOMsimilarity(adjacency,TOMType="signed");
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

save(dissTOM,geneTree,file="signedDissTOM_tree.RData")

# We like large modules, so we set the minimum module size relatively high:
#usually kept at 30
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

# Calculate eigengenes
MEList = moduleEigengenes(datt, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = flashClust(as.dist(MEDiss), method = "average");

save(dynamicMods,dynamicColors,MEs,METree,geneTree,file="1stPassModules.RData")

load(file="1stPassModules.RData")

sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
MEDissThres = 0.5
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datt, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

quartz()
# plotting the fabulous ridiculogram
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")

# how many genes in each module?
table(moduleColors)

#moduleColors
#black           blue          brown           cyan       darkgrey    darkmagenta darkolivegreen     darkorange 
#662            406            380            178            103             41            160             98 
#darkred          green         grey60      lightcyan     lightgreen    lightyellow        magenta   midnightblue 
#116            806            141            164            792            132            345            245 
#purple      steelblue          white         yellow 
#857             83             96            333 


# Save module colors and labels for use in subsequent parts
save(MEs, geneTree, moduleLabels, moduleColors, file = "networkdata_signed.RData")


###################
# plotting correlations with traits:
load(file = "networkdata_signed.RData")
load(file = "wgcna_input.RData");
datt=t(degs)

# Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, samp2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(samp2),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-1,1),
main = paste("Module-trait relationships"))


# scatterplots of gene significance (correlation-based) vs kME

load(file = "networkdata_signed.RData")
load(file = "wgcna_input.RData");
samp2
table(moduleColors)
whichTrait="Wc"

datt=t(degs)
nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(samples[,whichTrait]);
names(selTrait) = whichTrait
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
par(mfrow=c(3,3))
counter=0
for(module in modNames[1:length(modNames)]){
counter=counter+1
if (counter>9) {
	quartz()
	par(mfrow=c(3,3))
	counter=1
}
column = match(module, modNames);
moduleGenes = moduleColors==module;
#trr="heat resistance"
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, 1]),
xlab = paste(module,"module membership"),
ylab = paste("GS for", whichTrait),
col = "grey",mgp=c(2.3,1,0))
}


################
# eigengene-heatmap plot

load(file = "networkdata_signed.RData")
load(file = "wgcna_input.RData");
datt=t(degs)
samples

which.module="blue"
datME=MEs
datExpr=datt
quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,rlabels=F,rcols=which.module,
main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
ylab="eigengene expression",xlab="sample")

length(datExpr[1,moduleColors==which.module ]) 

#################
# saving selected modules for Fisher-based GO analysis

load(file = "networkdata_signed.RData")
load(file = "wgcna_input.RData");

whichModule="blue"
allkME =as.data.frame(signedKME(t(degs), MEs))

table(moduleColors==whichModule)
genes=row.names(degs)[moduleColors==whichModule]
#write.table(genes, file="genes")

fishers=data.frame("gene"=row.names(degs),"inModule"=as.numeric(row.names(degs) %in% genes))
nonModGenes=fishers[!fishers$gene %in% genes,]
blue_k=allkME[row.names(allkME) %in% genes,]
Mod_isos=row.names(blue_k)
#change the numbers to select which module you want to keep 
blue_only=blue_k[,c(-1,-2,-3,-4,-5, -6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16, -18, -19, -20)]
Mod_isos=data.frame(Mod_isos)
blue_only=data.frame(blue_only)
blue2=cbind(Mod_isos, blue_only)
colnames(blue2)=c("gene", "inModule")
#combine non-module genes and module genes (with kME values)
blue_df=rbind(blue2, nonModGenes)

#table(fishers$inModule)
write.csv(blue_df,file=paste(whichModule,".csv",sep=""),row.names=F,quote=F)


################
# plotting heatmap for named top-kME genes

load(file = "networkdata_signed.RData")
load(file = "wgcna_input.RData");
allkME =as.data.frame(signedKME(t(degs), MEs))
gg=read.table("amil_iso2gene.tab",sep="\t")
library(pheatmap)

###########################################################################

whichModule="blue"
top=164
quartz()
datME=MEs
datExpr=t(degs)
modcol=paste("kME",whichModule,sep="")
sorted=degs[order(allkME[,modcol],decreasing=T),]
head(sorted)
# selection top N names genes, attaching gene names
gnames=c();counts=0;hubs=c()
for(i in 1:length(sorted[,1])) {
	if (row.names(sorted)[i] %in% gg$V1) { 
		counts=counts+1
		gn=gg[gg$V1==row.names(sorted)[i],2]
		gn=paste(gn,row.names(sorted)[i],sep=".")
		if (gn %in% gnames) {
			gn=paste(gn,counts,sep=".")
		}
		gnames=append(gnames,gn) 
		hubs=data.frame(rbind(hubs,sorted[i,]))
		if (counts==top) {break}
	}
} 
row.names(hubs)=gnames

write.table(gnames, file="blue_module_genes")

#########################################################################

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

pheatmap(hubs,scale="row",col=contrasting,border_color=NA,treeheight_col=0,cluster_rows=T, cluster_cols=F)

##############################################################################









