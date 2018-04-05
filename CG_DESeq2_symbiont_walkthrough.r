#DESeq2 GBR common garden samples Symbiont gene expression 
source("http://bioconductor.org/biocLite.R")
install.packages("BiocInstaller",repos="http://www.bioconductor.org/packages/3.1/bioc")
biocLite("DESeq2")
biocLite("DESeq")
library(DESeq2)
library(arrayQualityMetrics)
library(VennDiagram)
library(genefilter)
library(gplots)
library(vegan)
library(rgl)
library(ape)
library(empiricalFDR.DESeq2)
library(ggplot2)
library(dendextend)
library(whisker)

#read counts
reads=read.table("allcounts_symC.txt", header=T, row.names=1)

names(reads)= sub(".sam.counts", "", names(reads))
ncol(reads)
nrow(reads)

#too few reads mapped for these samples and remove D individuals
reads$o8a=NULL
reads$o8b=NULL
reads$o8c=NULL
reads$oM1a=NULL
reads$oM1b=NULL
reads$oM1c=NULL
reads$o4a=NULL
reads$o4b=NULL
reads$o4c=NULL

reads$w27a=NULL
reads$w27b=NULL
reads$w27c=NULL
reads$w28a=NULL
reads$w28b=NULL
reads$w28c=NULL
reads$w30a=NULL
reads$w30b=NULL
reads$w30c=NULL
reads$w31a=NULL
reads$w31b=NULL
reads$w31c=NULL

#sum replicates 

yur=cbind(o0=rowSums(reads[,c('o0a','o0b','o0c')]),o10=rowSums(reads[,c('o10a','o10b','o10c')]), o1=rowSums(reads[,c('o1a', 'o1b', 'o1c')]), o2=rowSums(reads[,c('o2a', 'o2b', 'o2c')]), o3=rowSums(reads[,c('o3a', 'o3b', 'o3c')]), o6=rowSums(reads[,c('o6a', 'o6b', 'o6c')]), o7=rowSums(reads[,c('o7a', 'o7b', 'o7c')]), o9=rowSums(reads[,c('o9a', 'o9b', 'o9c')]), w22=rowSums(reads[,c('w22a', 'w22b', 'w22c')]), w23=rowSums(reads[,c('w23a', 'w23b', 'w23c')]), w24=rowSums(reads[,c('w24a', 'w24b', 'w24c')]), w25=rowSums(reads[,c('w25a', 'w25b', 'w25c')]), w26=rowSums(reads[,c('w26a', 'w26b', 'w26c')]))

# how many genes have mean count>2?
means=apply(yur,1,mean)
table(means>2)

# removing all genes with mean count less than 2 
yur=yur[means>2,]
nrow(yur)
#4642 genes kept 

yur=data.frame(yur)
ind=colnames(yur)
location=c(rep("Orpheus", 8), rep("Wilkie", 5))

samples=cbind(ind, location)
samples=data.frame(samples)

#model - location
dds_sym=DESeqDataSetFromMatrix(colData=samples, countData=yur, design=~location)

### make a PCA plot ###

vsd_PCA=varianceStabilizingTransformation(dds_sym)

quartz()

plotPCA = function (x, intgroup = "location", ntop = 500) 
{
  rv = rowVars(assay(x))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select, ]))
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]), 
                     1, paste, collapse = " : "))
  if (nlevels(fac) > 3) 
    #colours = brewer.pal(nlevels(fac), "Paired")
    colours=c("cornflowerblue", "azure4", "chocolate2", "firebrick", "darkolivegreen2", "deepskyblue2", "dimgrey", "goldenrod", "darkorchid1", 
              "gray1", "mediumpurple1", "magenta", "mistyrose", "midnightblue", "orange", "seagreen1", "royalblue", "yellow", "tan3", "turquoise2")
  else colours = c("red", "dodgerblue", "pink")
  xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
         pch = 2, cex = 1, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
                                                                                     text = list(levels(fac)), rep = FALSE)))
}

plotPCA(vsd_PCA, intgroup="location", ntop=500)


#get final deseq2 results
final_sym = DESeq(dds_sym)
vsd_sym=getVarianceStabilizedData(final_sym)
final_results_sym=results(final_sym)

### simulate counts for fdr ####
sim.sym=simulateCounts(final_sym)
simSym=DESeq(sim.sym)
simSym_results=results(simSym)

save(simSym_results, final_results_sym, samples, vsd_sym, file="DESeq2_symC.RData")

library(DESeq2)
library(empiricalFDR.DESeq2)
library(pheatmap)
load("DESeq2_symC.RData")

#### determine 10% fdr cutoff ##### 
table(final_results_sym$padj<0.1) # filtering based DEGs
table(is.na(final_results_sym$padj)) # this is how many needed to be tossed
fdrt=fdrTable(final_results_sym$pvalue,simSym_results$pvalue)
par(mfrow=c(1,2))
fdrBiCurve(fdrt,main="Site of Origin")
quartz()
efdrS=empiricalFDR(fdrt,plot=T,main="Site of Origin")
mtext(paste("10% FDR at p =",signif(efdrS,2)),cex=0.8)
table(final_results_sym$pval<efdrS) # empirical FDR based DEGs
final_results_sym$efdr=(final_results_sym$pval<efdrS)
final_results_sym$efdr[is.na(final_results_sym$efdr)]=FALSE
#FDR10% with no outliers is at p=0.0013

quartz()
degs_efdr_sym=which(final_results_sym$efdr)
length(degs_efdr_sym) #52
pheatmap(cor(vsd_sym[degs_efdr_sym,]),border_color=NA) 
pheatmap(cor(vsd_sym[vsd_sym,]), border_color=NA)

#principal coordinates 
dev.off()
ad.pcoa=pcoa(dist(t(vsd_sym),method="manhattan")/1000)
scores=ad.pcoa$vectors

#by location
quartz()
plot(scores[,1], scores[,2],pch=as.numeric(as.factor(samples$ind)),col=as.numeric(as.factor(samples$location)))
ordiellipse(scores, samples$location, label=T)
ordispider(scores, samples$ind)

#by genotype 
quartz()
plot(scores[,1], scores[,2],pch=as.numeric(as.factor(samples$ind)),col=as.numeric(as.factor(samples$ind)))
ordispider(scores, samples$ind,label=T)

#analysis of variance
adonis(t(vsd_sym)~samples$location,data=samples,method="manhattan")

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#samples$location  1   1119486 1119486  1.4985 0.11989  0.013 *
#  Residuals        11   8217954  747087         0.88011         
#Total            12   9337440                 1.00000         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



### sort top DEGs 

sym_fdr=subset(final_results_sym, pvalue < 0.0013)
sym_FDR=sym_fdr[order(sym_fdr$pvalue), , drop=FALSE]
sym_FDR_sorted=sym_FDR[2]
degs=row.names(sym_FDR_sorted)
degs_sym=row.names(sym_FDR_sorted)
sorted=vsd_sym[degs,]
write.table(degs, file="degs_sym.txt")
write.table(sym_FDR, file="sym_degs.txt", sep="\t")


gg=read.table("contig2gene.tab", sep="\t")
library(pheatmap)
library(gplots)
library(RColorBrewer)

top=52

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

contrasting = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting2 = colorRampPalette(rev(c("chocolate1","chocolate1","#FEE090","grey10", "cyan3","cyan")))(100)
contrasting3 = colorRampPalette(rev(c("chocolate1","#FEE090","grey10", "cyan3","cyan","cyan")))(100)

quartz()
pheatmap(hubs,scale="row",col=contrasting,border_color=NA,treeheight_col=0, cluster_rows=T, cluster_cols=F)

### save data for GO analysis 

#signed -log(p-values): origin
sym.p=data.frame("gene"=row.names(final_results_sym[!is.na(final_results_sym$pvalue),]))
sym.p$lpv=-log(final_results_sym[!is.na(final_results_sym$pvalue),"pvalue"]+1e-15,10)
sym.fc=data.frame("gene"=row.names(final_results_sym[!is.na(final_results_sym$pvalue),]))
sym.fc$lfc=final_results_sym[!is.na(final_results_sym$pvalue),"log2FoldChange"]
direction=as.numeric(sym.fc[,"lfc"]>0)
direction[direction==0]=-1
table(direction)
sym.p$lpv=sym.p$lpv*direction
write.csv(sym.p,file="symC_lpv.csv",row.names=F,quote=F)
save(sym.p,file="SymC_lpv.RData")

## save data for KEGG analysis 

# in the next line, choose contig2kegg.tab fom the transcriptome annotation files
keggTerms=read.table(file.choose(),sep="\t")
keggTerms=keggTerms[keggTerms$V1 %in% sym.p$gene,]
sym.p=sym.p[sym.p$gene %in% keggTerms$V1,]
keggTerms=keggTerms[match(sym.p$gene,keggTerms$V1),]
source("color4kegg.R")
kcols=color4kegg(sym.p$lpv)
keggColors=data.frame(cbind("term"=as.character(keggTerms$V2),"color"=kcols))
head(keggColors)

write.table(keggColors,quote=F,row.names=F,file="sym.kegg",sep="\t")
# go to http://www.genome.jp/kegg/tool/map_pathway2.html , upload mydata.kegg (check "use uncolored diagrams") 

###############################################################################
### KOG analysis ##### 

kog_genes=read.table("symC_iso2kogClass.tab", sep="\t")

library(KOGMWU)
data(gene2kog)

sym.lfc=final_results_sym[2]
sym.lfc=data.frame(sym.lfc)
genes=rownames(sym.lfc)
genes=data.frame(genes)
l2f=sym.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
write.csv(l2f, file="l2f.csv", row.names=F, quote=F)
sym.kog=kog.mwu(l2f, kog_genes , Alternative="t")
sym.kog 

