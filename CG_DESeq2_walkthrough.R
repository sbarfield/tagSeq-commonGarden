#DEseq2 for GBR common garden TagSeq (Acropora millepora) from Orpheus and Wilke Reefs

#installation of DESeq2
source("http://bioconductor.org/biocLite.R")
install.packages("BiocInstaller",repos="http://www.bioconductor.org/packages/3.1/bioc")
biocLite("DESeq2")
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
install.packages("vegan")

#read counts
counts=read.table("allcounts_new.txt", header=T, row.names=1)

#shorten names
names(counts)= sub(".sam.counts", "", names(counts))
ncol(counts)
nrow(counts)

#remove outliers- undersequenced samples
counts$o8a=NULL
counts$w22a=NULL
counts$w25c=NULL

names(counts)

# how many genes have mean count>2?
means=apply(counts,1,mean)
table(means>2)

# removing all genes with mean count less than 2 
counts=counts[means>2,]
nrow(counts)
counts=data.frame(counts)
#will be left with 19547 genes

#specify groups for analysis 
group=c(rep("Orph", 32), rep("WilkieC", 13), rep("WilkieD", 12))

#genotype
ind=c(rep("o0", 3), rep("o10",3), rep("o1",3), rep("o2",3), rep("o3",3), rep("o4",3), rep("o6",3), rep("o7",3), rep("o8",2), rep("o9",3), rep("oM1",3), rep("w22",2), rep("w23",3), rep("w24",3), rep("w25",2), rep("w26",3), rep("w27",3), rep("w28",3), rep("w30",3), rep("w31",3))

#ind.n
ind.n=c(rep(1,3), rep(2,3), rep(3,3), rep(4,3), rep(5,3), rep(6,3), rep(7,3), rep(8,2), rep(9,3), rep(10,3), rep(11,3), rep(1,2), rep(2,3), rep(3,3), rep(4,2), rep(5,3), rep(1,3), rep(2,3), rep(3,3), rep(4,3))

#replicate
replicate=c(rep(c(1,2,3), 8), 1,2, 1,2,3,1,2,3,1,2,1,2,3,1,2,3,1,2,1,2,3, rep(c(1,2,3), 4))

#location 
location=c(rep("Orpheus", 32), rep("Wilkie", 25))

#symbiont type 
symbiont=c(rep("cladeC", 45),  rep("cladeD", 12))

conditions=cbind(group, ind, ind.n, replicate, location, symbiont)
conditions=data.frame(conditions)

#model by individual then contrast by group
ddsi = DESeqDataSetFromMatrix(countData=counts, colData=conditions, design=~ind)
ddsi=DESeq(ddsi)
individual=results(ddsi)


wco=results(ddsi,contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieC"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="Orph"])),sep="")
),
listValues=c(0.5,-0.5)
)

wdo=results(ddsi,contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieD"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="Orph"])),sep="")
),
listValues=c(0.5,-0.5)
)

wcwd=results(ddsi,contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieC"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieD"])),sep="")
),
listValues=c(0.5,-0.5)
)

OvsW=results(ddsi, contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$location=="Orpheus"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$location=="Wilkie"])),sep="")
),
listValues=c(0.5,-0.5)
)

CvsD=results(ddsi, contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$symbiont=="cladeC"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$symbiont=="cladeD"])),sep="")
),
listValues=c(0.5,-0.5)
)



# plot PCA with this data
vsd=varianceStabilizingTransformation(ddsi)

quartz()

plotPCA = function (x, intgroup = "group", ntop = 500) 
{
  rv = rowVars(assay(x))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select, ]))
  fac = factor(apply(as.data.frame(colData(x)[, intgroup, drop = FALSE]), 
                     1, paste, collapse = " : "))
  if (nlevels(fac) > 4) 
    #colours = brewer.pal(nlevels(fac), "Paired")
    colours=c("cornflowerblue", "azure4", "chocolate2", "firebrick", "darkolivegreen2", "deepskyblue2", "dimgrey", "goldenrod", "darkorchid1", 
              "gray1", "mediumpurple1", "magenta", "mistyrose", "midnightblue", "orange", "seagreen1", "royalblue", "yellow", "tan3", "turquoise2")
  else colours = c("red", "dodgerblue", "green")
  xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
         pch = 2, cex = 1, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
                                                                                     text = list(levels(fac)), rep = FALSE)))
}

plotPCA(vsd, intgroup="group", ntop=500)


#variance stabilized gene counts 
vsd=getVarianceStabilizedData(ddsi)


#log2 gene counts (rlog function)

#log2 transformed data
rlog_group=rlog(ddsi)
rlog_group=assay(rlog_group)
save(rlog_group, file="rlog_group.RData")


#simulated datasets 
sim=simulateCounts(ddsi)

#run DESeq on simulated data set
simNew=DESeq(sim)

#get simulated results

sim_wco=results(simNew, contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieC"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="Orph"])),sep="")
),
listValues=c(0.5,-0.5)
)

sim_wdo=results(simNew, contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieD"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="Orph"])),sep="")
),
listValues=c(0.5,-0.5)
)

sim_wcwd=results(simNew, contrast=list(
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieC"])),sep=""),
  paste("ind",as.character(unique(conditions$ind[conditions$group=="WilkieD"])),sep="")
),
listValues=c(0.5,-0.5)
)

save(conditions, vsd, wco, wdo, wcwd, file="DESeq2_results.RData")
save(conditions, sim_wco, sim_wdo, sim_wcwd, simInd, file="DESeq2_sim_results.RData")
library(DESeq2)
library(empiricalFDR.DESeq2)
library(pheatmap)

#Wilkie C vs Orph
table(wco$padj<0.1) # filtering based DEGs
table(is.na(wco$padj)) # this is how many needed to be tossed
fdrt=fdrTable(wco$pvalue,sim_wco$pvalue)
par(mfrow=c(1,2))
fdrBiCurve(fdrt,main="WCO")
quartz()
efdrS=empiricalFDR(fdrt,plot=T,main="WCO")
mtext(paste("10% FDR at p =",signif(efdrS,2)),cex=0.8)
table(wco$pval<efdrS) # empirical FDR based DEGs
wco$efdr=(wco$pval<efdrS)
wco$efdr[is.na(wco$efdr)]=FALSE
#10% FDR is at p = 0.021

#Wilkie D vs Orph
table(wdo$padj<0.1) # filtering based DEGs
table(is.na(wdo$padj)) # this is how many needed to be tossed
fdrt=fdrTable(wdo$pvalue,sim_wdo$pvalue)
par(mfrow=c(1,2))
fdrBiCurve(fdrt,main="WDO")
quartz()
efdrS=empiricalFDR(fdrt,plot=T,main="WDO")
mtext(paste("10% FDR at p =",signif(efdrS,2)),cex=0.8)
table(wdo$pval<efdrS) # empirical FDR based DEGs
wdo$efdr=(wdo$pval<efdrS)
wdo$efdr[is.na(wdo$efdr)]=FALSE
#10% FDR is at p = 0.021

#Wilkie C vs Wilkie D
table(wcwd$padj<0.1) # filtering based DEGs
table(is.na(wcwd$padj)) # this is how many needed to be tossed
fdrt=fdrTable(wcwd$pvalue,sim_wcwd$pvalue)
par(mfrow=c(1,2))
fdrBiCurve(fdrt,main="WCWD")
quartz()
efdrS=empiricalFDR(fdrt,plot=T,main="WCWD")
mtext(paste("10% FDR at p =",signif(efdrS,2)),cex=0.8)
table(wcwd$pval<efdrS) # empirical FDR based DEGs
wcwd$efdr=(wcwd$pval<efdrS)
wcwd$efdr[is.na(wcwd$efdr)]=FALSE
#10% FDR is at p = 0.022

#Individual
table(individual$padj<0.1) # filtering based DEGs
table(is.na(individual$padj)) # this is how many needed to be tossed
fdrt=fdrTable(individual$pvalue,simInd$pvalue)
par(mfrow=c(1,2))
fdrBiCurve(fdrt,main="Ind")
quartz()
efdrS=empiricalFDR(fdrt,plot=T,main="Ind")
mtext(paste("10% FDR at p =",signif(efdrS,2)),cex=0.8)
table(individual$pval<efdrS) # empirical FDR based DEGs
individual$efdr=(individual$pval<efdrS)
individual$efdr[is.na(individual$efdr)]=FALSE
#10% FDR is at p = 0.022


#plot sample heat maps 
quartz()
degs_efdr_WCO=which(wco$efdr)
length(degs_efdr_WCO) 
pheatmap(cor(vsd_group[degs_efdr_WCO,]),border_color=NA) 

degs_efdr_WDO=which(wdo$efdr)
length(degs_efdr_WDO) 

degs_efdr_WCWD=which(wcwd$efdr)
length(degs_efdr_WCWD) 

pheatmap(cor(vsd), corder_color=NA)

### principal coordinates analysis 
library(ape)
ad.pcoa=pcoa(dist(t(vsd),method="manhattan")/1000)
scores=ad.pcoa$vectors

#code by group 
quartz()
plot(scores[,1], scores[,2],pch=as.numeric(as.factor(conditions$ind)),col=as.numeric(as.factor(conditions$group)))
ordispider(scores, conditions$ind)
ordiellipse(scores, conditions$group, label=T)

#analysis of variance 
adonis(t(vsd)~group,data=conditions,method="manhattan")

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#Df  SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
#group      2  154498154 77249077  2.2287 0.07625  0.001 ***
#  Residuals 54 1871730104 34661669         0.92375           
#Total     56 2026228258                  1.00000           
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#sorting and saving DEGS that pass 10% FDR cutoff 

WCO_fdr=subset(wco, pvalue < 0.021)
WCO_FDR=WCO_fdr[order(WCO_fdr$pvalue), , drop=FALSE]
WCO_FDR_sorted=WCO_FDR[2]
degs=row.names(WCO_FDR_sorted)
degs_WCO=row.names(WCO_FDR_sorted)
sorted=vsd[degs,]
write.table(degs, file="degs_WCO_NEW.txt")
write.table(WCO_FDR, file="wco_degs.txt", sep="\t")

WDO_fdr=subset(wdo, pvalue < 0.021)
WDO_FDR=WDO_fdr[order(WDO_fdr$pvalue), , drop=FALSE]
WDO_FDR_sorted=WDO_FDR[2]
degs=row.names(WDO_FDR_sorted)
degs_WDO=row.names(WDO_FDR_sorted)
sorted=vsd[degs,]
write.table(degs, file="degs_WDO_NEW.txt")
write.table(WDO_FDR, file="wdo_degs.txt", sep="\t")

WCWD_fdr=subset(wcwd, pvalue < 0.022)
WCWD_FDR=WCWD_fdr[order(WCWD_fdr$pvalue), , drop=FALSE]
WCWD_FDR_sorted=WCWD_FDR[2]
degs=row.names(WCWD_FDR_sorted)
degs_WCWD=row.names(WCWD_FDR_sorted)
sorted=vsd[degs,]
write.table(degs, file="degs_WCWD_NEW.txt")
write.table(WCWD_FDR, file="wcwd_degs.txt", sep="\t")


#make a heat map with top DEGs (with attached gene names)

gg=read.table("amil_iso2gene.tab",sep="\t")
library(pheatmap)
library(gplots)
library(RColorBrewer)

top=50

select=order(rowMeans(counts(cds)), decreasing=TRUE)[1:100]
hmcol=colorRampPalette(c("red", "white", "blue"))
heatmap.2(exprs(sorted), col=hmcol, trace="none")

my_palette <- colorRampPalette(c("gold", "black", "turquoise"))(n = 100)

heatmap(sorted, col = my_palette, Colv=NA)

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


#constructing a venn diagram to compare various models 
#selecting 10% efdr DEGS 

DEGS_wco=row.names(wco) [which(wco$pvalue<0.021 & !is.na(wco$pvalue))]
DEGS_wdo=row.names(wdo) [which(wdo$pvalue<0.021 & !is.na(wdo$pvalue))]
DEGS_wcwd=row.names(wcwd) [which(wcwd$pvalue<0.022 & !is.na(wcwd$pvalue))]

candidates=list("Wilkie C vs Orpheus"=DEGS_wco, "Wilkie D vs Orpheus"=DEGS_wdo, "Wilkie C vs D"=DEGS_wcwd)
quartz()
venn(candidates)

# pretty venn diagram
prettyvenn=venn.diagram(
  x = candidates,
  filename=NULL,
  col = "transparent",
  fill = c("coral", "cyan3", "lightgreen"),
  alpha = 0.5,
  label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
  cex = 2.0,
  fontfamily = "sans",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue", "darkgreen"),
  cat.cex = 1.0,
  cat.fontfamily = "sans",
  cat.dist = c(0.06, 0.06, 0.03),
  cat.pos = 1
);
quartz()
grid.draw(prettyvenn)

### proportional venn diagram ## 
require(venneuler)
v <- venneuler(c(A=1197, B=945, "A&B"=1197, C=1372, "A&C"=867, "B&C"=727, "A&B&C"=251))
quartz()
plot(v)


#KOG-MWU test for functional enrichment 
kog_genes=read.table("amil_defog_iso2kogClass.tab", sep="\t")

library(KOGMWU)
data(gene2kog)

#Wilkie C vs Orpheus 
wco.lfc=wco[2]
wco.lfc=data.frame(wco.lfc)
genes=rownames(wco.lfc)
genes=data.frame(genes)
l2f=wco.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
wco.kog=kog.mwu(l2f, gene2kog, Alternative="t")
wco.kog 


#Wilkie D vs Orpheus 
wdo.lfc=wdo[2]
wdo.lfc=data.frame(wdo.lfc)
genes=rownames(wdo.lfc)
genes=data.frame(genes)
l2f=wdo.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
wdo.kog=kog.mwu(l2f, gene2kog, Alternative="t")
wdo.kog


#Wilkie C vs Wilkie D 
wcwd.lfc=wcwd[2]
wcwd.lfc=data.frame(wcwd.lfc)
genes=rownames(wcwd.lfc)
genes=data.frame(genes)
l2f=wcwd.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
wcwd.kog=kog.mwu(l2f, gene2kog, Alternative="t")
wcwd.kog


ow.lfc=OvsW[2]
ow.lfc=data.frame(ow.lfc)
genes=rownames(ow.lfc)
genes=data.frame(genes)
l2f=ow.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
ow.kog=kog.mwu(l2f, gene2kog, Alternative="t")
ow.kog


cd.lfc=CvsD[2]
cd.lfc=data.frame(cd.lfc)
genes=rownames(cd.lfc)
genes=data.frame(genes)
l2f=cd.lfc$log2FoldChange
l2f=data.frame(l2f)
l2f=cbind(genes, l2f)
colnames(l2f)=c("gene", "lfc")
cd.kog=kog.mwu(l2f, gene2kog, Alternative="t")
cd.kog


#### saving data for GO analysis 

#signed (-log(pvalue)) - wilkie c vs orpheus 
wco.fc=data.frame("gene"=row.names(wco[!is.na(wco$pvalue),]))
wco.fc$lfc=wco[!is.na(wco$pvalue),"log2FoldChange"]
wco.p=data.frame("gene"=row.names(wco[!is.na(wco$pvalue),]))
wco.p$lpv=-log(wco[!is.na(wco$pvalue),"pvalue"]+1e-15,10)
direction=as.numeric(wco.fc[,"lfc"]>0)
direction[direction==0]=-1
table(direction)
wco.p$lpv=wco.p$lpv*direction
write.csv(wco.p,file="wco_lpv.csv",row.names=F,quote=F)
save(wco.p,file="wco_lpv.RData")

#log fold change 
wco.fc=data.frame("gene"=row.names(wco[!is.na(wco$log2FoldChange),]))
wco.fc$lfc=wco[!is.na(wco$log2FoldChange),"log2FoldChange"]
wco.p=data.frame("gene"=row.names(wco[!is.na(wco$log2FoldChange),]))
wco.p$lfc=wco$log2FoldChange
write.csv(wco.p,file="wco_l2f.csv",row.names=F,quote=F)
save(wco.p,file="wco_l2f.RData")

wco.fc=data.frame("gene"=row.names(final_results_WCO[!is.na(final_results_WCO$log2FoldChange),]))
wco.fc$lfc=final_results_WCO[!is.na(final_results_WCO$log2FoldChange),"log2FoldChange"]
wco.p=data.frame("gene"=row.names(final_results_WCO[!is.na(final_results_WCO$log2FoldChange),]))
wco.p$lfc=final_results_WCO$log2FoldChange
write.csv(wco.p,file="wco_l2f_old.csv",row.names=F,quote=F)
save(wco.p,file="wco_l2f_old.RData")

#signed (-log(pvalue)) - wilkie C vs. wilkie D
wcwd.fc=data.frame("gene"=row.names(final_results_WCWD[!is.na(final_results_WCWD$pvalue),]))
wcwd.fc$lfc=final_results_WCWD[!is.na(final_results_WCWD$pvalue),"log2FoldChange"]
wcwd.p=data.frame("gene"=row.names(final_results_WCWD[!is.na(final_results_WCWD$pvalue),]))
wcwd.p$lpv=-log(final_results_WCWD[!is.na(final_results_WCWD$pvalue),"pvalue"]+1e-15,10)
direction=as.numeric(wcwd.fc[,"lfc"]>0)
direction[direction==0]=-1
table(direction)
wcwd.p$lpv=wcwd.p$lpv*direction
write.csv(wcwd.p,file="wcwd_lpv.csv",row.names=F,quote=F)
save(wcwd.p,file="wcwd_lpv.RData")

#log fold change 
wcwd.fc=data.frame("gene"=row.names(wcwd[!is.na(wcwd$log2FoldChange),]))
wcwd.fc$lfc=wcwd[!is.na(wcwd$log2FoldChange),"log2FoldChange"]
wcwd.p=data.frame("gene"=row.names(wcwd[!is.na(wcwd$log2FoldChange),]))
wcwd.p$lfc=wcwd$log2FoldChange
write.csv(wcwd.p,file="wcwd_l2f.csv",row.names=F,quote=F)
save(wcwd.p,file="wcwd_l2f.RData")

#signed (-log(pvalue)) - wilkie D vs orpheus
wdo.fc=data.frame("gene"=row.names(final_results_WDO[!is.na(final_results_WDO$pvalue),]))
wdo.fc$lfc=final_results_WDO[!is.na(final_results_WDO$pvalue),"log2FoldChange"]
wdo.p=data.frame("gene"=row.names(final_results_WDO[!is.na(final_results_WDO$pvalue),]))
wdo.p$lpv=-log(final_results_WDO[!is.na(final_results_WDO$pvalue),"pvalue"]+1e-15,10)
direction=as.numeric(wdo.fc[,"lfc"]>0)
direction[direction==0]=-1
table(direction)
wdo.p$lpv=wdo.p$lpv*direction
write.csv(wdo.p,file="wdo_lpv.csv",row.names=F,quote=F)
save(wdo.p,file="wdo_lpv.RData")

#log fold change 
wdo.fc=data.frame("gene"=row.names(wdo[!is.na(wdo$log2FoldChange),]))
wdo.fc$lfc=wdo[!is.na(wdo$log2FoldChange),"log2FoldChange"]
wdo.p=data.frame("gene"=row.names(wdo[!is.na(wdo$log2FoldChange),]))
wdo.p$lfc=wdo$log2FoldChange
write.csv(wdo.p,file="wdo_l2f.csv",row.names=F,quote=F)
save(wdo.p,file="wdo_l2f.RData")

#OvsW log fold change for go
ow.fc=data.frame("gene"=row.names(OvsW[!is.na(OvsW$log2FoldChange),]))
ow.fc$lfc=OvsW[!is.na(OvsW$log2FoldChange),"log2FoldChange"]
ow.p=data.frame("gene"=row.names(OvsW[!is.na(OvsW$log2FoldChange),]))
ow.p$lfc=OvsW$log2FoldChange
write.csv(ow.p,file="ow_l2f.csv",row.names=F,quote=F)
save(ow.p,file="ow_l2f.RData")

#CvsD log fold change 
cd.fc=data.frame("gene"=row.names(CvsD[!is.na(CvsD$log2FoldChange),]))
cd.fc$lfc=CvsD[!is.na(CvsD$log2FoldChange),"log2FoldChange"]
cd.p=data.frame("gene"=row.names(CvsD[!is.na(CvsD$log2FoldChange),]))
cd.p$lfc=CvsD$log2FoldChange
write.csv(cd.p,file="cd_l2f.csv",row.names=F,quote=F)
save(cd.p,file="cd_l2f.RData")


