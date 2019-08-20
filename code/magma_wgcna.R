#!/usr/bin/env Rscript

# This script uses MAGMA to calculate gene-level significance values from GWAS summary statistics and then assesses module enrichment (code from [Gandal Science 2018](https://github.com/mgandal/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap/blob/master/code/04_networkEnrichment/4a_GWAS-Enrichment_MAGMA.R#L10) )

# This is a copy made 180407 of the original file /projects/tp/tmp-bmi-brain/src/magma_wgcna.R
# Changes: 
# 1) in cor.test, have set exact=F

# Missing
# 1) PT mapping
# 2) Save correlation values, pvalues and fdrs

# Usage

# time Rscript /projects/jonatan/tmp-bmi-brain/src/magma_wgcna.R maca "brain_pericyte" _kMEs.csv 180407
# time Rscript /projects/jonatan/tmp-bmi-brain/src/magma_wgcna.R maca "neuronal_stem_cell" _kMEs.csv 180407 

study_label = "campbell"
cell_label = "n12n13_7"
#cell_label = "brain_pericyte"
file_suffix = "_kMEs.csv"
#file_suffix = "_kMEs_PPI.csv"
output_label ="180409"

# Install some missing libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite("GO.db")
# biocLite("WGCNA")

# Command line arguments
rm(list=ls()); options(stringsAsFactors = F)
args = commandArgs(trailingOnly=TRUE)
if (length(args)<4) {
  stop("At least four arguments must be supplied (study label (e.g. fgf1_arc), cell label (e.g. astrocyte) file suffix (e.g. _hubgenes.csv) and output_label (e.g. 180315)). File sep is optional (default is ,)", call.=FALSE)
}  
study_label = args[1]
cell_label = args[2]
file_suffix = args[3]
output_label = args[4]
file_sep = ','

# Load libraries and set working directory

library(WGCNA); library(reshape); library(ggplot2); library(reshape2); library(plyr)
project_path = "/projects/jonatan/tmp-bmi-brain/"
magma_gwas_path = paste(project_path,"data/magma/",sep="/")
figs_path = paste(project_path,"figs/",sep="/")
mapping_hs_filepath = "/projects/tp/tmp-bmi-brain/data/mapping/gene_annotation_hsapiens.txt.gz"
mapping_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/Mus_musculus.GRCm38.90.gene_name_version2ensembl.txt.gz"
mapping_mm_synonyms_filepath = "/data/genetic-mapping/ncbi/Mus_musculus.gene_info_symbol2ensembl.gz"
mapping_hs_mm_filepath = "/projects/timshel/sc-genetics/sc-genetics/data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"

# Load WGCNA results
log_not_mapped_filepath = paste0(project_path,"log/",output_label,"_magma_wgcna_not_mapped_",study_label,"_",cell_label,".tab")
modulekME = read.csv(paste0(project_path,"data/wgcna/",study_label,"/",cell_label,file_suffix,sep=""),row.names=1, check.names = FALSE, sep = file_sep)

# Load MAGMA genes and remap to Ensembl gene IDs
d = dir(path=magma_gwas_path, pattern="[.]genes.out", recursive = T)
gwas = vector(mode="list")
for(i in 1:length(d)) {
  gwas[[i]] = read.table(paste(magma_gwas_path, d[[i]],sep=""),head=T, check.names = FALSE)
}
names(gwas) = gsub(".genes.out", "", d)

# Match and replace with ENSG 
#genes = union(gwas[[1]]$GENE, gwas[[2]]$GENE)
#for(i in 3:length(gwas)) genes = union(genes, gwas[[i]]$GENE)

# Remapping from human Entrez to human Ensembl gene IDs
mapping_hs_entrez2ensembl = read.csv(gzfile(mapping_hs_filepath),sep="\t",header=T)
for(i in 1:length(gwas)) {
	idx = match(gwas[[i]]$GENE, mapping_hs_entrez2ensembl$entrezgene)
	mapping = data.frame(entrez=gwas[[i]]$GENE, ensembl=mapping_hs_entrez2ensembl$ensembl_gene_id[idx])
	gwas[[i]]$gene_name = mapping$ensembl
}

# Remapping the WGCNA data to human ensembl IDs (using synonyms)
# Step 1: direct mapping
mapping_direct = read.table(gzfile(mapping_mm_filepath),sep="\t",header=T)
mapping = data.frame(symbol=rownames(modulekME), ensembl.mouse=mapping_direct$ensembl_gene_id[ match(rownames(modulekME), mapping_direct$gene_name_optimal) ])

# Step 2: map remaing using synonyms
mapping_synonyms = read.csv(gzfile(mapping_mm_synonyms_filepath),sep="\t",header=T)
mapping$ensembl.mouse[ which(is.na(mapping$ensembl.mouse)) ] = mapping_synonyms$ensembl[ match( mapping$symbol[which(is.na(mapping$ensembl.mouse)) ] ,mapping_synonyms$symbol) ]

# Step 3: orthology mapping
#mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
#mapping_mm_orthologs <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), mart=mart)
mapping_orthology = read.csv(gzfile(mapping_hs_mm_filepath),sep="\t",header=T)
mapping$ensembl.human = mapping_orthology$ensembl_gene_id[ match(mapping$ensembl.mouse,mapping_orthology$mmusculus_homolog_ensembl_gene) ]
#mapping$ensembl.human[mapping$ensembl.human == ""] = NA
df_not_mapped = mapping[is.na(mapping$ensembl.human),]
write.table(df_not_mapped,log_not_mapped_filepath,quote=F,sep="\t",row.names=F)
modulekME$symbol = mapping$symbol
modulekME$ensembl = mapping$ensembl.human
modulekME = na.omit(modulekME)
tmp = within(modulekME, rm("symbol","ensembl"))

# Average duplicated gene IDs
modulekME_ens <-aggregate(tmp, by=list(modulekME$ensembl),FUN=mean, na.rm=TRUE)
rownames(modulekME_ens) = modulekME_ens$Group.1
modulekME_ens = within(modulekME_ens, rm("Group.1"))

# Calculate spearman's correlation between gene module membership and GWAS gene significance
colors = colnames(modulekME_ens)
table.kme.cor.p = table.kme.cor.r<- matrix(NA,nrow=length(unique(colors)),ncol=length(gwas))
rownames(table.kme.cor.r) = rownames(table.kme.cor.p) = unique(colors)
colnames(table.kme.cor.r) = colnames(table.kme.cor.p) = names(gwas) 

for(m in unique(colors)) {
  for(i in 1:length(gwas)) {
    #col = paste("kME", m, sep="")
    col = m
    genes = intersect(rownames(modulekME_ens),gwas[[i]]$gene_name)
    x = -log10(gwas[[i]]$P[match(genes, gwas[[i]]$gene_name)])
    y = modulekME_ens[match(genes,rownames(modulekME_ens)), col]
    cor = cor.test(x,y,method="spearman", exact=F)
    table.kme.cor.r[m,i] = cor$estimate
    table.kme.cor.p[m,i] = cor$p.value
  }
}

table.kme.cor.p.fdr = p.adjust(table.kme.cor.p, method="fdr")
dim(table.kme.cor.p.fdr) = dim(table.kme.cor.p);  dimnames(table.kme.cor.p.fdr) = dimnames(table.kme.cor.p)

d = -log10(table.kme.cor.p.fdr) * sign(table.kme.cor.r) 
#pdf("SampleGraph.pdf",width=7,height=5)
#sizeGrWindow(9,7)
#par(mfrow = c(2,2))
#par(mar = c(4, 5, 4, 6));
labeledHeatmap(d,textMatrix = signif(table.kme.cor.r,1), xLabels = colnames(d), yLabels = rownames(d),invertColors = T, colors = blueWhiteRed(1000), main="GWAS - kME correlation", cex.text = 0.6)
#dev.off()

#dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))[c("tan","blue","yellow","purple","turquoise","green","greenyellow","salmon"),]
dat = as.data.frame(table.kme.cor.p.fdr*sign(table.kme.cor.r))
dat[dat<0]=1 #Only look for positive enrichment
dat = -log10(dat)
dat$module = gsub("kME","",rownames(dat))
#dat$module = gsub("ME","",rownames(dat))
dat2 = melt(dat)
dat2$variable=as.character(dat2$variable)

#p=ggplot(melt(dat),aes(x=variable,y=value,fill=colors)) + 
p=ggplot(melt(dat),aes(x=variable,y=value,fill="blue")) + 
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
  geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1))

p=ggplot(melt(dat),aes(x=variable,y=value,fill=module)) + 
  geom_bar(stat="identity",position=position_dodge(),color="black") +
  scale_fill_manual(values=sort(unique(dat2$module)))+ theme_classic() +
  geom_abline(intercept=-log10(0.05),slope=0,lty=2) + labs(x="",y="log10(P.fdr)") +
  theme(axis.text.x=element_text(angle=50, size=10, hjust=1))
p

ggsave(p, filename = paste0(figs_path,"/",output_label,"_wgcna_magma_",study_label,"_",cell_label,".pdf") ,width=45,height=12)
