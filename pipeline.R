# Last Modified: 6 Jan 2022
##################################################################################################################################################################################################
# libraries
library(DiffBind)
library(tidyverse)
library(parallel)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.refGene)
library(EnsDb.Mmusculus.v79)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(diffloop)
##################################################################################################################################################################################################
# load metadata file
samples <- read.csv('metadata.csv')

# load samples
dbObj <- dba(sampleSheet=samples)
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)

# PCA analysis on the samples
dba.plotPCA(dbObj,  attributes=DBA_TREATMENT)
plot(dbObj)

# make a contrast matrix
dbObj <- dba.contrast(dbObj, categories=DBA_TREATMENT, minMembers = 3)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
dba.show(dbObj, bContrasts=T)

# choose package for differential analysis edgeR/DESeq2
dba.plotPCA(dbObj, contrast=1, method=DBA_EDGER, attributes=DBA_TREATMENT, label=DBA_ID)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)

# draw MA plot
dba.plotMA(dbObj, method=DBA_EDGER)
dba.plotMA(dbObj, bXY=TRUE)
dba.plotBox(dbObj)

# calculate p-valus
pvals <- res_edger <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)

# write files for downstream analysis
out <- as.data.frame(res_edger)
write.table(out, file="diffBind/KI_Ischemia_vs_WT_Ischemia_edger.txt", sep="\t", quote=F, row.names=F)

# Create bed files for each keeping only significant peaks (p < 0.05)
ki_ischemia_enrich <- out %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(ki_ischemia_enrich, file="diffBind/ki_ischemia_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)

wt_ischemia_enrich <- out %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

# Write to file
write.table(wt_ischemia_enrich, file="diffBind/wt_ischemia_enriched.bed", sep="\t", quote=F, row.names=F, col.names=F)
##################################################################################################################################################################################################
# load data
samplefiles <- list.files("chipSeeker", pattern= ".bed", full.names=T)
samplefiles <- as.list(samplefiles)
names(samplefiles) <- c("KI_Ischemia", "KI_Sham", "WT_Ischemia", "WT_Sham")
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene

peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")

ki_ischemia <- data.frame(peakAnnoList$KI_Ischemia)
gene_id_ki_ischemia <- ki_ischemia$geneId
annotations_edb_ki_ischemia <- AnnotationDbi::select(EnsDb.Mmusculus.v79,keys = gene_id_ki_ischemia,columns = c("GENENAME"),keytype = "ENTREZID")
annotations_edb_ki_ischemia$ENTREZID <- as.character(annotations_edb_ki_ischemia$ENTREZID)
ki_ischemia %>% 
  left_join(annotations_edb_ki_ischemia, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="chipSeeker/ki_ischemia_annotation.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis 
ego_ki_ischemia <- enrichGO(gene = gene_id_ki_ischemia, 
                            keyType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego_ki_ischemia)
write.csv(cluster_summary, "chipSeeker/clusterProfiler_ki_ischemia.csv")

# Dotplot visualization
dotplot(ego_ki_ischemia, showCategory=10)

ekegg_ki_ischemia <- enrichKEGG(gene = gene_id_ki_ischemia,
                    organism = 'mmu',
                    pvalueCutoff = 0.05)

dotplot(ekegg_ki_ischemia)


ki_sham <- data.frame(peakAnnoList$KI_Sham)
gene_id_ki_sham <- ki_sham$geneId
annotations_edb_ki_sham <- AnnotationDbi::select(EnsDb.Mmusculus.v79,keys = gene_id_ki_sham,columns = c("GENENAME"),keytype = "ENTREZID")
annotations_edb_ki_sham$ENTREZID <- as.character(annotations_edb_ki_sham$ENTREZID)
ki_sham %>% 
  left_join(annotations_edb_ki_sham, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="chipSeeker/ki_sham_annotation.txt", sep="\t", quote=F, row.names=F)


# Run GO enrichment analysis 
ego_ki_sham <- enrichGO(gene = gene_id_ki_sham, 
                            keyType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego_ki_sham)
write.csv(cluster_summary, "chipSeeker/clusterProfiler_ki_sham.csv")

# Dotplot visualization
dotplot(ego_ki_sham, showCategory=10)

ekegg_ki_sham <- enrichKEGG(gene = gene_id_ki_sham,
                                organism = 'mmu',
                                pvalueCutoff = 0.05)

dotplot(ekegg_ki_sham)

wt_ischemia <- data.frame(peakAnnoList$WT_Ischemia)
gene_id_wt_ischemia <- wt_ischemia$geneId
annotations_edb_wt_ischemia <- AnnotationDbi::select(EnsDb.Mmusculus.v79,keys = gene_id_wt_ischemia,columns = c("GENENAME"),keytype = "ENTREZID")
annotations_edb_wt_ischemia$ENTREZID <- as.character(annotations_edb_wt_ischemia$ENTREZID)
wt_ischemia %>% 
  left_join(annotations_edb_wt_ischemia, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="chipSeeker/wt_ischemia_annotation.txt", sep="\t", quote=F, row.names=F)


# Run GO enrichment analysis 
ego_wt_ischemia <- enrichGO(gene = gene_id_wt_ischemia, 
                            keyType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego_wt_ischemia)
write.csv(cluster_summary, "chipSeeker/clusterProfiler_wt_ischemia.csv")

# Dotplot visualization
dotplot(ego_wt_ischemia, showCategory=10)

ekegg_wt_ischemia <- enrichKEGG(gene = gene_id_wt_ischemia,
                                organism = 'mmu',
                                pvalueCutoff = 0.05)

dotplot(ekegg_wt_ischemia)

wt_sham <- data.frame(peakAnnoList$WT_Sham)
gene_id_wt_sham <- wt_sham$geneId
annotations_edb_wt_sham <- AnnotationDbi::select(EnsDb.Mmusculus.v79,keys = gene_id_wt_sham,columns = c("GENENAME"),keytype = "ENTREZID")
annotations_edb_wt_sham$ENTREZID <- as.character(annotations_edb_wt_sham$ENTREZID)
wt_sham %>% 
  left_join(annotations_edb_wt_sham, by=c("geneId"="ENTREZID")) %>% 
  write.table(file="chipSeeker/wt_sham_annotation.txt", sep="\t", quote=F, row.names=F)

# Run GO enrichment analysis 
ego_ki_ischemia <- enrichGO(gene = gene_id_ki_ischemia, 
                keyType = "ENTREZID", 
                OrgDb = org.Mm.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

# Run GO enrichment analysis 
ego_wt_sham <- enrichGO(gene = gene_id_wt_sham, 
                            keyType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

# Output results from GO analysis to a table
cluster_summary <- data.frame(ego_wt_sham)
write.csv(cluster_summary, "chipSeeker/clusterProfiler_wt_sham.csv")

# Dotplot visualization
dotplot(ego_wt_sham, showCategory=10)

ekegg_wt_sham <- enrichKEGG(gene = gene_id_wt_sham,
                                organism = 'mmu',
                                pvalueCutoff = 0.05)

dotplot(ekegg_wt_sham)


# Create a list with genes from each sample
genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)

# Run KEGG analysis
compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "mouse",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
