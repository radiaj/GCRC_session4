#--------------------------------------------------------------------
# PART_2 R script: the basics of pathway analysis in R 
# By Radia Johnson, Ph.D.
# Date: Jan. 26th, 2017
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Run GAGE for pathway analysis
#--------------------------------------------------------------------
# Running gage from DESEQ2 based on https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
library("gage")
data(kegg.gs)

library("gageData")
data(hnrnp.cnts)
cnts=hnrnp.cnts

head(cnts) # Entrez ids used for the gene ids

# Run the DESeq2 DEG analysis
library("DESeq2")
grp.idx <- rep(c("KO", "CTRL"), each=4)
coldat=DataFrame(grp=factor(grp.idx, levels=c("CTRL", "KO")))

dds <- DESeqDataSetFromMatrix(cnts, colData=coldat, design = ~ grp)
dds <- DESeq(dds)
deseq2.res <- results(dds)
deseq2.fc=deseq2.res$log2FoldChange
names(deseq2.fc)=rownames(deseq2.res)
exp.fc=deseq2.fc


# Run the GAGE pathway analysis
fc.kegg.p <- gage(exp.fc, gsets = kegg.gs, ref = NULL, samp = NULL)
sel <- fc.kegg.p$greater[, "q.val"] < 0.1 &
 !is.na(fc.kegg.p$greater[, "q.val"])
path.ids <- rownames(fc.kegg.p$greater)[sel]
sel.l <- fc.kegg.p$less[, "q.val"] < 0.1 &
 !is.na(fc.kegg.p$less[,"q.val"])
path.ids.l <- rownames(fc.kegg.p$less)[sel.l]

# Look at on overview of the GAGE pathway results
str(fc.kegg.p)
head(fc.kegg.p$greater)
head(fc.kegg.p$less)

# Look at the top UPREGULATED pathways 
head(path.ids)

# Look at the top DOWNREGULATED pathways 
head(path.ids.1)

#--------------------------------------------------------------------
# Use pathview to visualize the top 2 downregulated pathways
#--------------------------------------------------------------------
library("pathview")

# Get the KEGG id to search the database (8 first characters)
path.ids2 <- substr(c(path.ids, path.ids.l), 1, 8)
out.suffix="deseq2_hnrnp"

pv.out.list <- sapply(path.ids2[1:3], function(pid) pathview(
 gene.data = exp.fc, pathway.id = pid,
 species = "hsa", out.suffix=out.suffix))

