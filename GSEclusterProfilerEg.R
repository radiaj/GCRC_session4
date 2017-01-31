#--------------------------------------------------------------------
# PART_2 R script: the basics of pathway analysis in R 
# By Radia Johnson, Ph.D.
# Date: Jan. 30th, 2017
#--------------------------------------------------------------------
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

#--------------------------------------------------------------------
# Pathway analysis with "ClusterProfiler"
#--------------------------------------------------------------------
# based on https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library("clusterProfiler")

# Load the fold changes from DESeq2 analysis and order in decreasing order
geneList = sort(exp.fc, decreasing = TRUE) # log FC is shown
head(geneList)

gene <- geneList[abs(geneList) >= 1] # Log 2 FC 
head(gene)

# GO over-representation test 
ego <- enrichGO(names(gene), 'org.Hs.eg.db', ont="BP", pvalueCutoff=1, qvalueCutoff=1)
head(ego) 
#head(summary(ego)) # older versions

# To see more results 
head(ego[, 1:5], 20) 
#head(summary(ego)[, 1:5], 20) # older versions

# Save a copy of the results
write.csv(ego, "hnrnp_EGO_results.csv") # Windows
#write.csv(summary(ego), "hnrnp_EGO_resultsUNIX.csv") # older versions

# Gene Set Enrichment Analysis of Gene Ontology 
egoGSEA <- gseGO(geneList, OrgDb='org.Hs.eg.db', ont="BP", pvalueCutoff=1, nPerm = 2) # Used few permutations to save time but I recommend running 100 or 1000 permutations (long computing time - recommend running of the cluster)
head(egoGSEA)
#head(summary(egoGSEA)) # older versions

#-------------------------------------------
# Visualizing clusterProfiler results
#-------------------------------------------
# For more information see http://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
barplot(ego, showCategory=10)
dotplot(ego)

# EnrichMap for ego results
enrichMap(ego, vertex.label.cex=0.7, layout=igraph::layout.kamada.kawai)

cnetplot(ego, categorySize="pvalue", foldChange=geneList)

# To visualize the cnetplot as an interactive plot
cnetplot(ego, categorySize="pvalue", foldChange=geneList, fixed=FALSE) # WORKS BETTER IN WINDOWS or the R console directly
# To convert ps to png https://cloudconvert.com/eps-to-png

# To plot a GSEAplot
gseaplot(egoGSEA, geneSetID = "GO:0001101")

#--------------------------------------------------------------------
# Reactome Pathway Analysis with the ReactomePA package
#--------------------------------------------------------------------
library("ReactomePA")

reacPA <- enrichPathway(gene=names(gene), organism = "human", pvalueCutoff = 0.05,
  pAdjustMethod = "BH", qvalueCutoff = 0.2, universe=names(geneList))
head(reacPA)

# Visualize the results
barplot(reacPA)

# Run GSEA analysis on Reactome database
y <- gsePathway(geneList, organism = "human", exponent = 1, nPerm = 1000,
  minGSSize = 10, maxGSSize = 500, pvalueCutoff =1,
  pAdjustMethod = "BH", verbose = TRUE, seed = FALSE, by = "fgsea")
res <- as.data.frame(y)
head(res)

# Visualize the results with enrichMap
enrichMap(y, layout=igraph::layout.kamada.kawai, vertex.label.cex = 0.7)

#-------------------------------------------
# Using the fgsea package
# For more information see http://www.bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
#-------------------------------------------
library("fgsea")
data(examplePathways)
data(exampleRanks)

head(examplePathways)
head(exampleRanks)

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=1000)
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]], exampleRanks) 

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)


