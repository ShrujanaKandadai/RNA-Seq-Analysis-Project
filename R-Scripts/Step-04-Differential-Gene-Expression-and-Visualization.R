#reading in the study design file 
design <- read.csv(gzfile("GSE127831_studydesign.csv.gz"))
write.table(design, "GSE127831_studydesign.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#HS vs. CL lesions as a factor 
disease <- factor(targets$disease)

#preparing a model matrix for CL vs HS 
design1 <- model.matrix(~0 + disease)
colnames.design1 <- levels(disease)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(archs4.dgelist.filtered.norm, design1, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design1)

# Contrast matrix ----
contrast.matrix.disease <- makeContrasts("infection vs. control"= diseasecutaneous - diseasecontrol,
                                         levels=design1)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix.disease)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

# now plot
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", ID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot)

#extract the DEGs
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=1)
sampleLabels <- targets$sample
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
datatable(diffGenes.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: DEGs in cutaneous leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#clustering - using unsupervised method: correlation
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

#we'll look at these clusters in more detail later
module.assign <- cutree(clustRows, k=2)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
myheatcolors2 <- colorRampPalette(c("purple", "white", "darkgreen"))(50)
myheatcolors1 <- bluered(75)
# Produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(myheatcolors1), scale='row', labRow=rownames(diffGenes),
          density.info="none", trace="none",  
          cexRow = 0.8,                         # Slightly smaller gene labels
          cexCol = 0.7,                         # Smaller column labels
          srtCol = 45,                          # Rotate column labels
          adjCol = c(1, 1),                     # Right-justify rotated labels
          offsetCol = 0.5,                      # Better spacing
          margins = c(10, 20))                  # Increase bottom margin) 

# View modules of co-regulated genes ----
# view your color assignments for the different clusters
names(module.color) <- names(module.assign) 

module.assign.df <- tibble(
  geneID = names(module.assign),
  module = as.vector(module.assign)
)

module.assign.df <- module.assign.df %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#FF9900",
    module == 2 ~ "#FF0099"))


ggplot(module.assign.df) +
  aes(module) +
  geom_bar(aes(fill=moduleColor)) +
  theme_bw()

#using modules to create a heatmap of genes that are upregulated in disease 
#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
modulePick <- 2 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)
#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub <- hclust(as.dist(1-cor(t(myModule), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule, 
          Rowv=as.dendrogram(hrsub), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors1), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], 
          cexRow = 0.8,                         # Slightly smaller gene labels
          cexCol = 0.7,                         # Smaller column labels
          srtCol = 45,                          # Rotate column labels
          adjCol = c(1, 1),                     # Right-justify rotated labels
          offsetCol = 0.5,                      # Better spacing
          margins = c(10, 20))                  # Increase bottom margin) 

#using modules to create a heatmap of genes that are downregualted in disease 
modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=rev(myheatcolors1), scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick],
          cexRow = 0.8,                         # Slightly smaller gene labels
          cexCol = 0.7,                         # Smaller column labels
          srtCol = 45,                          # Rotate column labels
          adjCol = c(1, 1),                     # Right-justify rotated labels
          offsetCol = 0.5,                      # Better spacing
          margins = c(10, 20))                  # Increase bottom margin


#Heatmap of Top 100 genes upregulated in CL vs. HS)}

sortedupFC <- myTopHits[order(-myTopHits$logFC),]
sortedupFC <- sortedupFC[!duplicated(sortedupFC$ID) & !is.na(sortedupFC$ID), ]  # avoid duplicates
rownames(sortedupFC) <- sortedupFC$ID
sorteduptop <- sortedupFC[-grep("IG", sortedupFC$ID),]
#Top 100 genes upregulated:
top100up <- sorteduptop[1:100,]$ID
TopUPtable100 <- archs4.filtered.norm.log2.cpm[c(top100up),]
TopUPmatrixcoding100 <- as.matrix(TopUPtable100)
colormapX <- colorRampPalette(colors=c("darkgreen","white","purple"))(50)
HeatmapUP100 <- heatmap.2(TopUPmatrixcoding100,
                          scale = "row", key=TRUE,
                          keysize = 1, key.title = NA,
                          col=colormapX, dendrogram = "none", Rowv = F,
                          margins=c(5,25),
                          labCol = NA, labRow = rownames(TopUPtable100),
                          main = "",
                          density.info="none", trace="none",
                          cexRow=0.8, cexCol=1)
# Save heatmap of top 100 upregulated coding genes to PNG file
png("Top100_UP_heatmap.png", width=1200, height=1800) #width/height in pixels
heatmap.2(TopUPmatrixcoding100,
          scale = "row", key=TRUE,
          keysize = 1, key.title = NA,
          col=colormapX, dendrogram = "none", Rowv = F,
          margins=c(5,25),
          labCol = NA, labRow = rownames(TopUPmatrixcoding100),
          main = "",
          density.info="none", trace="none",
          cexRow=0.7, cexCol=1)
dev.off()


