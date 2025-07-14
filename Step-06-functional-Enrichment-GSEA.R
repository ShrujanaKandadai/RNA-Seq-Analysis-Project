#functionalEnrichment part 3 (GSEA table)}

msigdbr_species()
hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
hs_gsea %>%
  dplyr::distinct(gs_collection, gs_subcollection) %>%
  dplyr::arrange(gs_collection, gs_subcollection)

# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in and pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
mydata.gsea <- mydata.gsea[!duplicated(names(mydata.gsea))]
# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
myGSEA.df <- as_tibble(myGSEA.res@result) #using @ symbol - a slot in an S4 class object

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)

#functionalEnrichment part 4 (enrich plot)}

# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = c(94,91), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[c(94,90)]
) #can also turn off this title

#functionalEnrichment part 4 (bubble plot)}

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
myGSEA.sig <- myGSEA.df %>% 
  filter(p.adjust < 0.01)

# Get top 10 pathways per phenotype
top_10_each <- myGSEA.sig %>%
  group_by(phenotype) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)
# Now: get unique pathways (to keep common y-axis)
top_IDs <- unique(top_10_each$ID)
# Subset full data for those IDs (to get both phenotypes for each)
top_df <- myGSEA.df %>%
  filter(ID %in% top_IDs)

# Optional: shorten pathway names for display
top_df$ID <- gsub("REACTOME_", "", top_df$ID)
top_df$ID <- gsub("_", " ", top_df$ID)

# Plot
ggplot(top_df, aes(x = phenotype, y = ID)) +
  geom_point(aes(size = setSize, color = NES, alpha = -log10(p.adjust))) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_alpha(range = c(0.4, 1)) +
  theme_bw() +
  labs(
    title = "Top GSEA Pathways by Phenotype",
    x = "Phenotype", y = "Pathway"
  ) +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )
