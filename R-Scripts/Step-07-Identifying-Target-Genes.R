## Variability in gene expression

#identifying highly variable gene expression signatures in CL lesions Part 1 }

# Identification of highly variable gene expression signatures in CL lesions.
exprs_voom <- v.DEGList.filtered.norm$E
# 1. Get group-wise means
group_means <- data.frame(
  gene = rownames(exprs_voom),
  disease_mean = rowMeans(exprs_voom[, targets$disease == "cutaneous"]),
  healthy_mean = rowMeans(exprs_voom[, targets$disease == "control"])
)

# 2. Compute log2 fold change
group_means$logFC <- group_means$disease_mean - group_means$healthy_mean


# 3. Filter for genes:
# - high disease expression (e.g. > 5)
# - low healthy expression (e.g. < 1)
# - large FC (e.g. logFC > 4)
filtered_genes <- group_means %>%
  filter(disease_mean > 5, healthy_mean < 1, logFC > 4) %>%
  arrange(desc(logFC))

# View top candidates
head(filtered_genes, 10)

#removing the genes that start with "IG" to identify other differentially expressed genes that can be potential treatment targets 
# Remove genes that start with "IG"
filtered_genes_noIG <- filtered_genes %>% filter(!str_detect(gene,"^IG"))

#subset the expression matrix for genes of interest
# Keep only the first matching row for each gene
top5_genes <- head(filtered_genes_noIG, 5)
expr_subset <- exprs_voom[rownames(exprs_voom) %in% top5_genes$gene, ]

# Step 3: Convert to long format
df_long <- as.data.frame(expr_subset) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

# Step 4: Add sample metadata
df_long <- df_long %>%
  left_join(targets, by = "sample")

fc_table <- df_long %>%
  group_by(gene, disease) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = disease, values_from = mean_expr) %>%
  mutate(log2FC = cutaneous - control,
         FC = round(2^log2FC)) %>%
  arrange(desc(FC))  # Optional: sort

fc_labels <- paste0("FC ", fc_table$FC)

label_df <- df_long %>%
  filter(disease == "cutaneous") %>%
  group_by(gene) %>%
  summarise(y = max(expression) + 0.5) %>%
  mutate(label = fc_labels)

pd <- position_dodge(width = 0.75)
ggplot(df_long, aes(x = gene, y = expression, fill = disease)) +
  # Boxplot with consistent width and dodge
  geom_boxplot(position = pd, outlier.shape = NA, alpha = 0.8, width = 0.5, color = "black") +
  
  # Jitter points, aligned using jitterdodge
  geom_jitter(aes(color = disease),
              size = 1.8, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  
  # FC labels above the boxes
  geom_text(data = label_df, aes(x = gene, y = y, label = label),
            inherit.aes = FALSE,
            vjust = -0.5, size = 4.5, fontface = "bold") +
  
  # Color & fill mapping
  scale_fill_manual(values = c("control" = "gray30", "cutaneous" = "red")) +
  scale_color_manual(values = c("control" = "gray30", "cutaneous" = "red")) +
  
  # Labels and theme settings
  labs(y = "CPM (log2 scale)", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "italic", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  labs( title = "Top 5 most variable gene expression in CL vs. HS")


#identifying highly variable gene expression signatures in CL lesions Part 2 }

expr_matrix <- v.DEGList.filtered.norm$E
# Define list of known cytotoxicity genes
cytotoxic_genes <- c("GZMA", "GZMB", "PRF1", "NKG7", "GNLY", 
                     "CD8A", "CD8B", "FASLG", "IFNG", "TNF", 
                     "KLRK1", "CASP3", "BID", "IL2", "TRAIL","IL1B")

# Subset those present in your filtered gene list
cytotoxic_targets <- filtered_genes_noIG %>%
  filter(gene %in% cytotoxic_genes)

# View them sorted by logFC (descending)
cytotoxic_targets <- cytotoxic_targets %>%
  arrange(desc(logFC))
top5_cytotoxic_genes <- cytotoxic_targets %>%
  arrange(desc(logFC)) %>%
  slice(1:5) %>%
  pull(gene)

expr_subset_cyto <- expr_matrix[rownames(expr_matrix) %in% top5_cytotoxic_genes, ] 

df_long_cyto <- as.data.frame(expr_subset_cyto) %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression")

# Step 4: Add sample metadata
df_long_cyto <- df_long_cyto %>%
  left_join(targets, by = "sample")

# Calculate mean expression per gene per group
fc_table_cyto <- df_long_cyto %>%
  group_by(gene, disease) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = disease, values_from = mean_expr) %>%
  mutate(log2FC = cutaneous - control,
         FC = round(2^log2FC)) %>%
  arrange(desc(FC))  # Optional: sort

fc_labels_cyto <- paste0("FC ", fc_table_cyto$FC)
label_df_cyto <- df_long_cyto %>%
  filter(disease == "cutaneous") %>%
  group_by(gene) %>%
  summarise(y = max(expression) + 0.5) %>%
  mutate(label = fc_labels_cyto)

pd <- position_dodge(width = 0.75)
ggplot(df_long_cyto, aes(x = gene, y = expression, fill = disease)) +
  # Boxplot with consistent width and dodge
  geom_boxplot(position = pd, outlier.shape = NA, alpha = 0.8, width = 0.5, color = "black") +
  
  # Jitter points, aligned using jitterdodge
  geom_jitter(aes(color = disease),
              size = 1.8, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
  
  # FC labels above the boxes
  geom_text(data = label_df_cyto, aes(x = gene, y = y, label = label),
            inherit.aes = FALSE,
            vjust = -0.5, size = 4.5, fontface = "bold") +
  
  # Color & fill mapping
  scale_fill_manual(values = c("control" = "gray30", "cutaneous" = "red")) +
  scale_color_manual(values = c("control" = "gray30", "cutaneous" = "red")) +
  
  # Labels and theme settings
  labs(y = "CPM (log2 scale)", x = NULL) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray30") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.text.x = element_text(face = "italic", size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold")
  )+
  labs( title = "Top 5 most variable cytotoxic genes expressed in CL vs. HS")


#identifying highly variable gene expression signatures in CL lesions Part 2 }
# Identify columns to be extracted from ARCHS4 database
my.sample.locations1 <- which(all.samples.human %in% mySamples[8:28]) # first time you've seen the %in% operator.

# extract gene symbols from the metadata
genes <- h5read(archs4.human, "meta/genes/symbol")

# Extract expression data from ARCHS4 ----
expression1 <- h5read(archs4.human, "data/expression",
                      index=list(my.sample.locations1, NULL))
# transpose to get genes as rows and samples as columns
expression1 <- t(expression1)

rownames(expression1) <- genes
colnames(expression1) <- all.samples.human[my.sample.locations1]
colSums(expression1) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist_outcomes <- DGEList(expression1)
archs4.cpm_outcomes <- cpm(archs4.dgelist_outcomes)
colSums(archs4.cpm_outcomes)

# Filter and normalize the extracted data ----
table(rowSums(archs4.dgelist_outcomes$counts==0)==21)
keepers1 <- rowSums(archs4.cpm_outcomes>1)>=7
archs4.dgelist.filtered_outcomes <- archs4.dgelist_outcomes[keepers1,]
dim(archs4.dgelist.filtered_outcomes)
archs4.dgelist.filtered.norm_outcomes <- calcNormFactors(archs4.dgelist.filtered_outcomes, method = "TMM")

archs4.filtered.norm.log2.cpm_outcomes <- cpm(archs4.dgelist.filtered.norm_outcomes, log=TRUE)

# Extract sample metadata from ARCHS4 to create a study design file ----
# extract the sample source
sample_source_name <- h5read(archs4.human, "meta/samples/source_name_ch1")
# extract sample title
sample_title <- h5read(archs4.human, name="meta/samples/title")
# extract sample characteristics
sample_characteristics<- h5read(archs4.human, name="meta/samples/characteristics_ch1")

# let's try putting this all together in a study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      Sample_source = sample_source_name[my.sample.locations],
                      Sample_characteristics = sample_characteristics[my.sample.locations])

#based on what we extracted from ARCHS4 above, lets customize and clean-up this study design file
studyDesign <- tibble(Sample_title = sample_title[my.sample.locations],
                      genotype = rep(c("HS", "CL"), times= c(7,21)),
                      treatment = rep(c("NInf", "Inf"), times= c(7,21)))

#capture experimental variables as factors from this study design
genotype <- factor(studyDesign$genotype)
treatment <- factor(studyDesign$treatment)
sampleName <- studyDesign$Sample_title

targets.onlypatients <- targets[8:28,]

#failure vs cure as a factor 
outcome <- factor(targets.onlypatients$Treat..outcome)

#preparing a model matrix for CL vs HS 
design1 <- model.matrix(~0 + outcome)
colnames.design1 <- levels(outcome)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm_outcomes <- voom(archs4.dgelist.filtered.norm_outcomes, design1, plot = TRUE)
# fit a linear model to your data
fit1 <- lmFit(v.DEGList.filtered.norm_outcomes, design1)

# Contrast matrix ----
contrast.matrix.outcomes <- makeContrasts("failure vs. cure"= outcomefailure - outcomecure,
                                          levels=design1)

# extract the linear model fit -----
fits1 <- contrasts.fit(fit1, contrast.matrix.outcomes)
#get bayesian stats for your linear model fit
ebFit1 <- eBayes(fits1)
#write.fit(ebFit, file="lmfit_results.txt")

# TopTable to view DEGs -----
myTopHitsOutcomes <- topTable(ebFit1, adjust ="BH", coef=1, number=40000, sort.by="logFC")

# convert to a tibble
myTopHitsOutcomes.df <- myTopHitsOutcomes %>%
  as_tibble(rownames = "geneID")

# Add p-value significance categories
myTopHitsOutcomes_sig.df <- subset(myTopHitsOutcomes.df, P.Value < 0.05)
sigfailure <- myTopHitsOutcomes_sig.df$ID

myTopHitsOutcomes_sig.df$FC <- 2^(myTopHitsOutcomes_sig.df$logFC)
myTopHitsOutcomes_sig.df$neglog10P <- -log10(myTopHitsOutcomes_sig.df$P.Value)

cytotoxic_subset <- subset(myTopHitsOutcomes_sig.df, ID %in% cytotoxic_genes)

# View as table
# Display as interactive table
datatable(cytotoxic_subset,
          options = list(pageLength = 10, scrollX = TRUE),
          rownames = FALSE)

ggplot(myTopHitsOutcomes_sig.df, aes(x = FC, y = neglog10P)) +
  geom_point(color = "gray70") +
  
  geom_point(data = subset(myTopHitsOutcomes_sig.df, ID %in% cytotoxic_genes),
             aes(x = FC, y = neglog10P), color = "blue", size = 3) +
  
  geom_text(data = subset(myTopHitsOutcomes_sig.df, ID %in% cytotoxic_genes),
            aes(label = ID), color = "black", vjust = -1, size = 3.5) +
  
  xlim(0, 4) + ylim(0, 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  
  labs(
    x = "Fold Change (Failure vs. Cure)",
    y = "-log10(P.Value)",
    title = "Volcano Plot Highlighting Cytotoxic Genes"
  ) +
  theme_minimal(base_size = 14)

# Make sure this is the full DEG result (not just p < 0.05)
# myTopHitsOutcomes.df must contain columns: ID, FC, P.Value
myTopHitsOutcomes.df$neglog10P <- -log10(myTopHitsOutcomes.df$P.Value)
myTopHitsOutcomes.df$FC <- 2^(myTopHitsOutcomes.df$logFC)

library(ggplot2)

ggplot(myTopHitsOutcomes.df, aes(x = FC, y = neglog10P)) +
  # Plot all genes (gray)
  geom_point(color = "gray70") +
  
  # Highlight cytotoxic genes with p < 0.05 in blue
  geom_point(data = subset(myTopHitsOutcomes_sig.df, ID %in% cytotoxic_genes),
             aes(x = FC, y = neglog10P), color = "blue", size = 3) +
  
  # Label the cytotoxic genes
  geom_text(data = subset(myTopHitsOutcomes_sig.df, ID %in% cytotoxic_genes),
            aes(label = ID), color = "red", vjust = -1, size = 3.5) +
  
  # Axes limits and p-value threshold lines
  xlim(0, 4) + ylim(0, 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  
  labs(
    x = "Fold Change (Failure vs. Cure)",
    y = "-log10(P.Value)",
    title = "Volcano Plot Highlighting Cytotoxic Genes (P < 0.05)"
  ) +
  theme_minimal(base_size = 14)

