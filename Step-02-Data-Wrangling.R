table(rowSums(archs4.dgelist$counts==0)==28)
keepers <- rowSums(archs4.cpm>1)>=7
archs4.dgelist.filtered <- archs4.dgelist[keepers,]
dim(archs4.dgelist.filtered)
archs4.dgelist.filtered.norm <- calcNormFactors(archs4.dgelist.filtered, method = "TMM")

archs4.filtered.norm.log2.cpm <- cpm(archs4.dgelist.filtered.norm, log=TRUE)

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

#Filtering was carried out to remove lowly expressed genes. 
#Genes with less than 1 count per million (CPM) in at least 7 or more samples filtered out.  
#This reduced the number of genes from `r nrow(archs4.dgelist)` to `r nrow(archs4.dgelist.filtered)`. 

library(tidyverse)
library(DT)

# use dplyr 'mutate' function to add new columns based on existing data
targets <- read_tsv("GSE127831_studydesign.txt")
sampleLabels <- targets$sample
colnames(archs4.filtered.norm.log2.cpm) <- c(sampleLabels)
archs4.filtered.norm.log2.cpm.df <- as_tibble(archs4.filtered.norm.log2.cpm, rownames = "geneID")
archs4.filtered.norm.log2.cpm.df
mydata.df <- archs4.filtered.norm.log2.cpm.df %>% 
  mutate(healthy.AVG = (HS1 + HS2 + HS3 + HS4 + HS5 + HS6 +HS7)/7,
         disease.AVG = (CL1 + CL2 + CL3 + CL4 + CL5 + CL6 + CL7 + CL8 + CL9 + CL10 + CL11 + CL12 + CL13 + CL14 + CL15 + CL16 + CL17 + CL18 + CL19 + CL20 + CL21)/21,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (disease.AVG - healthy.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

#now look at this modified data table
mydata.df

datatable(mydata.df[,c(1,30:32)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

