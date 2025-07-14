library(tidyverse)#for data wrangling
library(rhdf5)#for retreiving data from ARCHS4 database
library(edgeR)#for DGE analysis
library(gplots)#for heatmap.2
library(enrichplot)

archs4.human <- "human_gene_v2.latest.h5"
# use the h5 list (h5ls) function from the rhdf5 package to look at the contents of these databases
hdim <- h5ls(archs4.human)

# data for 67,186 HUMAN genes across 819,856 samples
all.samples.human <- h5read(archs4.human, name="meta/samples/geo_accession")
dim(all.samples.human)

#the sample identifiers are obtained from the GEO database for GSE 127831
mySamples <- c( "GSM3639530",	 #skin_HS1
                "GSM3639531",	#skin_HS2
                "GSM3639532",	#skin_HS3
                "GSM3639533",	#skin_HS4
                "GSM3639534",	#skin_HS5
                "GSM3639535",	#skin_HS6
                "GSM3639536",	#skin_HS7
                "GSM3639537",	#skin_CL1
                "GSM3639538",	#skin_CL2
                "GSM3639539",	#skin_CL3
                "GSM3639540",	#skin_CL4
                "GSM3639541",	#skin_CL5
                "GSM3639542",	#skin_CL6
                "GSM3639543",	#skin_CL7
                "GSM3639544",	#skin_CL8
                "GSM3639545",	#skin_CL9
                "GSM3639546",	#skin_CL10
                "GSM3639547",	#skin_CL11
                "GSM3639548",	#skin_CL12
                "GSM3639549",	#skin_CL13
                "GSM3639550",	#skin_CL14
                "GSM3639551",	#skin_CL15
                "GSM3639552",	#skin_CL16
                "GSM3639553",	#skin_CL17
                "GSM3639554",	#skin_CL18
                "GSM3639555",	#skin_CL19
                "GSM3639556",	#skin_CL20
                "GSM3639557")	#skin_CL21

# Identify columns to be extracted from ARCHS4 database
my.sample.locations <- which(all.samples.human %in% mySamples) # first time you've seen the %in% operator.

# extract gene symbols from the metadata
genes <- h5read(archs4.human, "meta/genes/symbol")

# Extract expression data from ARCHS4 ----
expression <- h5read(archs4.human, "data/expression",
                     index=list(my.sample.locations, NULL))
# transpose to get genes as rows and samples as columns
expression <- t(expression)

rownames(expression) <- genes
colnames(expression) <- all.samples.human[my.sample.locations]
colSums(expression) #this shows the sequencing depth for each of the samples you've extracted
archs4.dgelist <- DGEList(expression)
archs4.cpm <- cpm(archs4.dgelist)
colSums(archs4.cpm)
