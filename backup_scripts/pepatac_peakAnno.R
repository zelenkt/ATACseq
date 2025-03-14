#!/usr/bin/Rscript

##### LAST UPDATED 02/20/2023 Tomas Zelenka #####

##### LOAD LIBRARIES #####
suppressMessages(
if  (!require("DESeq2")) {
    install.packages("DESeq2", dependencies = TRUE)
    library(DESeq2)
    }
)

suppressMessages(
if  (!require("Rsubread")) {
    install.packages("Rsubread", dependencies = TRUE)
    library(Rsubread)
    }
)

suppressMessages(
if  (!require("getopt")) {
    install.packages("getopt", dependencies = TRUE)
    library(getopt)
    }
)
library(ChIPseeker)
library(clusterProfiler)
library(ggupset)
library(ggimage)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DOSE)
library(ReactomePA)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
files <- list("C4AP" = "/Users/4476125/Documents/datasets/06_ATAC_analysis_ChIPseeker/input_DAGs/C4AP/all_differential_peaks_noFilt_for_annotation.bed", "C8AP" = "/Users/4476125/Documents/datasets/06_ATAC_analysis_ChIPseeker/input_DAGs/C8AP/all_differential_peaks_noFilt_for_annotation.bed", "C8AT" = "/Users/4476125/Documents/datasets/06_ATAC_analysis_ChIPseeker/input_DAGs/C8AT/all_differential_peaks_noFilt_for_annotation.bed", "merged_allCD48" = "/Users/4476125/Documents/datasets/06_ATAC_analysis_ChIPseeker/input_DAGs/merged_allCD48/all_differential_peaks_noFilt_for_annotation.bed", "merged_onlyCD8" = "/Users/4476125/Documents/datasets/06_ATAC_analysis_ChIPseeker/input_DAGs/merged_onlyCD8/all_differential_peaks_noFilt_for_annotation.bed")

peakAnno <- annotatePeak(files[[5]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

# from the list of peakAnno I can also extract many aspects - such as geneId that is used for GO enrichment above...also ensembl, gene description and gene names, genomic features, respectively 

# this creates a dataframe of all useful peak annotation data
C4AP_annotated_all_noFilt <- data.frame(peakAnno@anno@elementMetadata@listData)

write.table(C4AP_annotated_all_noFilt, 
           file = 'C4AP_annotated_all_noFilt.txt', 
           quote = F, 
           sep = '\t')



##### GET ARGUMENTS AND SET VARIABLES #####

spec = matrix(c(
    'condition' , 'c', 1, "character", "Test condition (default: 'genotype')",
    'baseline'  , 'b', 1, "character", "Control group identifier (default: 'WT')",
    'samples'       , 's', 1, "character", "ATAC-seq / ChIP-seq sample list file (default: 'test-samples.csv')",
    'readcounts'       , 'r', 1, "character", "Raw read counts file (default: 'raw_read_counts_pepatac_consensus_peaks.txt')"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# OPTIONAL ARGS
if ( is.null(opt$condition ) ) { opt$condition = "genotype"                                             }
if ( is.null(opt$baseline  ) ) { opt$baseline  = "WT"                                                   }
if ( is.null(opt$samples       ) ) { opt$samples       = paste0(path.expand("/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/samples_C4P.csv"))        }
if ( is.null(opt$readcounts       ) ) { opt$readcounts       = paste0(path.expand("raw_read_counts_pepatac_consensus_peaks.txt"))        }


##### BEGIN SCRIPT #####
    myCounts <- as.matrix(read.table(opt$readcounts, sep = "\t", row.names=1, header=1))
    #myCounts <- read.table(opt$readcounts, sep = "\t", row.names=1, header=1)
    #myCounts <- as.matrix(read_tsv(opt$readcounts, sep = "\t", row.names=1, header=1))
#read_tsv()
    #colnames(myCounts) <- colnames(myCounts)[2:ncol(myCounts)]

## IMPORT SAMPLE INFO TABLE
    sampleInfo <- read.table(opt$samples, sep = ",", row.names=1, header=1)

## REMOVE UNNECCESARY SAMPLEINFO SAMPLE DATA AND REORDER TO MATCH "myCounts"
    sortedIndices <- match(colnames(myCounts), rownames(sampleInfo))
    sampleInfo <- sampleInfo[sortedIndices,,drop=FALSE]

## INPUT DATA INTO DESeq2
    dds <- DESeqDataSetFromMatrix(countData = myCounts,
                                  colData = sampleInfo,
                                  design = as.formula(sprintf("~ %s",opt$condition)))

## SPECIFY REFERENCE LEVEL(GENOTYPE)
    dds[[opt$condition]] <- relevel(dds[[opt$condition]], ref=opt$baseline)

## MINIMAL DATA FILTERING
    dds <- dds[ rowSums(counts(dds)) > 1, ]

## PERFORM DIFFERENTIAL ANALYSIS
    dds <- DESeq(dds)
    res <- results(dds)
    resOrdered <- res[order(res$padj),]
    normalizedCounts <- counts(dds, normalized=TRUE)

## EXPORT FILES
    write.table(resOrdered, file = "DESeq2_analysis.txt", sep = "\t") # Differential analysis file
    write.table(normalizedCounts, file = "DESeq2_normalized_counts.txt", sep = "\t") # Normalized counts file

