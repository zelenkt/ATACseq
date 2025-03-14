#!/usr/bin/Rscript

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
##### 2019.12.16 EYH #####

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


##### GET ARGUMENTS AND SET VARIABLES #####

spec = matrix(c(
    'condition' , 'c', 1, "character", "Test condition (default: 'genotype')",
    'baseline'  , 'b', 1, "character", "Control group identifier (default: 'WT')",
    'samples'       , 's', 1, "character", "ATAC-seq / ChIP-seq sample list file (default: 'test-samples.csv')",
    'readcounts'       , 'r', 1, "character", "Raw read counts file (default: '../summary/PEPATAC_tomas_mm10_peaks_coverage.tsv')"
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# OPTIONAL ARGS
if ( is.null(opt$condition ) ) { opt$condition = "genotype"                                             }
if ( is.null(opt$baseline  ) ) { opt$baseline  = "WT"                                                   }
if ( is.null(opt$samples       ) ) { opt$samples       = paste0(path.expand("/share/lab_avram/HPC_Cluster/user/tomas/pepatac/results1_macs_all/DESeq2/test-samples.csv"))        }
if ( is.null(opt$readcounts       ) ) { opt$readcounts       = paste0(path.expand("pepatac_raw_read_counts.txt"))        }


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

