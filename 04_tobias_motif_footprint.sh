#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=42
#SBATCH --mem-per-cpu=10GB
#SBATCH --qos=medium
#SBATCH -t 24:00:00
#SBATCH --job-name=Tobias
#SBATCH --output=Tobias.%j.out
#SBATCH --error=Tobias.%j.err

##### LAST UPDATED 07-09-2024 Tomas Zelenka  #####
##### OPTIONAL ARGUMENTS #####

# https://github.com/loosolab/TOBIAS/
# https://github.com/loosolab/TOBIAS_snakemake


# try to make my own fresh installation with conda

# ml Anaconda3/2024.02-1

# conda create --prefix /home/4476125/bin/conda_env/tobias_tz2 python=3.7
# conda install tobias -c bioconda

# # # this original version is likely not compatible with uropa or some other packages
# # # create conda environment in my home directory
# # conda create --prefix /home/4476125/bin/conda_env/tobias python=3.6

# # to remove unwanted env
# conda remove -p /home/4476125/tobias_tz --all

# # to move environments first export the packages, then move it and and re-create using the yaml file
# conda env export -p /home/4476125/kallisto/kallisto > kallisto.yaml
# mv kallisto* /home/4476125/bin/conda_env/
# conda env create -p /home/4476125/bin/conda_env/kallisto -f kallisto.yaml
# conda activate --prefix /home/4476125/bin/conda_env/kallisto
# source activate /home/4476125/bin/conda_env/kallisto

# conda env export -p /home/4476125/scvelo > scvelo.yaml
# mv scvelo* /home/4476125/bin/conda_env/
# conda env create -p /home/4476125/bin/conda_env/scvelo -f scvelo.yaml
# source activate /home/4476125/bin/conda_env/scvelo






# to run
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/atacseq_pipeline/04_tobias_motif_footprint.sh

# comparing isolated SEs
cd /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_01-07-2025/tobias2.cat1/
source activate /share/lab_avram/HPC_Cluster/user/tomas/bin/conda_env/tobias_tz2


name="sum.superenhancers.Filt.tobias1.cat3"
# name="sum.superenhancers.Filt.tobias1.cat1"
# name="sum.superenhancers.Filt.tobias2.cat1"


motifs="/share/lab_avram/HPC_Cluster/user/tomas/lib/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
genome="/share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/fasta/default/mm10.fa"

awk 'NR>1 {print $1 "\t" $2 "\t" $3}' ${name}.txt | tr -d '"' | awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' | sort -k1,1 -k2,2n > ${name}.bed

# maybe better approach I found here https://github.com/loosolab/TOBIAS/wiki/test-data:
uropa --bed ${name}.bed --gtf /share/lab_avram/HPC_Cluster/user/tomas/lib/refgenie_ensembl_mm10.gtf --show_attributes gene_id gene_name --feature_anchor start --distance 20000 10000 --feature gene
cut -f 1-6,16-17 ${name}_finalhits.txt | head -n 1 > ${name}_annotated_header.txt
cut -f 1-6,16-17 ${name}_finalhits.txt | tail -n +2 > ${name}_annotated.bed

TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATW/aligned_mm10/mC8ATW_sort_dedup.bam --genome $genome --peaks ${name}_annotated.bed --cores 40

TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATK/aligned_mm10/mC8ATK_sort_dedup.bam --genome $genome --peaks ${name}_annotated.bed --cores 40

TOBIAS FootprintScores --signal mC8ATW_sort_dedup_corrected.bw --regions ${name}_annotated.bed --output ${name}_mC8ATW_footprints.bw --cores 40
TOBIAS FootprintScores --signal mC8ATK_sort_dedup_corrected.bw --regions ${name}_annotated.bed --output ${name}_mC8ATK_footprints.bw --cores 40

TOBIAS BINDetect --motifs $motifs --signals ${name}_mC8ATK_footprints.bw ${name}_mC8ATW_footprints.bw --genome $genome --peaks ${name}_annotated.bed  --peak_header ${name}_annotated_header.txt --outdir BINDetect_output_annot_${name} --cond_names mC8ATK mC8ATW --cores 40




# prepare motif2gene mapping file
grep -i '>' JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | tr -d '>' > jaspar_filt.txt

# rearrange the jaspar file in excel which is faster - to remove the occurrences with two TF on the same line - duplicate them on two lines...then make first column lowercase with this command
awk '{print tolower($1) "\t" $0}' jaspar_filt_rearranged.txt > jaspar_filt_rearranged2.txt
# take the annotation file from pipseeker annotation that was used in scRNAseq analysis in ~/bin/pipseeker-v2.1.4-linux/pipseeker-gex-reference-GRCm39-2022.04
awk '{print tolower($2) "\t" $2 "\t" $1}' pipseeker-gex-reference-GRCm39-2022.04_geneInfo.tab > pipseeker-gex-reference-GRCm39-2022.04_geneInfo_modif.tab

awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' pipseeker-gex-reference-GRCm39-2022.04_geneInfo_modif.tab jaspar_filt_rearranged2.txt | grep -v -i "#NA" | awk -F'\t' -v OFS="\t" '{print $8,$2,$3,$4,$7}' > jaspar_geneID_annotated.txt


awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024.rna_allClustTogether.txt jaspar_geneID_annotated.txt > jaspar_annotated_by_scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024.txt
# new columns are geneID    gene  pct.1-pct.2  avg_log2FC  pct.1   pct.2  pval   adjPval    description
# this file was further modified in excel

name="sum.superenhancers.Filt.tobias1.cat1"
awk -F'\t' -v OFS="\t" '{print $3,$0}' tobias1.cat1/BINDetect_output_annot_${name}/bindetect_results.txt > tobias_res_${name}.txt
# cols in results: output_prefix	name	motif_id	cluster	total_tfbs	mC8ATK_mean_score	mC8ATK_bound	mC8ATW_mean_score	mC8ATW_bound	mC8ATK_mC8ATW_change	mC8ATK_mC8ATW_pvalue


# NOTE this order will only consider one of the pairs so for example for this Ahr::Arnt it would only take Arnt and not Ahr so I need to reverse the order and plase tobias_res first
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' jaspar_annotated_by_scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024_modif.txt tobias_res_${name}.txt > final_tobias_res_RNAannot_${name}.txt

awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tobias_res_${name}.txt jaspar_annotated_by_scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024_modif.txt > final_reversed_tobias_res_RNAannot_${name}.txt

# I have to add the following heading at the end of existing heading: motifID	motifGeneID	singleGeneMotif	geneMotif	refGene	RNAid	RNAgene	pct.1-pct.2	avg_log2FC	pct.1	pct.2	pval	adjPval	description


# first I'm testing just the whole SE regions...in the next step I can only intersect with Bcl11b peaks? though that would probably only show bcl11b motifs...what about differential K27ac peaks intersecting with SE or TE


TOBIAS PlotAggregate --TFBS tobias1.cat3/BINDetect_output_annot_sum.superenhancers.Filt.tobias1.cat3/BACH2_MA1470.2/beds/BACH2_MA1470.2_all.bed  --signals ../for_motif_tobias_7-12-2024/mC8ATW_sort_dedup_corrected.bw ../for_motif_tobias_7-12-2024/mC8ATK_sort_dedup_corrected.bw --output BACH2_MA1470_2_footprint_comparison_all_atacseq-all.pdf --share_y both --plot_boundaries --signal-on-x


name="sum.superenhancers.Filt.tobias1.cat3"
TOBIAS FilterFragments --bw ../for_motif_tobias_7-12-2024/mC8ATW_sort_dedup_corrected.bw --regions test_data/merged_peaks.bed

TOBIAS PlotAggregate --TFBS tobias1.cat3/BINDetect_output_annot_sum.superenhancers.Filt.tobias1.cat3/BACH2_MA1470.2/beds/BACH2_MA1470.2_all.bed  --signals ../for_motif_tobias_7-12-2024/mC8ATW_sort_dedup_corrected.bw ../for_motif_tobias_7-12-2024/mC8ATK_sort_dedup_corrected.bw --output BACH2_MA1470_2_footprint_comparison_all_tobias1.cat3.pdf --share_y both --plot_boundaries --signal-on-x

cd /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_01-07-2025/tobias1.cat1/
TOBIAS CreateNetwork --TFBS annotated_tfbs/* --origin ../motif2gene_mapping.txt 
# still don't have the right motif2gene mapping file and will not create it now...also I would probably have to modify the original JASPAR FILE so that not motifs are in duplicates which is perhaps the reason why it added the motif ID to the gene name Alx1_MA0854.2 so this would be a problem and it does not correspond to the example data in /share/lab_avram/HPC_Cluster/user/tomas/bin/TOBIAS_snakemake/data/annotated_tfbs
# hear I started extracting genes and motifs but it's wrong for i in *; do cd $i; awk '{print $4 "\t" $13}' beds/${i}_mC8ATK_bound.bed | sort | uniq | head ; cd ../; done | head






name="BINDetect_output_annot"
awk -F'\t' -v OFS="\t" '{print $3,$0}' BINDetect_output_annot/bindetect_results.txt > tobias_res_${name}.txt
# cols in results: output_prefix	name	motif_id	cluster	total_tfbs	mC8ATK_mean_score	mC8ATK_bound	mC8ATW_mean_score	mC8ATW_bound	mC8ATK_mC8ATW_change	mC8ATK_mC8ATW_pvalue

# # NOTE this order will only consider one of the pairs so for example for this Ahr::Arnt it would only take Arnt and not Ahr so I need to reverse the order and plase tobias_res first
# awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' jaspar_annotated_by_scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024_modif.txt tobias_res_${name}.txt > final_tobias_res_RNAannot_${name}.txt
awk -F'\t' -v OFS="\t" 'FNR==NR{a[$1]=$0; next} { if ($1 in a) {print $0,a[$1]} else {print $0,"#NA"} }' tobias_res_${name}.txt jaspar_annotated_by_scRNAseq_adapted_markers.seurat_obj_reclust2_12-17-2024_modif.txt > final_reversed_tobias_res_RNAannot_${name}.txt


echo 'SP8 Foxn1' > TF_names.txt
TOBIAS PlotChanges --bindetect BINDetect_output_annot/bindetect_results.txt --TFS TF_names.txt


TOBIAS ClusterMotifs --motifs $motifs --threshold 0.4 --type png




 

# # original working runs
# cd /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_7-12-2024/
# source activate /home/4476125/bin/conda_env/tobias_tz2

# genome="/share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/fasta/default/mm10.fa"
# peaks="/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_7-12-2024/merged_all_atac_C8AT_peaks.bed"
# # TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATW/aligned_mm10/mC8ATW_sort_dedup.bam --genome $genome --peaks $peaks
# # TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATK/aligned_mm10/mC8ATK_sort_dedup.bam --genome $genome --peaks $peaks

# TOBIAS FootprintScores --signal mC8ATW_sort_dedup_corrected.bw --regions $peaks --output mC8ATW_sort_dedup_footprints.bw --cores 40
# TOBIAS FootprintScores --signal mC8ATK_sort_dedup_corrected.bw --regions $peaks --output mC8ATK_sort_dedup_footprints.bw --cores 40

# I downloaded motifs from here..JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt, vertebrates single batch jaspar file https://jaspar.elixir.no/downloads/  https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
  


# # maybe better approach I found here https://github.com/loosolab/TOBIAS/wiki/test-data:
# uropa --bed merged_peaks.bed --gtf transcripts_chr4.gtf --show_attributes gene_id gene_name --feature_anchor start --distance 20000 10000 --feature gene 
# cut -f 1-6,16-17 merged_peaks_finalhits.txt | head -n 1 > merged_peaks_annotated_header.txt
# cut -f 1-6,16-17 merged_peaks_finalhits.txt | tail -n +2 > merged_peaks_annotated.bed


# # # this is probably not needed but they use annotated peak file....to generate annotated merged-peaks file I use UROPA which is an annotation tool developed by the same group, modify the json file to modify the name of the input file or gtf
# uropa -i merged_peaks_uropa_config.json --prefix uropa_results -threads 40 -l uropa.log
# cut -f 1-4,7-13,16-19 uropa_results_finalhits.txt | awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,".",$8,$13,$14}' > uropa_results_finalhits_sub.txt #Get a subset of columns
# head -n 1 uropa_results_finalhits_sub.txt > all_merged_annotated_header.txt #header
# tail -n +2 uropa_results_finalhits_sub.txt > all_merged_annotated.bed  #bedlines
# rm uropa_results_finalhits_sub.txt
# peak_chr	peak_start	peak_end	peak_id	peak_score	peak_strand	gene_id	gene_name
# awk -F'\t' -v OFS="\t" '{print $1,$2,$3,$4,".",$8,$13,$14}' all_merged_annotated.bed | head


# # the peak file must not have any obscure chromosomes so filter them out
# # # awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrM/ || $1~/^chrY/) {print $1}' all_merged_annotated_backup.bed | sort | uniq
# awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' all_merged_annotated_backup.bed | sort -k1,1 -k2,2n > all_merged_annotated_filt.bed 
# awk '!($1 ~ /^chrUn/ || $1 ~ /random$/ || $1~/^chrM/ || $1~/^chrY/) {print $0}' merged_all_atac_C8AT_peaks.bed | sort -k1,1 -k2,2n > merged_all_atac_C8AT_peaks_filt.bed 


# # for Testing
# cd /share/lab_avram/HPC_Cluster/user/tomas/bin/TOBIAS_snakemake/data/
# TOBIAS BINDetect --motifs motifs.jaspar --signals Bcell_footprints.bw Tcell_footprints.bw --genome genome.fa.gz --peaks merged_peaks_annotated.bed --peak_header merged_peaks_annotated_header.txt --outdir BINDetect_output_TEST --cond_names Bcell Tcell --cores 40 --verbosity 5

# # THIS WAS USED
# TOBIAS BINDetect --motifs /share/lab_avram/HPC_Cluster/user/tomas/lib/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals mC8ATK_sort_dedup_footprints.bw mC8ATW_sort_dedup_footprints.bw --genome $genome --peaks merged_all_atac_C8AT_peaks_filt.bed  --peak_header merged_all_atac_C8AT_peaks_header.txt --outdir BINDetect_output_nonAnnotPeaks --cond_names mC8ATK mC8ATW --cores 40

# # # this is an option with annotated peaks but there are no major differences in the results so I wonder if it is necessary to input annotated peaks, perhaps not
# TOBIAS BINDetect --motifs /share/lab_avram/HPC_Cluster/user/tomas/lib/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt --signals mC8ATK_sort_dedup_footprints.bw mC8ATW_sort_dedup_footprints.bw --genome $genome --peaks all_merged_annotated_filt.bed  --peak_header all_merged_annotated_header.txt --outdir BINDetect_output_annot --cond_names mC8ATK mC8ATW --cores 40


# TOBIAS PlotHeatmap --TFBS BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_mC8ATW_bound.bed BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_mC8ATW_unbound.bed BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_mC8ATK_bound.bed BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_mC8ATK_unbound.bed --signals mC8ATW_sort_dedup_corrected.bw mC8ATK_sort_dedup_corrected.bw --signal_labels WT KO --share_colorbar --output heatmap_bcl11b.png --sort_by -1

# TOBIAS PlotHeatmap --TFBS BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_all.bed --signals mC8ATW_sort_dedup_corrected.bw mC8ATK_sort_dedup_corrected.bw --output heatmap2_bcl11b.png --sort_by -2

# TOBIAS PlotAggregate --TFBS BINDetect_output_nonAnnotPeaks/Bcl11B_MA1989.2/beds/Bcl11B_MA1989.2_all.bed  --signals mC8ATW_sort_dedup_corrected.bw mC8ATK_sort_dedup_corrected.bw --output aggregate_bcl11b.png --share_y both --plot_boundaries --signal-on-x



# CONTINUE

# # testing network
# cd /share/lab_avram/HPC_Cluster/user/tomas/bin/TOBIAS_snakemake/data/
# TOBIAS CreateNetwork --TFBS annotated_tfbs/* --origin motif2gene_mapping.txt 
# it works but to make graphical plots I need to import it to R - I played with basics here /Users/4476125/Documents/datasets/08_rnaseq_atacseq_cutrun/data/results/tobias_networks ..I will also want to add RNAseq information and overlay it as color code
# see Fig. 2G as an example how they used the networks from TOBIAS: https://www.nature.com/articles/s41467-024-53295-1 (G) A TF-TF network is built of RUNX1 and the other top 5% TF genes in RT+ fibroblasts. Sizes of nodes represent the level of the network starting with RUNX1. Directed edges indicate TF binding sites in the respective gene.

# think of the best way how to create the motif2gene mapping...whether using all genes or just the genes associated with my peaks...ideally all genes ..then probably extract it from  /share/lab_avram/HPC_Cluster/user/tomas/lib/refgenie_ensembl_mm10.gtf

# also create this heatmap and re-do the volcano with labeling of my interest https://www.nature.com/articles/s41467-020-18035-1#Fig2


# for real run I need to generate motif2gene mapping file and then group all the bed files from annotated BINDetect run...also consider rerunning the peak annotation using the most recent UROPA approach that I found above
# # the creation of network requires output from BINDetect that is gene-annotated - thus the original peaks need to be gene-annotated using the UROPA tool above
# TOBIAS CreateNetwork --TFBS test_data/annotated_tfbs/* --origin test_data/motif2gene_mapping.txt 









# https://github.com/loosolab/TOBIAS/wiki/PlotTracks this allows to show directly genomic tracks


# Bcl11B_MA1989.2


# KLF4_MA0039.4
# NFATC4_MA1525.1
# NR4A2_MA0160.1

# NFKB1_MA0105.4
# NFKB2_MA0778.1



# TOBIAS PlotAggregate --TFBS test_data/BATF_all.bed  --signals test_data/Bcell_corrected.bw test_data/Tcell_corrected.bw --output BATFJUN_footprint_comparison_all.pdf --share_y both --plot_boundaries --signal-on-x




# install also snakemake pipeline https://github.com/loosolab/TOBIAS_snakemake
# conda env create -f environments/snakemake.yaml
# conda activate tobias_snakemake_env

# /home/4476125/bin/TOBIAS_snakemake/data





# # cd /Volumes/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_7-12-2024 
# # cat mC8ATW_peaks_rmBlacklist.narrowPeak mC8ATK_peaks_rmBlacklist.narrowPeak | bedtools sort | bedtools merge > merged_all_atac_C8AT_peaks.bed

# # # I have problems with python version...these tricks also does not work
# # srun --export=SBATCH_EXPORT=NONE ~/bin/tobias_motif_footprint.sh
# # sbatch --export ALL ~/bin/tobias_motif_footprint.sh



# # https://www.hdsu.org/chipatac2020/07_ATAC_Footprinting.html

# ml purge
# #module load Python/3.8.6-GCCcore-10.2.0
# # ml Python/3.6.6-foss-2018b 
# # matplotlib/3.0.3-foss-2019a-Python-3.7.2
# python --version
# which python

# source /home/4476125/tobias/bin/activate

# export PYTHONPATH="/home/4476125/tobias/lib/python3.6/site-packages" 

# ml Python/3.6.6-foss-2018b 

# python --version
# which python


# cd /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_7-12-2024/

# genome="/share/lab_avram/HPC_Cluster/annotations_genomes/alias/mm10/fasta/default/mm10.fa "
# peaks="/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/for_motif_tobias_7-12-2024/merged_all_atac_C8AT_peaks.bed"

# TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATW/aligned_mm10/mC8ATW_sort_dedup.bam --genome $genome --peaks $peaks

# # TOBIAS ATACorrect --bam /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/merged/results_pipeline/mC8ATK/aligned_mm10/mC8ATK_sort_dedup.bam --genome $genome --peaks $peaks


date