#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=8:00:00
#SBATCH --job-name=DESeq2
#SBATCH --output=DESEq2.%j.out
#SBATCH --error=DESEq2.%j.err

##### LAST UPDATED 01/24/2023 Tomas Zelenka #####
#####  2021.11.23 Brian Connolly #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-r file] [-t string] [-c string]


Optional Arguments:
    -h  Show this help text
    -d  Set target directory (directory to write analysis files) (default: $PWD/DESeq2)
    -s  Set ATAC-seq / ChIP-seq sample list file (default: /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/test-samples.csv)
    -t  Set test condition (column from ATAC-seq / ChIP-seq sample list file) (i.e. genotype or cell_subset) (default: genotype)
    -c  Set control group (i.e. WT if test condition is genotype, or Th0 if test condition is cell_subset (default: WT)
    -r  Set Raw reads counts as an output from PEPATAC pipeline, from feature counts or from other sources

Example usage:
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/pepatac/results1_macs_all --export=s=/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/test-samples.csv ~/bin/pepatac_DESeq2.sh"


while getopts ':h:d:s:t:c:' option; do
    case "$option" in
    h)  echo "$usage"
        exit 0
        ;;
	d)  if [[ -d "$OPTARG" ]]; then
            d="$OPTARG"
            printf "Target directory manually set to %s\n" "$OPTARG"
        else
            printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
            exit 1
        fi
        ;;
    
    r)  if [[ -r "$OPTARG" ]]; then
            r="$OPTARG"
            printf "Raw read counts manually set to: %s\n" "$OPTARG"
        else
            printf "%s is empty or does not exist. Please specify a valid file with Raw read counts and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
        fi
        ;;
    s)  if [[ -s "$OPTARG" ]]; then
            s="$OPTARG"
            printf "ATAC-seq / ChIP-seq sample list file manually set to: %s\n" "$OPTARG"
        else
            printf "%s is empty or does not exist. Please specify valid ATAC-seq / ChIP-seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"
        fi
        ;;
    t)  t="$OPTARG"
        printf "Test condition manually set to: %s\n" "$OPTARG"
        ;;
    c)  c="$OPTARG"
        printf "Control group manually set to: %s\n" "$OPTARG"
        ;;
   \?)  printf "Invalid option: -%s.\n" "$OPTARG"
        echo "$usage"
        exit 1
        ;;
    :)  printf "Invalid option: -%s requires an argument.\n" "$OPTARG"
        echo "$usage"
        exit 1
        ;;
    esac
done
shift $((OPTIND -1))

date;pwd

##### LOAD MODULES #####
ml purge
module load R/4.1.0-foss-2020b

##### Use Avram Lab's shared R library #####
export R_LIBS="/share/lab_avram/HPC_Cluster/share/R/x86_64-pc-linux-gnu-library/4.1"

##### SET VARIABLES #####
target_dir="${d:-${PWD}/DESeq2}"
if [[ "$target_dir" == "$PWD/DESeq2" && ! -d "$target_dir" ]]; then # IF TARGET DIR DOES NOT EXIST
    mkdir -p "$target_dir"
elif [[ ! -d "$target_dir" ]]; then # THIS CAN ONLY BE TRUE IF SOMEONE CALLED THIS FUNCTION THROUGH SBATCH AND DIRECTORY EXPORTED AN INVALID "d" variable. (e.g. sbatch --export=d=/some/invalid/dir DESeq2.sh)
    printf "%s doesn't exist or is not a valid directory. Please specify a valid directory and rerun %s. Exiting...\n" "$target_dir" "$(basename $0)"
    exit 1
fi

readcounts="${r:-${PWD}/summary/PEPATAC_tomas_mm10_peaks_coverage.tsv}"  # samples="${r:-${PWD}/DESeq2/test-samples.csv}"
[[ -r "$readcounts" ]] || { printf "%s is empty or does not exist. Please specify  a valid file with Raw read counts and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY RAW READ COUTNS FILE IS VALID

# prepare the raw read count output from pepatac (the original TSV file needs to be converted like this otherwise it creates problems when preparing matrix)
colCount=$(awk '{print NF; exit}' ${readcounts}); for ((col=1; col<=$colCount; col++)); do awk -v col="$col" '{print $(col)}' ${readcounts} > tmp_columns_for_readcountcoverage_$col; done; paste tmp_columns_for_readcountcoverage* > $target_dir/pepatac_raw_read_counts.txt; rm tmp_columns_for_readcountcoverage*


samples="${s:-/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/test-samples.csv}"  # samples="${r:-${PWD}/DESeq2/test-samples.csv}"
[[ -s "$samples" ]] || { printf "%s is empty or does not exist. Please specify valid ATAC-seq / ChIP-seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY SAMPLE LIST FILE IS VALID

condition="${t:-genotype}"
read -ra condition_list <<<$(head -n1 "$samples") # GENERATE ARRAY OF POSSIBLE CONDITIONS
[[ "${condition_list[@]}" =~ "$condition" ]] || { printf "%s is not a valid condition.  Please select a condition from the following list: %s. Exiting...\n" "$condition" "${condition_list[*]}"; exit 1; } #VERIFY SELECTED CONDITION IS VALID

baseline="${c:-WT}"
baseline_check=$(awk -F$',' -v condition=$condition -v baseline=$baseline 'BEGIN {s=0} 
    NR==1{for(i=1; i<=NF; i++) if($i==condition) {a[i]++; break} }
    {for (i in a) if ($i==baseline) {s++; exit} }
    END {print s}' "$samples")
[[ $baseline_check > 0 ]] || { printf "There are no ATAC-seq / ChIP-seq samples of type \"%s\" under condition \"%s\". Please select a valid control group and rerun %s. Exiting...\n" "$baseline" "$condition" "$(basename $0)"; exit 1; } #VERIFY THAT AT LEAST ONE SAMPLE UNDER SELECTED CONDITION MATCHES BASELINE TYPE



##### BEGIN SCRIPT #####
printf "Condition = %s\n" "$condition"
printf "Baseline = %s\n" "$baseline"
printf "ATAC-seq / ChIP-seq Sample list file = %s\n" "$samples"
printf "Raw read counts file = %s\n" "$readcounts"

( cd "$target_dir" && Rscript "$HOME"/bin/pepatac_DESeq2.R -c "$condition" -b "$baseline" -s "$samples" -r "$readcounts" ) # RUN DESEQ2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR

# ( cd "$target_dir" && Rscript "$HOME"/bin/annotate_genes.R -i DESeq2_analysis.txt -o DESeq2_analysis_annotated.txt ) # RUN ANNOTATE_GENES_DESE2.R WITHIN SUBSHELL TO CHANGE DIRECTORY TO ANALYSIS DIR



mv *.err $target_dir/
mv *.out $target_dir/

date

