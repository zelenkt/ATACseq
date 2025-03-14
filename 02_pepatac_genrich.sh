#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=30GB
#SBATCH --time=8:00:00
#SBATCH --job-name=Genrich
#SBATCH --output=Genrich.%j.out
#SBATCH --error=Genrich.%j.err

##### LAST UPDATED 02/17/2023 Tomas Zelenka #####

##### ARGUMENTS AND USAGE #####
usage="$(basename $0) [-h] [-d dir] [-t string] 


Optional Arguments:
    -h  Show this help text
    -d  Set target directory (directory to write analysis files) (default: $PWD/genrich)
    -s  Set ATAC-seq / ChIP-seq sample list file (default: /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/test-samples.csv)
    -t  Set test condition (column from ATAC-seq / ChIP-seq sample list file) (i.e. genotype or cell_subset) (default: genotype)

Example usage:
    sbatch --chdir=/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/CD4_cells/peritoneal_lavage --export=s=/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/samples_C4P.csv ~/bin/pepatac_DESeq2.sh"


while getopts ':h:d:t:' option; do
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
ml SAMtools/1.15.1-GCC-11.3.0
ml BEDTools/2.30.0-GCC-11.2.0
ml Genrich/0.6.1

genome="hg38"
#genome="mm10"
target_dir="${d:-/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/melanoma_human_TILs/patient41/}"
#target_dir="${d:-/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/results/CD8_cells/tumor/}"
mkdir ${target_dir}genrich
cd ${target_dir}genrich

#samples_C8T.csv
#samples_TILs_p19.csv
samples="${s:-/share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023/tools/samples_TILs_p41.csv}"  # samples="${r:-${PWD}/DESeq2/test-samples.csv}"
[[ -s "$samples" ]] || { printf "%s is empty or does not exist. Please specify valid ATAC-seq / ChIP-seq sample list file and rerun %s. Exiting...\n" "$OPTARG" "$(basename $0)"; exit 1; } # VERIFY SAMPLE LIST FILE IS VALID

condition="${t:-genotype}"
read -ra condition_list <<<$(head -n1 "$samples") # GENERATE ARRAY OF POSSIBLE CONDITIONS
[[ "${condition_list[@]}" =~ "$condition" ]] || { printf "%s is not a valid condition.  Please select a condition from the following list: %s. Exiting...\n" "$condition" "${condition_list[*]}"; exit 1; } #VERIFY SELECTED CONDITION IS VALID

awk -F',' -v condition=$condition 'BEGIN {s=0} 
    NR==1{for(i=1; i<=NF; i++) if($i==condition) {a[i]++; break} }
    NR>1{print $i}' $samples | sort | uniq > tmp_list_uniq_condition.txt
conditionCount=$(cat tmp_list_uniq_condition.txt | wc -l)

for ((conditionOrder=1; conditionOrder<=$conditionCount; conditionOrder++))
    do conditionName=$(sed -n ${conditionOrder}p ${target_dir}genrich/tmp_list_uniq_condition.txt)
        touch tmp_list_of_namesort_dedup_pwd.txt > tmp_list_of_namesort_dedup_pwd.txt

        awk -F',' -v condition=$condition -v conditionName=$conditionName 'BEGIN {s=0} 
            NR==1{for(i=1; i<=NF; i++) if($i==condition) {a[i]++; break} }
            {if ($i==conditionName) print $1}' $samples | sort > tmp_counting_samples.txt
        sampleCount=$(cat tmp_counting_samples.txt | wc -l)

        for ((sampleOrder=1; sampleOrder<=$sampleCount; sampleOrder++))
            do sampleName=$(sed -n ${sampleOrder}p tmp_counting_samples.txt)
            
            if [[ ! -f ${target_dir}results_pipeline/${sampleName}/aligned_${genome}/${sampleName}_namesort_dedup.bam ]]; then
                printf "Sorting ${sampleName} SAM/BAM file by queryname...\n"

                samtools sort -n ${target_dir}results_pipeline/${sampleName}/aligned_${genome}/${sampleName}_sort_dedup.bam > ${target_dir}results_pipeline/${sampleName}/aligned_${genome}/${sampleName}_namesort_dedup.bam
                printf "Done.\n"
            else
                printf "Sorted ${sampleName} SAM/BAM file by queryname already exists\n" "$base"
            fi 
            echo "${target_dir}results_pipeline/${sampleName}/aligned_${genome}/${sampleName}_namesort_dedup.bam" >> tmp_list_of_namesort_dedup_pwd.txt
            done; rm tmp_counting_samples.txt

        NameSortBamList=$(sed -z 's/\n/,/g;s/,$/ /' tmp_list_of_namesort_dedup_pwd.txt | cat)

        Genrich -j -t $NameSortBamList -E /share/lab_avram/HPC_Cluster/annotations_genomes/alias/${genome}/blacklist/default/${genome}_blacklist.bed.gz -e chrY,chrM,chr1_GL456210_random,chr1_GL456211_random,chr1_GL456212_random,chr1_GL456213_random,chr1_GL456221_random,chr4_GL456216_random,chr4_JH584292_random,chr4_GL456350_random,chr4_JH584293_random,chr4_JH584294_random,chr4_JH584295_random,chr5_JH584296_random,chr5_JH584297_random,chr5_JH584298_random,chr5_GL456354_random,chr5_JH584299_random,chr7_GL456219_random,chrX_GL456233_random,chrY_JH584300_random,chrY_JH584301_random,chrY_JH584302_random,chrY_JH584303_random,chrUn_GL456239,chrUn_GL456367,chrUn_GL456378,chrUn_GL456381,chrUn_GL456382,chrUn_GL456383,chrUn_GL456385,chrUn_GL456390,chrUn_GL456392,chrUn_GL456393,chrUn_GL456394,chrUn_GL456359,chrUn_GL456360,chrUn_GL456396,chrUn_GL456372,chrUn_GL456387,chrUn_GL456389,chrUn_GL456370,chrUn_GL456379,chrUn_GL456366,chrUn_GL456368,chrUn_JH584304,chr1_KI270706v1_random,chr1_KI270707v1_random,chr1_KI270708v1_random,chr1_KI270709v1_random,chr1_KI270710v1_random,chr1_KI270711v1_random,chr1_KI270712v1_random,chr1_KI270713v1_random,chr1_KI270714v1_random,chr2_KI270715v1_random,chr2_KI270716v1_random,chr3_GL000221v1_random,chr4_GL000008v2_random,chr5_GL000208v1_random,chr9_KI270717v1_random,chr9_KI270718v1_random,chr9_KI270719v1_random,chr9_KI270720v1_random,chr11_KI270721v1_random,chr14_GL000009v2_random,chr14_GL000225v1_random,chr14_KI270722v1_random,chr14_GL000194v1_random,chr14_KI270723v1_random,chr14_KI270724v1_random,chr14_KI270725v1_random,chr14_KI270726v1_random,chr15_KI270727v1_random,chr16_KI270728v1_random,chr17_GL000205v2_random,chr17_KI270729v1_random,chr17_KI270730v1_random,chr22_KI270731v1_random,chr22_KI270732v1_random,chr22_KI270733v1_random,chr22_KI270734v1_random,chr22_KI270735v1_random,chr22_KI270736v1_random,chr22_KI270737v1_random,chr22_KI270738v1_random,chr22_KI270739v1_random,chrY_KI270740v1_random,chrUn_KI270302v1,chrUn_KI270304v1,chrUn_KI270303v1,chrUn_KI270305v1,chrUn_KI270322v1,chrUn_KI270320v1,chrUn_KI270310v1,chrUn_KI270316v1,chrUn_KI270315v1,chrUn_KI270312v1,chrUn_KI270311v1,chrUn_KI270317v1,chrUn_KI270412v1,chrUn_KI270411v1,chrUn_KI270414v1,chrUn_KI270419v1,chrUn_KI270418v1,chrUn_KI270420v1,chrUn_KI270424v1,chrUn_KI270417v1,chrUn_KI270422v1,chrUn_KI270423v1,chrUn_KI270425v1,chrUn_KI270429v1,chrUn_KI270442v1,chrUn_KI270466v1,chrUn_KI270465v1,chrUn_KI270467v1,chrUn_KI270435v1,chrUn_KI270438v1,chrUn_KI270468v1,chrUn_KI270510v1,chrUn_KI270509v1,chrUn_KI270518v1,chrUn_KI270508v1,chrUn_KI270516v1,chrUn_KI270512v1,chrUn_KI270519v1,chrUn_KI270522v1,chrUn_KI270511v1,chrUn_KI270515v1,chrUn_KI270507v1,chrUn_KI270517v1,chrUn_KI270529v1,chrUn_KI270528v1,chrUn_KI270530v1,chrUn_KI270539v1,chrUn_KI270538v1,chrUn_KI270544v1,chrUn_KI270548v1,chrUn_KI270583v1,chrUn_KI270587v1,chrUn_KI270580v1,chrUn_KI270581v1,chrUn_KI270579v1,chrUn_KI270589v1,chrUn_KI270590v1,chrUn_KI270584v1,chrUn_KI270582v1,chrUn_KI270588v1,chrUn_KI270593v1,chrUn_KI270591v1,chrUn_KI270330v1,chrUn_KI270329v1,chrUn_KI270334v1,chrUn_KI270333v1,chrUn_KI270335v1,chrUn_KI270338v1,chrUn_KI270340v1,chrUn_KI270336v1,chrUn_KI270337v1,chrUn_KI270363v1,chrUn_KI270364v1,chrUn_KI270362v1,chrUn_KI270366v1,chrUn_KI270378v1,chrUn_KI270379v1,chrUn_KI270389v1,chrUn_KI270390v1,chrUn_KI270387v1,chrUn_KI270395v1,chrUn_KI270396v1,chrUn_KI270388v1,chrUn_KI270394v1,chrUn_KI270386v1,chrUn_KI270391v1,chrUn_KI270383v1,chrUn_KI270393v1,chrUn_KI270384v1,chrUn_KI270392v1,chrUn_KI270381v1,chrUn_KI270385v1,chrUn_KI270382v1,chrUn_KI270376v1,chrUn_KI270374v1,chrUn_KI270372v1,chrUn_KI270373v1,chrUn_KI270375v1,chrUn_KI270371v1,chrUn_KI270448v1,chrUn_KI270521v1,chrUn_GL000195v1,chrUn_GL000219v1,chrUn_GL000220v1,chrUn_GL000224v1,chrUn_KI270741v1,chrUn_GL000226v1,chrUn_GL000213v1,chrUn_KI270743v1,chrUn_KI270744v1,chrUn_KI270745v1,chrUn_KI270746v1,chrUn_KI270747v1,chrUn_KI270748v1,chrUn_KI270749v1,chrUn_KI270750v1,chrUn_KI270751v1,chrUn_KI270752v1,chrUn_KI270753v1,chrUn_KI270754v1,chrUn_KI270755v1,chrUn_KI270756v1,chrUn_KI270757v1,chrUn_GL000214v1,chrUn_KI270742v1,chrUn_GL000216v2,chrUn_GL000218v1,chrEBV -o ${target_dir}genrich/${conditionName}_genrich_peaks.narrowPeak

        rm tmp_list_of_namesort_dedup_pwd.txt

    done


cat *_genrich_peaks.narrowPeak | sortBed | mergeBed -c 4,5,6,7,8,9 -o collapse,mean,distinct,mean,mean,mean | sortBed | uniq > combined_genrich_peaks.narrowPeak
rm tmp_list_uniq_condition.txt

###### if I want to you bedtools intersect, I can use the following condition:
#peakFileCount=$(ls *_genrich_peaks.narrowPeak | wc -l)
#if [ $peakFileCount ==  ]; then
#printf "Combining two peak datasets together...\n"
#bedtools command
#printf "Done.\n"
#else
#printf "Please make sure only two conditions are used \n" "$base"
#fi 

#peakFileCount=$(awk 'END { print NR }' tmp_list_uniq_condition.txt)

#file1=$(awk 'NR==1{print $1}' tmp_list_uniq_condition.txt)"_genrich_peaks.narrowPeak"
#file2=$(awk 'NR==2{print $1}' tmp_list_uniq_condition.txt)"_genrich_peaks.narrowPeak"



mv ../*.err $target_dir/genrich
mv ../*.out $target_dir/genrich

date;pwd


