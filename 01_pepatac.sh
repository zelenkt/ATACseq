#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=15GB
#SBATCH -t 48:00:00
#SBATCH --job-name=PEPATAC
#SBATCH --output=PEPATAC.%j.out
#SBATCH --error=PEPATAC.%j.err

##### LAST UPDATED 01.23.2023 Tomas Zelenka  #####
##### OPTIONAL ARGUMENTS #####

# move to the pepatac directory and then run the script as:
cd /share/lab_avram/HPC_Cluster/user/tomas/06_ATACseq_murine_tumors_Leo_Tomas_02-02-2023
# sbatch /share/lab_avram/HPC_Cluster/user/tomas/bin/atacseq_pipeline/01_pepatac.sh

ml purge
ml preseq/2.0
ml skewer/0.2.2
ml wigToBigWig/wigToBigWig
ml bedToBigBed/bedToBigBed
ml Bowtie2/2.3.5.1-GCC-8.3.0
ml samblaster/0.1.26-GCC-10.2.0
ml BEDTools/2.30.0-GCC-11.2.0
ml pigz/2.6-GCCcore-11.2.0
ml SAMtools/1.15.1-GCC-11.3.0
ml FastQC/0.11.9-Java-1.8
ml Trimmomatic/0.39-Java-1.8
ml picard/2.25.5-Java-13
ml Genrich/0.6.1
#ml Homer/4.11.1-foss-2021a probably I need to find a different position in the loading list as this makes some other packages disfunctional
ml R/4.2.1-foss-2022a

source /share/lab_avram/HPC_Cluster/user/tomas/bin/soft_versions/peptac/bin/activate

export R_LIBS_USER='/home/4476125/R'
export REFGENIE="/share/lab_avram/HPC_Cluster/annotations_genomes/genome_config.yaml"
export PYTHONPATH="/share/lab_avram/HPC_Cluster/user/tomas/bin/soft_versions/peptac/lib/python3.9/site-packages" 
export DIVCFG="/share/lab_avram/HPC_Cluster/user/tomas/pepatac/compute_config.yaml"



# it works better to run the dry run wiht -d which updates the sub scripts and then to run the sub script
#cd /share/lab_avram/HPC_Cluster/user/tomas/pepatac
# looper run tools/configuration_C8T.yaml -d
# looper runp tools/configuration_C8T.yaml -d
# looper report tools/configuration_C8T.yaml



# add the following in the sub script
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10GB
#SBATCH -t 48:00:00
#SBATCH --job-name=8ATW
#SBATCH --output=PEPATAC_mC8ATW.%j.out
#SBATCH --error=PEPATAC_mC8ATW.%j.err


# then run
# sbatch results/CD8_cells/tumor/submission/PEPATAC_C8ATK2.sub