#!/bin/bash
#SBATCH --job-name=sort_bam_files                       # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting >
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:min>
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to ou>
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard e>
#SBATCH --mem=32G                                       # memory e.g.: 100G ; 250G ; 100M etc..>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyronchang2@gmail.com		# enter your email>

#This script is for unsorted bam files sorted by samtools based on position.
#Load the modules samtools and ht-seq count
module load samtools/1.17

#Define the directories.

STAR_BAM_dir=/project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/Alignments_STAR

mkdir -p /project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/Sorted_BAM
sorted_bam_file_dir=/project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/Sorted_BAM   

#Use a for loop to execute the script to sort and index bam files.
# first loop through the bam files

for bamfiles in "$STAR_BAM_dir"/*.bam; do

Basename=$(basename "$bamfiles" .bam) #extract the name before .bam

sorted_bam_file="${sorted_bam_file_dir}/${Basename}_sorted.bam"

samtools sort "$bamfiles" -o "$sorted_bam_file"

samtools index "$sorted_bam_file"

done
#End of the line
