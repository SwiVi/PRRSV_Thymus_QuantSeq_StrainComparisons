#!/bin/bash
#SBATCH --job-name=HT-seq_count   			# name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting >
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:min>
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to ou>
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard e>
#SBATCH --mem=60G                                       # memory e.g.: 100G ; 250G ; 100M etc..>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyronchang2@gmail.com		# enter your email>

#This script is for unsorted bam files sorted by samtools based on position.
#Load the modules samtools and ht-seq count
module load htseq/2.0.2

# Define all the directories

GTF_file=/project/nadc_prrsv/Wiarda/GeneAnnotationFiles/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.gtf

mkdir -p /project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/HT_seq_output
output_dir=/project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/HT_seq_output

sorted_bam_dir=/project/nadc_prrsv/HTseqCount_SV671_Thymus_TC/Sorted_BAM

#extract the basename and define path for ht-seq output files using for loop

for sorted_bam_file in "$sorted_bam_dir"/*_sorted.bam; do
Basename=$(basename "$sorted_bam_file" _sorted.bam) #extract the name before .bam
htseq_output="${output_dir}/${Basename}.txt" #put the basename after the read count files

# perform Ht-seq count command
htseq-count -f bam -r pos -t gene -s no -m intersection-nonempty -i gene_id "$sorted_bam_file" "$GTF_file">"$htseq_output"
#Flag information
# -f format of the file
# -r pos align by position
# -s no no strand 
# -m method
# -i gene_id in gtf column; alternatively you can do gene_name, but it will only you ensemble ID.


#ht-seq output files (e.g..txt files)

done
 
#END of the line
