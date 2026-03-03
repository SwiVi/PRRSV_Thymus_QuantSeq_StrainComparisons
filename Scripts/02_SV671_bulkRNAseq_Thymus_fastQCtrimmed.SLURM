#!/bin/bash

#SBATCH --job-name=fastQC_trimmed                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                             # number of cores/tasks in this job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=32G   					# memory e.g.: 100G ; 250G ; 100M etc... (need mem node for some reqs)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jayne.wiarda@usda.gov					# enter your email address here USDA emails sometimes dont work, I use gmail

# ENTER COMMANDS HERE:

# Run fastQC on raw .fastq files 

# Load fastQC module:
module load fastqc/0.12.1

# Make directory for output files:
mkdir /project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/fastQC_trimmed

# Run fastQC:
for file in /project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/Trimmomatic/*Thymus*fastq.gz
do
fastqc -f fastq -t 6 --noextract --nogroup --outdir /project/nadc_prrsv/Wiarda/SV671_bulkRNAseq_Thymus/fastQC_trimmed ${file} 
done

#End of file