#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J run_slurm_job_processing_atac_multicov
#SBATCH -o run_slurm_job_processing_atac_multicov_%A.o
#SBATCH -e run_slurm_job_processing_atac_multicov_%A.e
#SBATCH --ntasks 1
#SBATCH --time 48:00:00
#SBATCH --mem=36G

###################################
start_time=$(date +"%A, %F %T %Z")
current_user=$(whoami)
hostname=$(hostname)
job_id=$SLURM_JOB_ID
job_name=$SLURM_JOB_NAME
###################################


module load bbc2/deeptools/deeptools-3.5.2
module load bbc2/bedtools/bedtools-2.30.0
module load bbc2/ucsc_tools/ucsc_tools-20231127
# Work directory
workdir="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39"
cd ${workdir}


# Read counts (bed file is from the diffbind)
# mkdir -p ${workdir}/count
# bedtools multicov \
# -bams ${workdir}/align/Ctrl-WT1.sorted.rmdup.bam ${workdir}/align/Ctrl-WT2.sorted.rmdup.bam ${workdir}/align/Ctrl-KO1.sorted.rmdup.bam ${workdir}/align/Ctrl-KO2.sorted.rmdup.bam ${workdir}/align/R16-WT1.sorted.rmdup.bam ${workdir}/align/R16-WT2.sorted.rmdup.bam ${workdir}/align/R16-KO1.sorted.rmdup.bam ${workdir}/align/R16-KO2.sorted.rmdup.bam ${workdir}/align/R36-WT1.sorted.rmdup.bam ${workdir}/align/R36-WT2.sorted.rmdup.bam ${workdir}/align/R36-KO1.sorted.rmdup.bam ${workdir}/align/R36-KO2.sorted.rmdup.bam \
# -bed  ${workdir}/src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO.bed  > \
# ${workdir}/count/rawcount_multicov_atac.bed


# Read counts at feature count TE region (from RNA-seq)
mkdir -p ${workdir}/count
bed_RNA_TE="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_RNA-seq-v2/mm39_final/src/r/Senp6_oc_rna_RProj-013125/runtime/featurecount/table/TE_cpm_ave_FC.bed"

bedtools multicov \
-bams ${workdir}/align/Ctrl-WT1.sorted.rmdup.bam ${workdir}/align/Ctrl-WT2.sorted.rmdup.bam ${workdir}/align/Ctrl-KO1.sorted.rmdup.bam ${workdir}/align/Ctrl-KO2.sorted.rmdup.bam ${workdir}/align/R16-WT1.sorted.rmdup.bam ${workdir}/align/R16-WT2.sorted.rmdup.bam ${workdir}/align/R16-KO1.sorted.rmdup.bam ${workdir}/align/R16-KO2.sorted.rmdup.bam ${workdir}/align/R36-WT1.sorted.rmdup.bam ${workdir}/align/R36-WT2.sorted.rmdup.bam ${workdir}/align/R36-KO1.sorted.rmdup.bam ${workdir}/align/R36-KO2.sorted.rmdup.bam \
-bed  $bed_RNA_TE  > \
${workdir}/count/rawcount_multicov_atac_at_TE.bed
