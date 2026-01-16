#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J run_slurm_job_processing_atac_enrichment_atac_TE
#SBATCH -o run_slurm_job_processing_atac_enrichment_atac_TE.o
#SBATCH -e run_slurm_job_processing_atac_enrichment_atac_TE.e
#SBATCH --ntasks 1
#SBATCH --time 48:00:00
#SBATCH --mem=36G

# Load modules
module load bbc2/deeptools/deeptools-3.5.2
module load bbc2/ucsc_tools/ucsc_tools-20231127
module load bbc2/bedtools/bedtools-2.30.0

# Set variables
table="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac_RProj-012525/runtime/table"
TEbed="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/mm39_rmsk_TE.bed"
timepoints=("0h" "16h" "36h")
subtypes=("LTR" "LINE" "SINE" "DNA")

# Create bed directory BEFORE writing anything
mkdir -p ${table}/bed
cd ${table}

# Make BED files
tail -n +2 Senp6_osteoclast_dar_0h_KOvsWT_p0.05_f0_gain.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_0h_KOvsWT_gain.bed
tail -n +2 Senp6_osteoclast_dar_16h_KOvsWT_p0.05_f0_gain.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_16h_KOvsWT_gain.bed
tail -n +2 Senp6_osteoclast_dar_36h_KOvsWT_p0.05_f0_gain.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_36h_KOvsWT_gain.bed

tail -n +2 Senp6_osteoclast_dar_0h_KOvsWT_p0.05_f0_loss.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_0h_KOvsWT_loss.bed
tail -n +2 Senp6_osteoclast_dar_16h_KOvsWT_p0.05_f0_loss.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_16h_KOvsWT_loss.bed
tail -n +2 Senp6_osteoclast_dar_36h_KOvsWT_p0.05_f0_loss.txt | cut -f2-4 > bed/Senp6_osteoclast_dar_36h_KOvsWT_loss.bed

# Intersection
for time in "${timepoints[@]}"; do
    bedtools intersect \
        -a ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_gain.bed \
        -b ${TEbed} -wa -wb \
        > ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_TE_gain.bed
done

for time in "${timepoints[@]}"; do
    bedtools intersect \
        -a ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_loss.bed \
        -b ${TEbed} -wa -wb \
        > ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_TE_loss.bed
done

# Subset by TE class
for time in "${timepoints[@]}"; do
    for subtype in "${subtypes[@]}"; do
        grep -w "${subtype}" ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_TE_gain.bed \
            > ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_gain_${subtype}.bed
    done
done

for time in "${timepoints[@]}"; do
    for subtype in "${subtypes[@]}"; do
        grep -w "${subtype}" ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_TE_loss.bed \
            > ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_loss_${subtype}.bed
    done
done
