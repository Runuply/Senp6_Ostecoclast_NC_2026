#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J run_slurm_job_processing_atac_bw_deeptools
#SBATCH -o run_slurm_job_processing_atac_bw_deeptools.o
#SBATCH -e run_slurm_job_processing_atac_bw_deeptools.e
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

#Software versions
# samtools-1.17
# hisat2-2.2.1
# stringtie-2.2.1
# picard-3.0.0
# R-4.0.2
# python-3.8.1
# bwa-0.7.17
# HOMER-4.11

module load bbc2/deeptools/deeptools-3.5.2
module load bbc2/ucsc_tools/ucsc_tools-20231127
module load bbc2/bedtools/bedtools-2.30.0
# Work directory
workdir="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39"
cd ${workdir}

#### MAKE BIGWIG FROM BAM AND PLOT GENEBODY HISTONE LEVELS ####
cd ${workdir}/align
# for i in *.sorted.rmdup.bam; do bamCoverage --bam $i -o ${workdir}/bigwig_RPGC/${i/.sorted.rmdup.bam/_RPGC.bw} -p 16 --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2723414844 --ignoreForNormalization chrX chrM chrY --extendReads;done

# # If you want to check the chromsome size, please
# cd /home/ye.liu/yang-secondary/ye/project/genomes/mm39
# awk '$1 ~ /^chr[0-9XY]+$/ {sum += $2} END {print sum}' chromosome_sizes.txt
# 2723414844 ## tested 01/24/25

#
#combine biological replicates for visualization
#cd ${workdir}/bigwig_RPGC

# # Step 1: Merge replicates for each group
#for condition in Ctrl R16 R36; do
 #  for genotype in WT KO; do
 #    bigWigMerge ${condition}-${genotype}1_RPGC.bw ${condition}-${genotype}2_RPGC.bw ${condition}-${genotype}_Combined.bedGraph
 #  done
 #done

 # Step 2: Sort the merged bedGraph files
 #for i in *_Combined.bedGraph; do
 #  sort -k1,1 -k2,2n $i > ${i/.bedGraph/.sorted.bedGraph}
 #done

# Step 3: Convert sorted bedGraph files to bigWig
 #for i in *_Combined.sorted.bedGraph; do
#   bedGraphToBigWig $i /home/ye.liu/yang-secondary/ye/project/genomes/mm39/chromosome_sizes.txt ${i/.sorted.bedGraph/.bw}
# done


#deeptools
cd /home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac_RProj-012525/runtime/table
tail -n +2 Senp6_osteoclast_dar_0h_KOvsWT.txt | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_0h_KOvsWT.bed
tail -n +2 Senp6_osteoclast_dar_16h_KOvsWT.txt | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_16h_KOvsWT.bed
tail -n +2 Senp6_osteoclast_dar_36h_KOvsWT.txt  | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_36h_KOvsWT.bed

tail -n +2 Senp6_osteoclast_dar_0h_KOvsWT_p0.05_f0.txt | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_0h_KOvsWT_upregulated.bed
tail -n +2 Senp6_osteoclast_dar_16h_KOvsWT_p0.05_f0.txt | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_16h_KOvsWT_upregulated.bed
tail -n +2 Senp6_osteoclast_dar_36h_KOvsWT_p0.05_f0.txt  | cut -f 2,3,4 > bed/Senp6_osteoclast_dar_36h_KOvsWT_upregulated.bed


# Set variables
TEbed="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/mm39_rmsk_TE.bed"
table="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac_RProj-012525/runtime/table"
timepoints=("0h" "16h" "36h")
subtypes=("LTR" "LINE" "SINE" "DNA")

# Step 1: Intersect DAR BED files with TE BED file
mkdir -p ${table}/bed
for time in "${timepoints[@]}"; do
    input_bed="${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated.bed"
    output_te_bed="${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_TE.bed"
    bedtools intersect -a ${input_bed} -b ${TEbed} -wa -wb > ${output_te_bed}
done


# Step 2: Subset by TE subtype
for time in "${timepoints[@]}"; do
    for subtype in "${subtypes[@]}"; do
        grep "${subtype}" ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_TE.bed > \
            ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_${subtype}.bed
    done
done

# Step 3: Generate matrix and heatmap with combined bw
cd ${workdir}/deeptools
for time in "${timepoints[@]}"; do
    # Set bigwig files based on timepoint
    if [ "$time" == "0h" ]; then
        bw_files="../bigwig_RPGC/Ctrl-WT_Combined.bw ../bigwig_RPGC/Ctrl-KO_Combined.bw"
    elif [ "$time" == "16h" ]; then
        bw_files="../bigwig_RPGC/R16-WT_Combined.bw ../bigwig_RPGC/R16-KO_Combined.bw"
    elif [ "$time" == "36h" ]; then
        bw_files="../bigwig_RPGC/R36-WT_Combined.bw ../bigwig_RPGC/R36-KO_Combined.bw"
    fi

    computeMatrix reference-point \
        --referencePoint center \
        --missingDataAsZero \
        -bs 100 \
        -p 12 \
        -R ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_LTR.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_LINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_SINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_DNA.bed \
        -S $bw_files \
        -b 3000 \
        -a 3000 \
        --averageTypeBins max \
        -out ATAC_diffPeaks_matrix_${time}.gz

    plotHeatmap \
        -m ATAC_diffPeaks_matrix_${time}.gz \
        --heatmapWidth 3 \
        --heatmapHeight 10 \
        --colorMap "BuGn" \
        --missingDataColor "grey" \
        --regionsLabel "LTR" "LINE" "SINE" "DNA" \
        --samplesLabel "WT" "KO" \
        --sortRegions ascend \
        --outFileName ATAC_diffPeaks_heatmap_${time}.pdf
done

# Step 4: Generate matrix and heatmap with reps

cd ${workdir}/deeptools
for time in "${timepoints[@]}"; do
    # Set bigwig files based on timepoint
    if [ "$time" == "0h" ]; then
        bw_files="../bigwig_RPGC/Ctrl-WT1_RPGC.bw ../bigwig_RPGC/Ctrl-WT2_RPGC.bw ../bigwig_RPGC/Ctrl-KO1_RPGC.bw ../bigwig_RPGC/Ctrl-KO2_RPGC.bw"
    elif [ "$time" == "16h" ]; then
        bw_files="../bigwig_RPGC/R16-WT1_RPGC.bw ../bigwig_RPGC/R16-WT2_RPGC.bw ../bigwig_RPGC/R16-KO1_RPGC.bw ../bigwig_RPGC/R16-KO2_RPGC.bw"
    elif [ "$time" == "36h" ]; then
        bw_files="../bigwig_RPGC/R36-WT1_RPGC.bw ../bigwig_RPGC/R36-WT2_RPGC.bw ../bigwig_RPGC/R36-KO1_RPGC.bw ../bigwig_RPGC/R36-KO2_RPGC.bw"
    fi

    computeMatrix reference-point \
        --referencePoint center \
        --missingDataAsZero \
        -bs 100 \
        -p 12 \
        -R ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_LTR.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_LINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_SINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_DNA.bed \
        -S $bw_files \
        -b 3000 \
        -a 3000 \
        --averageTypeBins max \
        -out ATAC_diffPeaks_matrix_${time}.gz

    plotHeatmap \
        -m ATAC_diffPeaks_matrix_${time}.gz \
        --heatmapWidth 3 \
        --heatmapHeight 10 \
        --colorMap "BuGn" \
        --missingDataColor "grey" \
        --regionsLabel "LTR" "LINE" "SINE" "DNA" \
        --samplesLabel "WT1" "WT2" "KO1" "KO2" \
        --sortRegions ascend \
        --outFileName ATAC_diffPeaks_heatmap_${time}_rep.pdf
done


# deeptools of ATAC-seq peaks ovelapped with gained intergenic or intronic TEs from RNA-seq (FC>2)

# Set variables
TEbed="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/mm39_rmsk_TE.bed"
table="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac_RProj-012525/runtime/table"
timepoints=("0h" "16h" "36h")
subtypes=("LTR" "LINE" "SINE" "DNA")

# Step 1: Intersect DAR-TE BED with the upregulated TEs (RNA)-here are the total TE transcriptome of file TE_cpm_ave_FC.bed from Featurecount
mkdir -p ${table}/bed
for time in "${timepoints[@]}"; do
    input_bed_a="${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_TE.bed"
    input_bed_b="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_RNA-seq-v2/mm39_final/src/r/Senp6_oc_rna_RProj-013125/runtime/featurecount/table/TE_cpm_ave_FC_2.bed"
    output_te_bed="${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_TE_RNAFC2.bed"
    bedtools intersect -a ${input_bed_a} -b ${input_bed_b} -wa -wb > ${output_te_bed}
done

# Step 2: Subset by TE subtype
for time in "${timepoints[@]}"; do
    for subtype in "${subtypes[@]}"; do
        grep "${subtype}" ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_TE_RNAFC2.bed > \
            ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_RNAFC2_${subtype}.bed
    done
done


# Step 3: Generate matrix and heatmap with combined bw
cd ${workdir}/deeptools
for time in "${timepoints[@]}"; do
    # Set bigwig files based on timepoint
    if [ "$time" == "0h" ]; then
        bw_files="../bigwig_RPGC/Ctrl-WT_Combined.bw ../bigwig_RPGC/Ctrl-KO_Combined.bw"
    elif [ "$time" == "16h" ]; then
        bw_files="../bigwig_RPGC/R16-WT_Combined.bw ../bigwig_RPGC/R16-KO_Combined.bw"
    elif [ "$time" == "36h" ]; then
        bw_files="../bigwig_RPGC/R36-WT_Combined.bw ../bigwig_RPGC/R36-KO_Combined.bw"
    fi

    computeMatrix reference-point \
        --referencePoint center \
        --missingDataAsZero \
        -bs 100 \
        -p 12 \
        -R ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_RNAFC2_LTR.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_RNAFC2_LINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_RNAFC2_SINE.bed \
           ${table}/bed/Senp6_osteoclast_dar_${time}_KOvsWT_upregulated_RNAFC2_DNA.bed \
        -S $bw_files \
        -b 3000 \
        -a 3000 \
        --averageTypeBins max \
        -out ATAC_diffPeaks_RNAFC2_matrix_${time}.gz

    plotHeatmap \
        -m ATAC_diffPeaks_RNAFC2_matrix_${time}.gz \
        --heatmapWidth 3 \
        --heatmapHeight 10 \
        --colorMap "BuGn" \
        --missingDataColor "grey" \
        --regionsLabel "LTR" "LINE" "SINE" "DNA" \
        --samplesLabel "WT" "KO" \
        --sortRegions ascend \
        --outFileName ATAC_diffPeaks_RNAFC2_heatmap_${time}.pdf
done

# cd /home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac/DiffBind/table/
# #Make bed files from the Diifbind output
# awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_16h_WTvsKO_p0.05.txt > bed/Senp6_osteoclast_dar_16h_WTvsKO_p0.05.bed
# awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_0h_WTvsKO.txt > bed/Senp6_osteoclast_dar_0h_WTvsKO.bed
# awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_16h_WTvsKO.txt > bed/Senp6_osteoclast_dar_16h_WTvsKO.bed
# awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_36h_WTvsKO.txt > bed/Senp6_osteoclast_dar_36h_WTvsKO.bed
#
# # Make KO upregulated bed files from Diffbind
# awk 'NR > 1 && $10 > 1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_0h_WTvsKO_p0.05.txt > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated.bed
# awk 'NR > 1 && $10 > 1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_16h_WTvsKO_p0.05.txt > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated.bed
# awk 'NR > 1 && $10 > 1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_36h_WTvsKO_p0.05.txt > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated.bed
#
# # Make bed files
# #bedtools intersect -a atac_peaks.bed -b mm39_rmsk_TE.bed -wa > atac_peaks_in_TE.bed
# TEbed="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/mm39_rmsk_TE.bed"
# bedtools intersect -a bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed
# bedtools intersect -a bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed
# bedtools intersect -a bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed
# # Make bed file, subset by TE subtype
# grep "LTR" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LTR.bed
# grep "LINE" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LINE.bed
# grep "SINE" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_SINE.bed
# grep "DNA" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_DNA.bed
#
# grep "LTR" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LTR.bed
# grep "LINE" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LINE.bed
# grep "SINE" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_SINE.bed
# grep "DNA" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_DNA.bed
#
# grep "LTR" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LTR.bed
# grep "LINE" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LINE.bed
# grep "SINE" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_SINE.bed
# grep "DNA" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_DNA.bed
#
# cd ${workdir}/deeptools

# computeMatrix reference-point \
#     --referencePoint center \
#     --missingDataAsZero \
#     -bs 100 \
#     -p 12 \
#     -R ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LTR.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_SINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_DNA.bed \
#     -S ../bigwig_RPGC/Ctrl-WT_Combined.bw ../bigwig_RPGC/Ctrl-WT_Combined.bw  \
#     -b 3000 \
#     -a 3000 \
#     --averageTypeBins max \
#     -out ATAC_diffPeaks_matrix_0h.gz
#
# plotHeatmap \
#       -m ATAC_diffPeaks_matrix_0h.gz \
#       --heatmapWidth 2 \
#       --heatmapHeight 20 \
#       --colorMap "BuGn" \
#       --missingDataColor "grey" \
#       --whatToShow "heatmap and colorbar" \
#       --regionsLabel "LTR" "LINE" "SINE" "DNA" \
#       --samplesLabel "WT" "KO" \
#       --outFileName ATAC_diffPeaks_heatmap_0h.pdf
