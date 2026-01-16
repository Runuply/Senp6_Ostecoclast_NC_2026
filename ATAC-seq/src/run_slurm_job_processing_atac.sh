#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J run_slurm_job_processing_atac_idr_2
#SBATCH -o run_slurm_job_processing_atac_idr_2.o
#SBATCH -e run_slurm_job_processing_atac_idr_2.e
#SBATCH --ntasks 8
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


# Setup environment
source ~/miniconda2/etc/profile.d/conda.sh
conda activate processing_rna_chip_atac_reads


# Load module, use this e.g. to search the tool, module av -w 1 bbc2 2>&1| grep 'fastqc'
# module load bbc2/samtools/samtools-1.17
# module load bbc2/hisat2/hisat2-2.2.1
# module load bbc2/trim_galore/trim_galore-0.6.10
# module load bbc2/picard/picard-3.0.0
# module load bbc2/bwa/bwa-0.7.17
# module load bbc2/macs2/macs2-2.2.7.1


# Work directory
workdir="/home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39"
cd ${workdir}


# Make new sub directory, you can make multiple at a time
#mkdir -p ../trimmedreads # align/unsort_bam align/sort_bam

##############     Processing ATAC-seq reads ###################

# Define the list of sample names
Samples=("A01" "A02" "A04" "A06" "A07" "A08" "A10" "A11")
#Samples=("A11" )

#Samples=("R36WT2")

#HISAT2_reference="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/hisat2_index/mm39"
GRCm39_GENCODE_rmsk_TE_gtf="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/GRCm39_GENCODE_rmsk_TE.gtf"
mm39_ncbiRefSeq_gtf="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/gene/mm39.ncbiRefSeq.gtf"
stringtie="/home/ye.liu/yang-secondary/ye/biotools/stringtie/stringtie"
PICARD="/varidata/research/software/BBC/picard/picard-3.0.0/picard.jar"
for Sample in "${Samples[@]}"; do

  # Define the input FASTQ files for R1 and R2
      Sample_fastq_R1="../raw/${Sample}_L000_R1_001.fastq.gz"
      Sample_fastq_R2="../raw/${Sample}_L000_R2_001.fastq.gz"

      # Define the output trimmed FASTQ files for R1 and R2
      Sample_fastq_R1_trimmed="../trimmedreads/${Sample}_L000_R1_001_val_1.fq.gz"
      Sample_fastq_R2_trimmed="../trimmedreads/${Sample}_L000_R2_001_val_2.fq.gz"
###Trim PE reads

#trim_galore --fastqc -o ../trimmedreads --paired $Sample_fastq_R1 $Sample_fastq_R2 > ${Sample}_trimgalore_summary.txt


#Align ATAC-seq PE reads with BWA
# bwa mem -t 8 /home/ye.liu/yang-secondary/ye/project/genomes/mm39/bwa_index/mm39 $Sample_fastq_R1_trimmed $Sample_fastq_R2_trimmed | samtools sort -@8 -o ${Sample}.sorted.bam
#Remove duplicate reads (picard-2.23.3)
# java -Xms8g -Xmx16g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates I=${Sample}.sorted.bam O=${Sample}.sorted.rmdup.bam M=${Sample}_dup_metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true

## Call histone peaks using MACS2
# macs2 callpeak -t ${Sample}.sorted.rmdup.bam -n $Sample -g 2.7e9 --outdir ./ -f BAMPE -s 50 --broad --broad-cutoff 0.1 #Broad peaks
# macs2 callpeak -t ${Sample}.sorted.rmdup.bam -n $Sample -g 2.7e9 --outdir ./ -f BAMPE -s 50 #Narrow peaks
done

# # Declare the ID mapping
# declare -A id_match=(
#   ["A01"]="R36-WT1"
#   ["A02"]="R36-WT2"
#   ["A04"]="R36-KO1"
#   ["A06"]="R36-KO2"
#   ["A07"]="Ctrl-WT1"
#   ["A08"]="Ctrl-WT2"
#   ["A10"]="Ctrl-KO1"
#   ["A11"]="Ctrl-KO2"
# )
#
# # Loop through the files in the folder
# for file in *; do
#   # Extract the file ID (e.g., "A01") from the file name
#   id=$(echo "$file" | cut -d '_' -f 1)
#
#   # Check if the ID exists in the mapping
#   if [[ -n ${id_match[$id]} ]]; then
#     # Replace the ID in the file name with the corresponding value from the mapping
#     new_file="${file//$id/${id_match[$id]}}"
#     echo "Renaming $file to $new_file"
#     mv "$file" "$new_file"
#   fi
# done

#Call IDR peaks
conda activate py36 ## idr-2.0.42 is installed in this environment, note, when i load this environment, and mac2 module will induce error, make sure de module the macs2

# idr --samples ${Sample}_Rep1_peaks.Peak ${Sample}_Rep2_peaks.Peak --input-file-type "PeakType" -o ${Sample}.IDR_narrowPeaks.txt --log-output-file ${Sample}_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose

# idr --samples Ctrl-KO1_peaks.narrowPeak Ctrl-KO2_peaks.narrowPeak --input-file-type narrowPeak -o Ctrl-KO-idr.bed --log-output-file  Ctrl-KO_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose
# idr --samples Ctrl-WT1_peaks.narrowPeak Ctrl-WT2_peaks.narrowPeak --input-file-type narrowPeak -o Ctrl-WT-idr.bed --log-output-file  Ctrl-WT_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose
cd idr
idr --samples ../macs2/R36-KO1_peaks.narrowPeak ../macs2/R36-KO2_peaks.narrowPeak --input-file-type narrowPeak -o R36-KO-idr.bed --log-output-file  R36-KO_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose
idr --samples ../macs2/R36-WT1_peaks.narrowPeak ../macs2/R36-WT2_peaks.narrowPeak --input-file-type narrowPeak -o R36-WT-idr.bed --log-output-file  R36-WT_narrowPeaks.IDR.log --idr-threshold 0.05 --plot --verbose



#Annotate genomic region of peaks with Homer
# annotatePeaks.pl ${Sample}_peak.bed mm39 -genomeOntology ${Sample}_genomeOntology -gtf gencode.v37.annotation.gtf > ${Sample}_annotatePeaks.bed

#Generate Average histone signals across Gene Body
# computeMatrix scale-regions --missingDataAsZero -bs 100 -p 12 -R gencode.v37.AutosomalTranscriptsOnly.forDeepTools.gtf -S $Sample_EV_Histone_Day0.bw $Sample_EV_Histone_Day5.bw $Sample_EV_Histone_Day23.bw $Sample_KO_Histone_Day0.bw $Sample_KO_Histone_Day5.bw $Sample_KO_Histone_Day23.bw -m 5000 -b 3000 -a 3000 -out $Histone_GeneBody_O786_combined.tab.gz
# plotHeatmap  --heatmapWidth 6 -m $Histone_GeneBody_O786_combined.tab.gz --perGroup -out $Histone_GeneBody_O786_combined.pdf --missingDataColor "grey" --colorMap "Color" --heatmapHeight 10
module load bbc2/bedtools/bedtools-2.30.0

cd /home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/src/r/Senp6_oc_atac/DiffBind/table/
#Make bed files from the Diifbind output
awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_16h_WTvsKO_p0.05.txt > bed/Senp6_osteoclast_dar_16h_WTvsKO_p0.05.bed

# Make KO upregulated bed files from Diffbind

awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_0h_WTvsKO_p0.05_up.txt > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated.bed
awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_16h_WTvsKO_p0.05_up.txt > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated.bed
awk 'NR>1 {print $3 "\t" $4 "\t" $5}' Senp6_osteoclast_dar_36h_WTvsKO_p0.05_up.txt > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated.bed


# Make bed files
#bedtools intersect -a atac_peaks.bed -b mm39_rmsk_TE.bed -wa > atac_peaks_in_TE.bed
TEbed="/home/ye.liu/yang-secondary/ye/project/genomes/mm39/rmsk/mm39_rmsk_TE.bed"
bedtools intersect -a bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed
bedtools intersect -a bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed
bedtools intersect -a bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated.bed -b $TEbed -wa -wb > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed
# Make bed file, subset by TE subtype

grep "LTR" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LTR.bed
grep "LINE" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LINE.bed
grep "SINE" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_SINE.bed
grep "DNA" bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_DNA.bed

grep "LTR" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LTR.bed
grep "LINE" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LINE.bed
grep "SINE" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_SINE.bed
grep "DNA" bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_DNA.bed

grep "LTR" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LTR.bed
grep "LINE" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LINE.bed
grep "SINE" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_SINE.bed
grep "DNA" bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_TE.bed > bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_DNA.bed


#!/bin/bash

# deeptools
module load bbc2/deeptools/deeptools-3.5.2
cd /home/ye.liu/yang-secondary/ye/project/Senp6_osteoclast_project/L2S6_ATAC-seq-v2/L2S6_ATAC-seq_prep3/mm39/deeptools

# Define the input BED file
BED_FILE="../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_p0.05.bed"

# Define the BigWig files as an array
BIGWIG_FILES=(
    "../bigwig_RPGC/Ctrl-KO1_RPGC.bw"
    "../bigwig_RPGC/Ctrl-KO2_RPGC.bw"
    "../bigwig_RPGC/Ctrl-WT1_RPGC.bw"
    "../bigwig_RPGC/Ctrl-WT2_RPGC.bw"
)

# Loop through each BigWig file
for BIGWIG in "${BIGWIG_FILES[@]}"; do
    # Extract the base name of the BigWig file (e.g., Ctrl-KO1_RPGC)
    BASE_NAME=$(basename "$BIGWIG" .bw)

    # Construct the output file name in the current directory
    OUTPUT_FILE="ATAC_diffPeaks_matrix_${BASE_NAME}.gz"

    # Run computeMatrix for the current BigWig file
    computeMatrix reference-point \
        --referencePoint center \
        --missingDataAsZero \
        -bs 100 \
        -p 12 \
        -R $BED_FILE \
        -S $BIGWIG \
        -b 3000 \
        -a 3000 \
        -out $OUTPUT_FILE

    echo "Processed: $BIGWIG -> $OUTPUT_FILE"
done

echo "All computeMatrix tasks completed!"



# Define the list of matrix files in the current directory
MATRIX_FILES=(
    "ATAC_diffPeaks_matrix_Ctrl-KO1_RPGC.gz"
    "ATAC_diffPeaks_matrix_Ctrl-KO2_RPGC.gz"
    "ATAC_diffPeaks_matrix_Ctrl-WT1_RPGC.gz"
    "ATAC_diffPeaks_matrix_Ctrl-WT2_RPGC.gz"
)

# Loop through each matrix file
for MATRIX in "${MATRIX_FILES[@]}"; do
    # Extract the base name of the matrix file (e.g., Ctrl-KO1_RPGC)
    BASE_NAME=$(basename "$MATRIX" .gz)

    # Generate the heatmap using plotHeatmap
    plotHeatmap \
        --heatmapWidth 3 \
        -m $MATRIX \
        --perGroup \
        -out "${BASE_NAME}_heatmap_only.pdf" \
        --missingDataColor "grey" \
        --colorMap "BuGn" \
        --heatmapHeight 10 \
         --whatToShow 'heatmap and colorbar' # no profile

    echo "Generated heatmap for $BASE_NAME"
done


####
#combine all the bw
computeMatrix reference-point \
    --referencePoint center \
    --missingDataAsZero \
    -bs 100 \
    -p 12 \
    -R ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LTR.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_SINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_DNA.bed \
    -S ../bigwig_RPGC/Ctrl-KO1_RPGC.bw ../bigwig_RPGC/Ctrl-KO2_RPGC.bw ../bigwig_RPGC/Ctrl-WT1_RPGC.bw ../bigwig_RPGC/Ctrl-WT2_RPGC.bw \
    -b 3000 \
    -a 3000 \
    --averageTypeBins max \
    -out ATAC_diffPeaks_matrix_0h_t.gz

plotHeatmap \
                -m ATAC_diffPeaks_matrix_0h_t.gz \
                --heatmapWidth 3 \
                --heatmapHeight 10 \
                --colorMap "BuGn" \
                --missingDataColor "grey" \
                --regionsLabel "LTR" "LINE" "SINE" "DNA" \
                --outFileName ATAC_diffPeaks_heatmap_0h_t.pdf

computeMatrix reference-point \
                --referencePoint center \
                --missingDataAsZero \
                -bs 100 \
                -p 12 \
                -R ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LTR.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_LINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_SINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_16h_WTvsKO_upregulated_DNA.bed \
                -S ../bigwig_RPGC/R16-KO1_RPGC.bw ../bigwig_RPGC/R16-KO2_RPGC.bw ../bigwig_RPGC/R16-WT1_RPGC.bw ../bigwig_RPGC/R16-WT2_RPGC.bw \
                -b 3000 \
                -a 3000 \
                --averageTypeBins max \
                -out ATAC_diffPeaks_matrix_16h_t.gz

plotHeatmap \
                -m ATAC_diffPeaks_matrix_16h_t.gz \
                --heatmapWidth 3 \
                --heatmapHeight 10 \
                --colorMap "BuGn" \
                --missingDataColor "grey" \
                --regionsLabel "LTR" "LINE" "SINE" "DNA" \
                --outFileName ATAC_diffPeaks_heatmap_16h_t.pdf

computeMatrix reference-point \
                                --referencePoint center \
                                --missingDataAsZero \
                                -bs 100 \
                                -p 12 \
                                -R ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LTR.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_LINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_SINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_36h_WTvsKO_upregulated_DNA.bed \
                                -S ../bigwig_RPGC/R36-KO1_RPGC.bw ../bigwig_RPGC/R36-KO2_RPGC.bw ../bigwig_RPGC/R36-WT1_RPGC.bw ../bigwig_RPGC/R36-WT2_RPGC.bw \
                                -b 3000 \
                                -a 3000 \
                                --averageTypeBins max \
                                -out ATAC_diffPeaks_matrix_36h_t.gz

plotHeatmap \
                                -m ATAC_diffPeaks_matrix_36h_t.gz \
                                --heatmapWidth 3 \
                                --heatmapHeight 10 \
                                --colorMap "BuGn" \
                                --missingDataColor "grey" \
                                --regionsLabel "LTR" "LINE" "SINE" "DNA" \
                                --outFileName ATAC_diffPeaks_heatmap_36h_t.pdf


    computeMatrix reference-point \
        --referencePoint center \
        --missingDataAsZero \
        -bs 100 \
        -p 12 \
        -R ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LTR.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_LINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_SINE.bed ../src/r/Senp6_oc_atac/DiffBind/table/bed/Senp6_osteoclast_dar_0h_WTvsKO_upregulated_DNA.bed \
        -S ../bigwig_RPGC/Ctrl-WT_Combined.bw ../bigwig_RPGC/Ctrl-WT_Combined.bw  \
        -b 3000 \
        -a 3000 \
        --averageTypeBins max \
        -out ATAC_diffPeaks_matrix_0h.gz


  plotHeatmap \
      -m ATAC_diffPeaks_matrix_0h_t.gz \
      --heatmapWidth 3 \
      --heatmapHeight 10 \
      --perGroup \
      --colorMap "BuGn" \
      --missingDataColor "grey" \
      --whatToShow "heatmap and colorbar" \
      --outFileName ATAC_diffPeaks_heatmap_0h.pdf


# --perGroup will make everything in one column
# --whatToShow "heatmap and colorbar" this will not give you frofile
# make sure there is nothing after `\`


  plotHeatmap \
      -m ATAC_diffPeaks_matrix_0h.gz \
      --heatmapWidth 2 \
      --heatmapHeight 20 \
      --colorMap "BuGn" \
      --missingDataColor "grey" \
      --whatToShow "heatmap and colorbar" \
      --zMin 0 \
      --zMax 30 \
      --regionsLabel "LTR" "LINE" "SINE" "DNA" \
      --samplesLabel "WT" "KO" \
      --outFileName ATAC_diffPeaks_heatmap_0h.pdf

  plotHeatmap \
          -m ATAC_diffPeaks_matrix_0h.gz \
          --heatmapWidth 2 \
          --heatmapHeight 20 \
          --colorMap "BuGn" \
          --missingDataColor "grey" \
          --whatToShow "heatmap and colorbar" \
          --regionsLabel "LTR" "LINE" "SINE" "DNA" \
          --samplesLabel "WT" "KO" \
          --outFileName ATAC_diffPeaks_heatmap_0h.pdf

######################################################################
end_time=$(date +"%A, %F %T %Z")

# Calculate runtime in seconds
start_timestamp=$(date -d "$start_time" +"%s")
end_timestamp=$(date -d "$end_time" +"%s")
runtime=$((end_timestamp - start_timestamp))

echo "############## Job Summary ##############"
echo "Author: $current_user @ $hostname"
echo "Job ID: $job_id"
echo "Job Name: $job_name"
echo "Start time: $start_time"
echo "End time: $end_time"
echo "Runtime: $(date -u -d @"$runtime" +"%H:%M:%S")"
echo "#########################################"

# Adding a funny message
funny_messages=("Job well done! You rocked it!"
                "Time to celebrate with some ice cream!"
                "Congratulations, you're awesome!"
                "Great job! Treat yourself to a well-deserved break!"
                "You nailed it! Take a victory lap!"
                "Hooray! Dance like no one's watching!"
                "Bravo! You're a superstar!"
                "Mission accomplished! You deserve a high-five!"
                "You did it! Now go and conquer the world!"
                "You're officially a job-completing expert!"
                "Way to go! Your brilliance shines through!"
                "Woot woot! Party time!")
random_message=${funny_messages[RANDOM % ${#funny_messages[@]}]}
echo "------------------------------"
echo "   \   ^_^  */"
echo "    \ (Z O) /"
echo "     \  V  /"
echo "      \   /"
echo "       | |"
echo "Funny Message: $random_message"
######################################################################
