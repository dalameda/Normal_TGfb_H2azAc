#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=ChIP_pipeline            # Job name
#SBATCH --output=ChIP_pipeline.%j.out       # File to write standard output
#SBATCH --error=ChIP_pipeline.%j.err        # File to write standard error
#SBATCH --mem=40G                               # Memory allocated to the job
#SBATCH --cpus-per-task=8                       # Number of CPU cores per task

# # LAUNCH like this:
# conda env create -f scripts/chip.yaml --name ChIP
# script="draft_2.sbs"
# sample_info=sample.info.20231013.txt
# LIST=($(cut -f1 data/${sample_info} | cat )); NUMFILES=${#LIST[@]}; ZBNUMFILES=$(($NUMFILES - 1))
# if [ $ZBNUMFILES -ge 0 ]; then sbatch --mem=40G --array=0-$ZBNUMFILES scripts/$script data/$sample_info; fi

# Load necessary software modules
module load bowtie/2.5.1                        # For aligning sequences to the genome
module load samtools/1.13                       # For manipulating SAM/BAM files
module load bedtools/2.30.0                     # For genome arithmetic
module load deeptools/3.5.1                     # For generating BigWig files
module load cutadapt/4.4                        # For trimming adapter sequences from reads
module load subread/2.0.6                       # For read counting (featureCounts)
module load miniconda/3                         # To manage the software environment


# Activate the Conda environment
source activate ChIP # chip.yaml for environment configuration

# # chip.yaml
# channels:
#   - conda-forge
#   - bioconda
#   - r
# dependencies:
#   - bioconda::macs2
#   - gawk
#   - conda-forge::r-r.utils==2.12.0
#   - bioconda::ucsc-bedgraphtobigwig
#   - bioconda::cutadapt
#   - bioconda::bioconductor-affy
#   - bioconda::homer
#   - r-mass
#   - r::r-essentials
#   - bioconductor-limma
#   - bioconda::bioconductor-deseq2
#   - bioconda::bioconductor-edger

# Define file paths and variables
sample_info="data/sample.info.20231013.txt"     # Metadata file with sample information
genome="mm10"                                   # Reference genome
genome_mix=/mnt2/fscratch/users/oncoh_011_ibima/dalameda/snakemake/chip/data/refGenomes/mm10.fa
genome_index=/mnt2/fscratch/users/oncoh_011_ibima/dalameda/snakemake/chip/data/refGenomes/mm10_index
# Specify the blacklist file location
blacklist="data/mm10.reordered.blacklist.bed"

# Prepare output directories
trim_dir="data/fastq"
map_dir="data/bam"
bam_dir="data/bam"
bigwig_dir="data/bigwig_files"
peak_dir="data/peak_files"
mkdir -p $trim_dir $map_dir $bam_dir $bigwig_dir $peak_dir

# Read sample information and assign variables for each sample
sample_list=($(cut -f1 $sample_info))
genome_list=($(cut -f2 $sample_info))
spike_genome_list=($(cut -f3 $sample_info))
Rd1list=($(cut -f4 $sample_info | cat ))
Rd2list=($(cut -f5 $sample_info | cat ))
conditions=$(cut -f6 $sample_info | sort | uniq)

# Extract details for the current sample based on SLURM task ID
sample=${sample_list[$SLURM_ARRAY_TASK_ID]}
genome=${genome_list[$SLURM_ARRAY_TASK_ID]}
spike_genome=${spike_genome_list[$SLURM_ARRAY_TASK_ID]}
fastqR1=${Rd1list[$SLURM_ARRAY_TASK_ID]}
fastqR2=${Rd2list[$SLURM_ARRAY_TASK_ID]}

# Step 1: Adapter trimming with Cutadapt
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  -o $trim_dir/${sample}_1.trimmed.fastq \
  -p $trim_dir/${sample}_2.trimmed.fastq $fastqR1 $fastqR2

# Align reads to the reference genome and spike-in sequences
# Note: Ensure genome and spike-in indices are prepared in advance

# Step 2: Mapping with Bowtie2
bowtie2 -p $SLURM_CPUS_PER_TASK --time -k 1 -x ${genome_index} \
  -1 $trim_dir/${sample}_1.trimmed.fastq \
  -2 $trim_dir/${sample}_2.trimmed.fastq \
  -S $map_dir/${sample}_sample.sam

# Convert SAM to BAM, sort, remove duplicates, and index
samtools view -F 0x4 -h $map_dir/${sample}_sample.sam | samtools sort -o $bam_dir/${sample}_sample_sorted.bam
samtools rmdup $bam_dir/${sample}_sample_sorted.bam $bam_dir/${sample}_sample_sorted_rmdup.bam
samtools index $bam_dir/${sample}_sample_sorted_rmdup.bam

# Step 3: Removing blacklisted regions
# remove unwanted chroms
awk '{print $1}' ${genome_fai} | grep -vE 'chrUn|random|chrM|chrEBV' - > chrom_list.txt
chroms=$(tr '\n' ' ' < chrom_list.txt)
samtools view -b -h $bam_dir/${sample}_sample_sorted_rmdup.bam $chroms -o $bam_dir/${sample}_sample_sorted_rmdup_chr.bam
rm chrom_list.txt
samtools index $bam_dir/${sample}_sample_sorted_rmdup_chr.bam
# remove blacklisted
bedtools intersect -v -abam $bam_dir/${sample}_sample_sorted_rmdup_chr.bam -b $blacklist > $bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam
samtools index $bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam

# Generate BigWig files for visualization with deepTools
bamCoverage --binSize 5 --normalizeUsing CPM --exactScaling \
  -b $bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam \
  -of bigwig -o $bigwig_dir/${sample}_sample.bw \
  --numberOfProcessors $SLURM_CPUS_PER_TASK

# Step 4: Peak Calling with MACS2 for the Spike-in Data
# Perform peak calling separately for each dataset (sample and spike-in) using MACS2.
# This step identifies enriched regions (peaks) in your ChIP-seq data.
macs2 callpeak -t $bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam \
-f BAM -g dm -n ${sample}_sample --outdir $peak_dir \
--broad --broad-cutoff 0.1

# Add the processed bam file to a list for later use in multi-sample analysis
echo "$bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam" >> $bam_dir/bam.list.txt

# Step 5: Annotate Peaks using HOMER
# Annotate the identified peaks with known genomic features. This step helps to understand
# the potential functional implications of the binding sites.
perl /mnt2/fscratch/users/oncoh_011_ibima/dalameda/conda_env/ChIP/share/homer/.//configureHomer.pl -install mm10
mkdir -p ${peak_dir}/HOMER/peakAnnotation/
annotatePeaks.pl \
        ${peak_dir}/${sample}_sample_peaks.broadPeak \
        ${genome} \
        -gid \
        -cpu ${SLURM_CPUS_PER_TASK} \
        -annStats ${peak_dir}/HOMER/peakAnnotation/${sample}_sample_broadPeak.annotateStats.txt \
        > ${peak_dir}/HOMER/peakAnnotation/${sample}_sample_broadPeak.annotatePeaks.txt

# Create tag directories for each sample using HOMER. This step is essential for
# quantitative analysis of ChIP-seq data, like creating tag density plots or finding motifs.
mkdir -p ${peak_dir}/HOMER/tagDirectory/
makeTagDirectory ${peak_dir}/HOMER/tagDirectory/${sample}_sample \
    $bam_dir/${sample}_sample_sorted_rmdup_chr_black.bam \
    -genome ${genome} \
    -checkGC

# Step 6: Differential Binding Analysis Preparation
# Prepare for differential binding analysis by identifying consensus peaks across replicates for each condition.
# This step is vital for pinpointing genomic regions consistently enriched across replicates,
# laying the groundwork for subsequent analyses that compare binding patterns between conditions.

# Create necessary directories to store output files from the consensus peak analysis
mkdir -p ${peak_dir}/Consensus/
mkdir -p ${peak_dir}/Consensus/sample_peaks.broadPeak/

# Initialize empty strings to accumulate lists of consensus peak files and broad peak files for all samples and conditions
cobound_list=""  # Will store paths to consensus peak files for each condition
broadpeak_list=""  # Will store paths to all individual broad peak files from MACS2

# Iterate over each experimental condition to process and identify consensus peaks
for condition in "${conditions[@]}"; do
    # Retrieve the list of sample IDs associated with the current condition from the sample information file
    replicas=($(awk -v cond="$condition" '$6 == cond {print $1}' ${sample_info}))
    replica_count=$(awk -v condition="$condition" '$6 == condition {print $0}' ${sample_info} | wc -l)
    replica_count=$((replica_count - 1))  # Adjust count for zero-based indexing, required for specifying the -cobound parameter in mergePeaks

    # Collect paths to peak files for all replicates under the current condition to be merged
    files_sp=""  # Accumulates peak file paths for merging
    for replica in "${replicas[@]}"; do
        files_sp+=" ${peak_dir}/${replica}_sample_peaks.broadPeak"
    done

    # Execute mergePeaks to identify peaks common across replicates of the same condition, generating a matrix and Venn diagram of shared peaks
    mergePeaks -d given $files_sp -matrix ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition} \
               -venn ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition}.venn.txt > ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition}

    # Refine the list of consensus peaks by ensuring they are cobound by all replicates
    # This step filters the original MACS2 output to retain only those peaks present in the consensus set for the condition
    mergePeaks -d given $files_sp -prefix ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition} -cobound ${replica_count}

    # Further process each file to extract consensus peaks, contributing to a comprehensive list of peaks per sample and condition
    for file in $files_sp; do
        awk -F'\t' 'NR==FNR{a[$1]; next} ($4 in a){print}' \
        ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition}.coBoundBy${replica_count}.txt \
        $file > ${peak_dir}/Consensus/sample_peaks.broadPeak/$(basename ${file%_sample_peaks.broadPeak}).coBoundBy${replica_count}.txt
        broadpeak_list+=" ${peak_dir}/Consensus/sample_peaks.broadPeak/$(basename ${file%_sample_peaks.broadPeak}).coBoundBy${replica_count}.txt"
    done

    # Compile a list of all consensus peaks across conditions for downstream analysis
    cobound_list+=" ${peak_dir}/Consensus/sample_peaks.broadPeak/${condition}.coBoundBy${replica_count}.txt"
done

# The output from this script includes:
# - Consensus peak files for each condition, suitable for differential binding analysis.
# - A comprehensive list of broadPeak files from MACS2, ready for further comparative studies.

# Step 7: Merge Replicate Peaks to Identify Consensus Peaks and Quantify Binding Intensities

# Aggregate peak information across all conditions and replicates around a comprehensive list of consensus peaks.
# Prepare these consensus peaks for quantitative analysis by formatting them into SAF and using featureCounts to quantify binding intensities. This process enables a comparative analysis of binding across samples or conditions based on peak intensity.

# Create a comprehensive list of consensus peaks by merging individual broadPeak files.
mergePeaks -d given $broadpeak_list > ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks.txt

# Prepare for differential binding analysis by formatting the consensus peaks into a SAF file.
# The SAF format is required by featureCounts for read quantification and includes columns for gene ID,
# chromosome, start and end positions, and strand information.
echo -e "GeneID\tChr\tStart\tEnd\tStrand" > ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks.saf

# Convert the consensus peaks into SAF format, extracting necessary fields from the merged peak file.
# This conversion is a preparatory step that structures the data for efficient read counting.
awk 'NR>1 {print $1"\t"$2"\t"$3"\t"$4"\t"$5}' ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks.txt >> ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks.saf

# Accumulate a list of BAM files for read counting. This list includes all BAM files generated from
# aligned reads across the entire dataset, encompassing all conditions and replicates.
while IFS= read -r bam; do
  files_bam+=" ${bam}"
done < $bam_dir/bam.list.txt

# Execute featureCounts to quantify reads aligning within consensus peaks for each sample.
# This command applies a fractional overlap of 0.2 to consider a read as contributing to peak intensity,
# ensuring that only reads significantly overlapping with peaks are counted. The output is a matrix of
# read counts per consensus peak across all samples, facilitating downstream differential binding analysis.
featureCounts -p -O --fracOverlap 0.2 -T ${SLURM_CPUS_PER_TASK} -a ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks.saf -F SAF -o ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks_counts.txt ${files_bam}

# This quantification step is pivotal for comparing binding intensities across samples or conditions, laying
# the foundation for identifying differentially bound regions reflective of biological differences driven by
# experimental conditions.

# Step 8: Annotate Consensus Peaks for Downstream Analysis
# Use HOMER to annotate the consensus peaks, providing biological context for the identified regions
annotatePeaks.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks_counts.txt mm10 > ${peak_dir}/Consensus/sample_peaks.broadPeak/rep.over.cond_total.peaks_counts.annot.txt

# Merge tagDirectories replicates

makeTagDirectory ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_sample/ -d ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_ChIP18_sample/ ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_ChIP19_sample/
makeTagDirectory ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_sample/ -d ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_ChIP18_sample/ ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_ChIP19_sample/

# # Files of Diff Bound peaks for Motif analysis
#
# ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.txt # Down binding TGFb_H2azAc vs Normal_H2azAc
# ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.txt # Up binding TGFb_H2azAc vs Normal_H2azAc
#
# # Background for Motif analysis
#
# ChIP-Seq.pf.dds_p0.05_lgFC0.75NOCHANGE_TGFbvsNormal.woDOWN.peaks.txt # All peaks TGFb_H2azAc vs Normal_H2azAc without DOWN binding Peaks
# ChIP-Seq.pf.dds_p0.05_lgFC0.75NOCHANGE_TGFbvsNormal.woUP.peaks.txt # All peaks binding TGFb_H2azAc vs Normal_H2azAc without UP binding Peaks

# DEG Peaks Motifs (from DEG results)

findMotifsGenome.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.txt \
$genome ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks_MotifAnalysis.custom.bg -size given -mask -bg ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75NOCHANGE_TGFbvsNormal.woDOWN.peaks.txt
findMotifsGenome.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.txt \
$genome ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks_MotifAnalysis.custom.bg -size given -mask -bg ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75NOCHANGE_TGFbvsNormal.woUP.peaks.txt

# DEG Peaks Histograms

annotatePeaks.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.txt \
 mm10 -size 4000 -hist 10 -d ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_sample/ \
 ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_sample/ > ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.HIST.txt
annotatePeaks.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.txt \
 mm10 -size 4000 -hist 10 -d ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_sample/ \
 ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_sample/ > ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.HIST.txt

# DEG Peaks Heatmaps

annotatePeaks.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.txt \
 mm10 -size 4000 -hist 10 -ghist -d ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_sample/ \
 ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_sample/ > ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75DOWN_TGFbvsNormal.peaks.HEAT.txt
annotatePeaks.pl ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.txt \
 mm10 -size 4000 -hist 10 -ghist -d ${peak_dir}/HOMER/tagDirectory/Normal_H2azAc_sample/ \
 ${peak_dir}/HOMER/tagDirectory/TGFb_H2azAc_sample/ > ${peak_dir}/Consensus/sample_peaks.broadPeak/ChIP-Seq.pf.dds_p0.05_lgFC0.75UP_TGFbvsNormal.peaks.HEAT.txt
