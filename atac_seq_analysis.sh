#fastp------------------------------------
#!/bin/bash

# Set the number of threads
THREADS=10

# Input and output directories
INPUT_DIR="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/Rawdata"
OUTPUT_DIR="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/1_fastp"

# List of sample folders
SAMPLES=("neu0_1" "neu0_2" "neu100_1" "neu100_2")

# Create output directory (if it does not exist)
mkdir -p "$OUTPUT_DIR"

# Iterate over each sample folder
for SAMPLE in "${SAMPLES[@]}"; do
    # Find R1 and R2 files in the sample folder
    R1_FILE=$(find "$INPUT_DIR/$SAMPLE" -name "*_R1*.fastq.gz")
    R2_FILE=$(find "$INPUT_DIR/$SAMPLE" -name "*_R2*.fastq.gz")

    # Check if files exist
    if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
        echo "Error: FASTQ files for $SAMPLE not found."
        continue
    fi

    # Run fastp for quality control
    fastp -i "$R1_FILE" -I "$R2_FILE" -o "${OUTPUT_DIR}/${SAMPLE}_R1.trimmed.fastq.gz" -O "${OUTPUT_DIR}/${SAMPLE}_R2.trimmed.fastq.gz" --detect_adapter_for_pe -j "${OUTPUT_DIR}/${SAMPLE}.fastp.json" -h "${OUTPUT_DIR}/${SAMPLE}.fastp.html" --thread $THREADS
done



#================================= bwa-mem2===============================
#align----------------------------
#!/bin/bash

# Set the number of threads
THREADS=10

# Input and output directories
FASTQ_DIR="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/1_fastp"
OUTPUT_DIR="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/2_bwa_mem2"
INDEX="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_zhouxuan_250116/3_bwa_mem2/GRCm39.genome.fa"

# List of samples
SAMPLES=("neu0_1" "neu0_2" "neu100_1" "neu100_2")

# Create output directory (if it does not exist)
mkdir -p "$OUTPUT_DIR"

# Iterate over each sample
for SAMPLE in "${SAMPLES[@]}"; do
  # Define input files and output file
  R1_FILE="${FASTQ_DIR}/${SAMPLE}_R1.trimmed.fastq.gz"
  R2_FILE="${FASTQ_DIR}/${SAMPLE}_R2.trimmed.fastq.gz"
  OUTPUT_FILE="${OUTPUT_DIR}/${SAMPLE}.sam"

  # Check if input files exist
  if [[ ! -f "$R1_FILE" || ! -f "$R2_FILE" ]]; then
    echo "Error: FASTQ files for $SAMPLE not found."
    continue
  fi

  # Run bwa-mem2 for alignment
  bwa-mem2 mem -t $THREADS "$INDEX" "$R1_FILE" "$R2_FILE" > "$OUTPUT_FILE"
done



conda activate ac4c
#===================================================================
#!/bin/bash

# Set input and output directories
input_dir="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/2_bwa_mem2"
output_dir="${input_dir}/rmMT"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Iterate over all .sam files in the input directory
for sam_file in "${input_dir}"/*.sam; do
  # Get the file name without path and extension
  base_name=$(basename "${sam_file}" .sam)
  
  # Generate intermediate bam file name
  sorted_bam="${input_dir}/${base_name}_sorted.bam"
  
  # Convert sam file to bam file and sort
  samtools view -S -b "${sam_file}" | samtools sort -o "${sorted_bam}"
  
  # Remove mitochondrial aligned reads and generate new bam file
  samtools view -h "${sorted_bam}" | grep -v chrM | samtools sort -O bam -o "${output_dir}/${base_name}.rmChrM.bam" -T .
  
  # Delete intermediate bam file
  rm "${sorted_bam}"
done

echo "Processing complete, all output files are saved to ${output_dir}"




#========================================#Mark duplicates
#!/bin/bash

# Set input and output directories
input_dir="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/2_bwa_mem2/rmMT"
output_dir="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/3_dedup"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Iterate over all .bam files in the input directory
for bam_file in "${input_dir}"/*.rmChrM.bam; do
  # Get the file name without path and extension
  base_name=$(basename "${bam_file}" .rmChrM.bam)
  
  # Define intermediate and output file names
  marked_bam="${output_dir}/${base_name}.marked.bam"
  metrics_file="${output_dir}/${base_name}.dup.metrics"
  filtered_bam="${output_dir}/${base_name}.filtered.bam"

  # Mark duplicate reads
  picard MarkDuplicates QUIET=true INPUT="${bam_file}" OUTPUT="${marked_bam}" METRICS_FILE="${metrics_file}" REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=.

  # View duplication rate
  echo "Duplicate metrics for ${base_name}:"
  head -n 8 "${metrics_file}" | cut -f 7,9 | grep -v ^# | tail -n 2

  # Filter reads
  samtools view -h -b -f 2 -F 1548 -q 30 "${marked_bam}" | samtools sort -o "${filtered_bam}"

  # Index the filtered BAM file
  samtools index "${filtered_bam}"
done

echo "Processing complete, all output files are saved to ${output_dir}"






conda activate macs2
#==========================
#!/bin/bash

# Set input and output directories
input_dir="/data02/tangyk/scrna/Prostate_Cancer/otherdata/Neutrophils/atac_seq_yadan_250120/3_dedup"
output_dir="/home01/tangyk/scrna/Prostate_Cancer/mouse_data/4_atac_seq_yadan_250120/4_macs2"

# Create output directory if it does not exist
mkdir -p "${output_dir}"

# Iterate over all .filtered.bam files in the input directory
for bam_file in "${input_dir}"/*.filtered.bam; do
  # Get the file name without path and extension
  base_name=$(basename "${bam_file}" .filtered.bam)
  
  # Define output directory
  sample_output_dir="${output_dir}/${base_name}"
  mkdir -p "${sample_output_dir}"
  
  # Run MACS2 for peak calling
  macs2 callpeak -t "${bam_file}" -n "${base_name}" --shift -100 --extsize 200 --nomodel -B --SPMR -g hs --outdir "${sample_output_dir}" 2> "${sample_output_dir}/${base_name}.macs2.log"
  
  echo "Peak calling completed for ${base_name}, results saved in ${sample_output_dir}"
done

echo "Peak calling for all samples is complete, results are saved to ${output_dir}"
