#!/bin/bash -x

set -euoE pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
DEFAULT=16
coreNum=${coreNum:-$DEFAULT}

if [ -z "${coreNum:-}" ]; then
  echo "coreNum was not set"
fi

samplerate=${samplerRate:-1}
samplereads=${sampleReads:-0}
# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#S3OUTPUTPATH = "${3}"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
#QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/SPAdes"

FASTQC_OUTPUT_1="${LOCAL_OUTPUT}/fastqc_raw"
FASTQC_OUTPUT_2="${LOCAL_OUTPUT}/fastqc_filtered"
ALIGNMENT_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
QC_FASTQ="${LOCAL_OUTPUT}/trimmed_fastq"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/Qualimap"
#FASTQ_NAME=${fastq1%/*}


mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${BAMQC_OUTPUT}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${FASTQC_OUTPUT_1}" "${FASTQC_OUTPUT_2}" "${ALIGNMENT_OUTPUT}"


hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT_1} \
"${RAW_FASTQ}/read1.fastq.gz" \
"${RAW_FASTQ}/read2.fastq.gz"

# Constant definitions for bbduk
adapterFile="adapters,phix"
trimQuality=${trimQuality:-30}
minLength=${minLength:-100}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
# Filter out reads that have an average entropy of under 0.5. A homopolymer such as GGGGGGGGGG would have entropy of zero
timem bbduk.sh -Xmx32g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 \
    in1="${RAW_FASTQ}/read1.fastq.gz" \
    in2="${RAW_FASTQ}/read2.fastq.gz" \
    out1="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    out2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" |\
    tee -a ${LOG_DIR}/bbduk.log


#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT_2} \
"${QC_FASTQ}/read1_trimmed.fastq.gz" \
"${QC_FASTQ}/read2_trimmed.fastq.gz"


# downsample 80% of trimmed reads
if [[ $samplerate -ne 1 ]]
then
  reformat.sh \
  minlength=145 \
  samplerate=${samplerate} \
  in="${QC_FASTQ}/read1_trimmed.fastq.gz" in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
  out="${QC_FASTQ}/read1_sampled.fastq.gz" out2="${QC_FASTQ}/read2_sampled.fastq.gz" |\
  tee -a ${LOG_DIR}/reformat.log.txt
elif [[ $samplereads -ne 0 ]]
then
  reformat.sh \
  minlength=145 \
  samplereadstarget=${samplereads} \
  in="${QC_FASTQ}/read1_trimmed.fastq.gz" in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
  out="${QC_FASTQ}/read1_sampled.fastq.gz" out2="${QC_FASTQ}/read2_sampled.fastq.gz" |\
  tee -a ${LOG_DIR}/reformat.log.txt
else
  reformat.sh \
  in="${QC_FASTQ}/read1_trimmed.fastq.gz" in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
  out="${QC_FASTQ}/read1_sampled.fastq.gz" out2="${QC_FASTQ}/read2_sampled.fastq.gz" |\
  tee -a ${LOG_DIR}/reformat.log.txt
fi



## Run SPAdes assembly
timem spades.py \
-t ${coreNum} \
-k 27,47,63,77,99,121 \
--careful \
--pe1-1 "${QC_FASTQ}/read1_sampled.fastq.gz" \
--pe1-2 "${QC_FASTQ}/read2_sampled.fastq.gz" \
-o ${ASSEMBLY_OUTPUT} |\
tee -a ${LOG_DIR}/spades_assembly.log

#filter assembly contigs by Length - at least 1000
#reformat.sh in="${ASSEMBLY_OUTPUT}/contigs.fasta" out="${ASSEMBLY_OUTPUT}/filtered_contigs_1kbp.fasta" minlength=1000

#BWA alignment reads to contigs
bwa index "${ASSEMBLY_OUTPUT}/contigs.fasta"
bwa mem -t ${coreNum} "${ASSEMBLY_OUTPUT}/contigs.fasta" "${QC_FASTQ}/read1_trimmed.fastq.gz" "${QC_FASTQ}/read2_trimmed.fastq.gz" | samtools sort -@ ${coreNum} -o "${ALIGNMENT_OUTPUT}/reads_aligned_contigs.sorted.bam" -

samtools index "${ALIGNMENT_OUTPUT}/reads_aligned_contigs.sorted.bam"
samtools flagstat "${ALIGNMENT_OUTPUT}/reads_aligned_contigs.sorted.bam" > "${LOG_DIR}/reads_aligned_contigs.flagstat.txt"


# get coverage and depth
pileup.sh in="${ALIGNMENT_OUTPUT}/reads_aligned_contigs.sorted.bam" out=${LOG_DIR}/contigs_coverage.txt overwrite=t 2>${LOG_DIR}/contigs_stats.txt

# bam qc
qualimap bamqc -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${ALIGNMENT_OUTPUT}/reads_aligned_contigs.sorted.bam" -c


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
