#!/bin/bash -x

set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

LOCAL=$(pwd)
coreNum=${coreNum:-16}
ShortTargetBase=${ShortTargetBase:-500000000}
LongTargetBase=${LongTargetBase:-500000000}
genomeSize=${genomeSize:-5000000}
coverage=${coverage:-400}


# s3 inputs from env variables
#fastq1="${1}"
#fastq2="${2}"
#longreads="${3}" input long reads in fastq format
#S3OUTPUTPATH="${4}"
#genomeSize="${5}" genome size in bp


# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"


LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
ASSEMBLY_OUTPUT="${LOCAL_OUTPUT}/UNICYCLER"
QUAST_OUTPUT="${LOCAL_OUTPUT}/quast"
FASTQC_OUTPUT="${LOCAL_OUTPUT}/fastqc_illumina"
FASTQC_OUTPUT2="${LOCAL_OUTPUT}/fastqc_nanopore"
RAW_NANOPORE="${LOCAL_OUTPUT}/nanopore"
#FASTQ_NAME=${fastq1%/*}
SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

LOCAL_DB_PATH=${LOCAL}/databases
LOCAL_DB_NAME="contigs"
OUTPUT_PREFIX="trimmed_reads_vs_${LOCAL_DB_NAME}"
BAMQC_OUTPUT="${LOCAL_OUTPUT}/bamqc_reads_vs_contigs"
BAMQC_OUTPUT2="${LOCAL_OUTPUT}/bamqc_longreads_vs_contigs"
ALIGN_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
POST_OUTPUT="${LOCAL_OUTPUT}/post-processing"
ANNO_OUTPUT="${LOCAL_OUTPUT}/annotation"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${POST_OUTPUT}"
mkdir -p "${ASSEMBLY_OUTPUT}" "${QUAST_OUTPUT}" "${FASTQC_OUTPUT}" "${FASTQC_OUTPUT2}" "${ANNO_OUTPUT}"
mkdir -p "${LOCAL_DB_PATH}" "${ALIGN_OUTPUT}"  "${BAMQC_OUTPUT}" "${BAMQC_OUTPUT2}" "${RAW_NANOPORE}"
trap '{ rm -rf ${OUTPUTDIR} ; exit 255; }' 1

hash_kmer=${hash_kmer:-51}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
aws s3 cp --quiet ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"
aws s3 cp --quiet ${longreads} "${RAW_FASTQ}/long.fasta"

# count the number of PE reads
#zcat "${RAW_FASTQ}/read1.fastq.gz" | echo $((`wc -l`/2)) > ${LOG_DIR}/illumina_count.txt
# count the long reads
#zcat "${RAW_FASTQ}/long.fastq" | echo $((`wc -l`/4)) > ${LOG_DIR}/nanopore_count.txt


## Collate called bases into a single fastq file
#aws s3 sync --quiet ${longreads} ${RAW_NANOPORE}
#find "${RAW_NANOPORE}/*.fastq" | xargs cat > "${RAW_FASTQ}/long.fastq"

#cp  ${fastq1} "${RAW_FASTQ}/read1.fastq.gz"
#cp  ${fastq2} "${RAW_FASTQ}/read2.fastq.gz"


###############################################################
# Pre-processing before assembly
echo "Reads Pre-processing before alignment"
###############################################################


# Constant definitions for bbduk, increased the quailty
adapterFile="adapters,phix"
trimQuality=${trimQuality:-25} #old 25
minLength=${minLength:-50}  #old 50
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

# Use bbduk to trim short reads, -eoom exits when out of memory
#in1="${QC_FASTQ}/read1_deduped.fastq.gz" \
#in2="${QC_FASTQ}/read2_deduped.fastq.gz" \
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
echo "Reads trimming/filtering" >> ${LOG_DIR}/bbtools.log.txt
echo "**************************" >> ${LOG_DIR}/bbtools.log.txt
bbduk.sh -Xmx16g tbo -eoom hdist=1 qtrim=rl ktrim=r \
    entropy=0.5 entropywindow=50 entropyk=5 entropytrim=rl \
    in="${RAW_FASTQ}/read1.fastq.gz" \
    in2="${RAW_FASTQ}/read2.fastq.gz" \
    out1="${QC_FASTQ}/read1_trimmed.fastq.gz" \
    out2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
    ref=${adapterFile} \
    k="${kmer_value}" \
    mink="${min_kmer_value}" \
    trimq="${trimQuality}" \
    minlen="${minLength}" \
    tossbrokenreads=t \
    refstats="${LOCAL_OUTPUT}/BBDuk/adapter_trimming_stats_per_ref.txt" \
    >> ${LOG_DIR}/bbtools.log.txt 2>&1

#Run fastqc for short reads
fastqc \
-t ${coreNum} \
-o ${FASTQC_OUTPUT} \
"${QC_FASTQ}/read1_trimmed.fastq.gz" \
"${QC_FASTQ}/read2_trimmed.fastq.gz"


# Run mapping
bbmap.sh ambiguous=all tossbrokenreads=t \
in="${QC_FASTQ}/read1_trimmed.fastq.gz"  \
in2="${QC_FASTQ}/read2_trimmed.fastq.gz" \
ref="${RAW_FASTQ}/long.fasta" out=stdout | \
samtools view -@ ${coreNum} -bh -o "${ALIGN_OUTPUT}/contigs_aligned.bam" - | \
tee -a ${LOG_DIR}/bbmap_mapping.log.txt

#-S "${ALIGN_OUTPUT}/contigs_aligned.sam"
# index contigs
#samtools faidx "${ASSEMBLY_OUTPUT}/assembly.fasta"
# create bam file
#samtools import "${ASSEMBLY_OUTPUT}/assembly.fasta.fai" "${ALIGN_OUTPUT}/contigs_aligned.sam" "${ALIGN_OUTPUT}/contigs_aligned.bam"
# sort bam file
samtools sort -@ ${coreNum} "${ALIGN_OUTPUT}/contigs_aligned.bam"  -o "${ALIGN_OUTPUT}/contigs_aligned.sorted.bam"
# index bam
samtools index  "${ALIGN_OUTPUT}/contigs_aligned.sorted.bam"
# generate contigs stats
samtools idxstats "${ALIGN_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.idxstats.txt"
samtools flagstat "${ALIGN_OUTPUT}/contigs_aligned.sorted.bam" > "${LOG_DIR}/contigs_aligned.flagstat.txt"

#bamqc
qualimap bamqc --java-mem-size=16G -nt ${coreNum} -outdir "${BAMQC_OUTPUT}" -bam "${ALIGN_OUTPUT}/contigs_aligned.sorted.bam" -c

grep "mean coverageData" ${BAMQC_OUTPUT}/genome_results.txt | cut -f2- -d"=" | sed s/X// > "${LOG_DIR}/coverage_depth.txt"


######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
# rm -rf "${OUTPUTDIR}"
