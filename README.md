README
====================

This short-read de novo assembly pipeline first filters and trims reads based on quialty, then assembles filtered reads using SPAdes, a popular open-source assembly tool for genome reconstruction from DNA sequencing reads.   


## The example batch script showing the Spades assembly example for the Nextflow framework with a seedfile as the input file.

<!--
```{bash}
nextflow run -resume main.nf --sampleReads 2000 --seedfile 's3://genomics-workflow-core/Results/Spades/seedfile.csv' --output_path 's3://genomics-workflow-core/Results/Spades/Test' -profile docker
```
-->

```{bash}
aws batch submit-job \
    --job-name nf-spades \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="https://github.com/FischbachLab/nf-spades, \
"--seedfile", "s3://genomics-workflow-core/Results/Spades/seedfile.csv", \
"--output_path", "s3://genomics-workflow-core/Results/Spades/Project_Name" "
```

## A seedfile example
### Note that the seedfile is a CSV (comma-separated values) file with header
### The format of the seedfile is sample_name,short_R1,short_R2
```{bash}
sampleName,R1,R2
20240411_TN-DB488_F03,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/20240411_TN-DB488_F03_R1.fastq.gz,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/2024041
1_TN-DB488_F03_R2.fastq.gz
20240411_TN-DB490_G03,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/20240411_TN-DB490_G03_R1.fastq.gz,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/2024041
1_TN-DB490_G03_R2.fastq.gz
20240411_TN-DB491_H03,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/20240411_TN-DB491_H03_R1.fastq.gz,s3://genomics-workflow-core/Results/Basespace/NextSeq/20240709_M101-WGS_HKT5LBGXW/2024041
1_TN-DB491_H03_R2.fastq.gz
```

## The final filtered contigs (minlength>1kb) file can be found at
```{bash}
s3://genomics-workflow-core/Results/Spades/Project_Name/Sample_Name/SPAdes/filtered_contigs.fasta
```