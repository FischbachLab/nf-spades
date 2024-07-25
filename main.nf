#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run Spades metagenome assembly pipeline

    Required Arguments:
      --seedfile      file      a file containing sample name, reads1 and reads2
      --sampleRate    num       sampling rate (0-1) [1]
      --sampleReads   num       exact number of output reads (or pairs) desired [0].
      --outdir   path      output a s3 path

    Options:
      -profile        docker run locally

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}


/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

//core = params.
//mem =
//srate=

//def output_path = "${params.output_path}"
//def output_path = "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/${params.output_prefix}"

//println output_path
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunksize'.
 * Finally assign the result channel to the variable 'fasta_ch'
 */



/*
 * Run SPAdes assembly
 /SPAdes/contigs.fasta
 publishDir "${output_path}", mode:'copy'

 */
process spades_assembly {

    tag "$sample"
    container params.container

    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }

    input:
	  tuple val(sample), val(reads1), val(reads2)

    output:
    //path "tmp_*/Sync/SPAdes/contigs.fasta" into spades_ch

    script:
    """
    export fastq1="${reads1}"
    export fastq2="${reads2}"
    export sampleRate="${params.sampleRate}"
    export sampleReads="${params.sampleReads}"
    export S3OUTPUTPATH="${params.outdir}/${sample}"
    bash run_spades.sh
    """
}



workflow {

  seedfile_ch = Channel
  .fromPath(params.seedfile)
  .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
  .splitCsv(header: ['sample', 'reads1', 'reads2'], sep: ',', skip: 1)
  .map{ row -> tuple(row.sample, row.reads1, row.reads2)}


  seedfile_ch | spades_assembly

}

/*
Run Quast assessment

process QUAST {

  cpus 2
  memory 8.GB
  container params.docker_container_quast

  publishDir "${output_path}/${sample}/Quast", mode:'copy'

  input:
    path assembly from spades_ch

  output:
    path "Quast/*"

  script:
  """
  quast.py ${assembly} -o Quast
  """
}
*/
