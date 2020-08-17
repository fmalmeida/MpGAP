process spades_sreads_assembly {
  publishDir "${params.outdir}/shortreads-only", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save all output
  tuple file("spades_${id}/contigs.fasta"), val(id), val('spades') // Gets contigs file

  script:

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads"
    x = "Performing a illumina-only assembly with SPAdes, using paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    id = sreads.getSimpleName()
    x = "Performing a illumina-only assembly with SPAdes, using single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
    x = "Performing a illumina-only assembly with SPAdes, using paired end reads"
  }

  """
  spades.py -o "spades_${id}" -t ${params.threads} \\
  ${params.spades_additional_parameters} $parameter
  """
}
