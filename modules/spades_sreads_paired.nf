process spades_sreads_paired_assembly {
  publishDir "${params.outdir}/shortreads-only/paired_reads", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Executing SPAdes with paired end reads"
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)

  output:
  tuple file("spades_${id}/contigs.fasta"), val(id), val('spades') // Gets contigs file
  file "*" // Save all output

  script:
  """
  spades.py -o "spades_${id}" -t ${params.threads} \\
  ${params.spades_additional_parameters} -1 $sread1 -2 $sread2
  """
}
