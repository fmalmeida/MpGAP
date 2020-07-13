process spades_sreads_both_assembly {
  publishDir "${params.outdir}/shortreads-only/paired_and_single_reads", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Executing SPAdes with both paired and single end reads"
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save all output
  tuple file("spades_${id}/contigs.fasta"), val(id), val('spades') // Gets contigs file

  script:
  """
  spades.py -o "spades_${id}" -t ${params.threads} \\
  ${params.spades_additional_parameters} -1 $sread1 -2 $sread2 -s $sreads
  """
}
