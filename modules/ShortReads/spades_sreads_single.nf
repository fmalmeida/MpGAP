process spades_sreads_single_assembly {
  publishDir "${params.outdir}/shortreads-only/single_reads", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Executing SPAdes with single end reads"
  cpus params.threads

  input:
  file(sreads)

  output:
  file "*" // Save all output
  tuple file("spades_${id}/contigs.fasta"), val(id), val('spades') // Gets contigs file

  script:
  id = sreads.getSimpleName()
  """
  spades.py -o "spades_${id}" -t ${params.threads} \\
  ${params.spades_additional_parameters} -s $sreads
  """
}
