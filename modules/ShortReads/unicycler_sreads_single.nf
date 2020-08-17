process unicycler_sreads_single_assembly {
  publishDir "${params.outdir}/shortreads-only/paired_reads", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Executing Unicycler with single end reads"
  cpus params.threads

  input:
  file(sreads)

  output:
  file "unicycler_${id}" // Save everything
  tuple file("unicycler_${id}/assembly.fasta"), val(id), val('unicycler') // Gets contigs file

  script:
  id = sreads.getSimpleName()
  """
  unicycler -s $sreads -o unicycler_${id} -t ${params.threads} \\
  --no_correct ${params.unicycler_additional_parameters} &> unicycler.log
  """
}
