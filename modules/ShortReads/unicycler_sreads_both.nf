process unicycler_sreads_both_assembly {
  publishDir "${params.outdir}/shortreads-only/paired_and_single_reads", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Executing Unicycler with both paired and single end reads"
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "unicycler_${id}" // Save everything
  tuple file("unicycler_${id}/assembly.fasta"), val(id), val('unicycler') // Gets contigs file

  script:
  """
  unicycler -1 $sread1 -2 $sread2 -s $sreads \\
  -o unicycler_${id} -t ${params.threads} \\
  ${params.unicycler_additional_parameters} &> unicycler.log
  """
}
