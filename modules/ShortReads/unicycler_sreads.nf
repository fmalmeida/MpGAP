process unicycler_sreads_assembly {
  publishDir "${params.outdir}/shortreads-only", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "unicycler_${id}" // Save everything
  tuple file("unicycler_${id}/assembly.fasta"), val(id), val('unicycler') // Gets contigs file

  script:

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads --no_correct"
    x = "Performing a illumina-only assembly with Unicycler, using paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads --no_correct"
    id = sreads.getSimpleName()
    x = "Performing a illumina-only assembly with Unicycler, using single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
    x = "Performing a illumina-only assembly with Unicycler, using paired end reads"
  }

  """
  unicycler $parameter -o unicycler_${id} -t ${params.threads} \\
  ${params.unicycler_additional_parameters} &> unicycler.log
  """
}
