process unicycler_hybrid {
  publishDir "${params.outdir}/hybrid", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus params.threads

  input:
  file lreads
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save everything
  tuple file("unicycler_${lreads.getSimpleName()}/assembly.fasta"), val("unicycler_${lreads.getSimpleName()}"), val('unicycler') // Gets contigs file

  script:

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads -l $lreads --no_correct"
    x = "Performing a hybrid assembly with Unicycler, using paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads -l $lreads --no_correct"
    x = "Performing a hybrid assembly with Unicycler, using single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2 -l $lreads"
    x = "Performing a hybrid assembly with Unicycler, using paired end reads"
  }

  """
  unicycler ${parameter} \\
  -o unicycler_${lreads.getSimpleName()} -t ${params.threads} \\
  ${params.unicycler_additional_parameters} &>unicycler.log
  """
}
