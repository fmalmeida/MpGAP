process unicycler_hybrid {
  publishDir "${params.outdir}/${lrID}/hybrid", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  tag { x }
  cpus params.threads

  input:
  file lreads
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save everything
  tuple file("unicycler/assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  // Check reads
  lrID  = lreads.getSimpleName()
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
  -o unicycler -t ${params.threads} \\
  ${params.unicycler_additional_parameters} &>unicycler.log
  """
}
