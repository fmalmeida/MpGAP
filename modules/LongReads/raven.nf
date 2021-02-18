process raven_assembly {
  publishDir "${params.outdir}/${lrID}/longreads_only", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with raven"

  input:
  file lreads

  output:
  file "raven" // Saves all files
  tuple file("raven/raven_contigs.fa"), val(lrID), val('raven') // Gets contigs file

  script:
  lrID = lreads.getSimpleName()
  """
  mkdir raven ;
  raven $lreads --threads ${params.threads} ${params.raven_additional_parameters} > raven/raven_contigs.fa ;
  """
}
