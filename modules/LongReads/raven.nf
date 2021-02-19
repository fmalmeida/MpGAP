process raven_assembly {
  publishDir "${params.outdir}/${lrID}/longreads_only/raven", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with raven"

  input:
  file lreads

  output:
  file "raven_contigs.fa" // Saves all files
  tuple file("raven_contigs.fa"), val(lrID), val('raven') // Gets contigs file

  script:
  lrID = lreads.getSimpleName()
  """
  raven $lreads --threads ${params.threads} ${params.raven_additional_parameters} > raven_contigs.fa ;
  """
}
