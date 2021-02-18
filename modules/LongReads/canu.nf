process canu_assembly {
  publishDir "${params.outdir}/${lrID}/longreads_only", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with Canu"

  input:
  file lreads

  output:
  file "canu/" // Saves all files
  tuple file("canu/*.contigs.fasta"), val(lrID), val('canu') // Gets contigs file

  script:
  lr = (params.lr_type == 'nanopore') ? '-nanopore-raw' : '-pacbio-raw'
  lrID = lreads.getSimpleName()
  """
  canu -p ${lrID} -d canu maxThreads=${params.threads}\
  genomeSize=${params.genomeSize} ${params.canu_additional_parameters} $lr $lreads
  """
}
