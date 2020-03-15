process canu_assembly {
  publishDir "${params.outdir}/longreads-only", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  file lreads

  output:
  file "canu_${lrID}/" // Saves all files
  tuple file("canu_${lrID}/*.contigs.fasta"), val(lrID), val('canu') // Gets contigs file

  script:
  lr = (params.lr_type == 'nanopore') ? '-nanopore-raw' : '-pacbio-raw'
  lrID = lreads.getSimpleName()
  """
  canu -p ${lrID} -d canu_${lrID} maxThreads=${params.threads}\
  genomeSize=${params.genomeSize} ${params.canu_additional_parameters} $lr $lreads
  """
}
