process canu_assembly {
  publishDir "${params.outdir}/${lrID}/${type}", mode: 'copy'
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

  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2'
  }
  """
  canu -p ${lrID} -d canu maxThreads=${params.threads}\
  genomeSize=${params.genomeSize} ${params.canu_additional_parameters} $lr $lreads

  # Rename contigs
  mv canu/${lrID}.contigs.fasta canu/canu_${lrID}_contigs.fasta
  """
}
