process wtdbg2_assembly {
  publishDir "${params.outdir}/${lrID}/${type}/wtdbg2", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with wtdbg2"

  input:
  file lreads

  output:
  file "*" // Saves all files
  tuple file("wtdbg2_assembly.fasta"), val(lrID), val('wtdbg2') // Gets contigs file

  script:
  lr                = (params.lr_type == 'nanopore') ? 'ont' : params.wtdbg2_technology
  lrID              = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  wtdbg2.pl -t ${params.threads} -x $lr -g ${params.genomeSize} -o ${lrID} ${params.wtdbg2_additional_parameters} $lreads

  # Rename contigs
  cp ${lrID}.cns.fa wtdbg2_assembly.fasta
  """
}
