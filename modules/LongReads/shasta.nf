process shasta_assembly {
  publishDir "${params.outdir}/${lrID}/${type}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with shasta"

  input:
  file lreads

  output:
  file "shasta/" // Saves all files
  tuple file("shasta/shasta_assembly.fasta"), val(lrID), val('shasta') // Gets contigs file

  script:
  lr        = (params.lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  in_reads  = (lreads.getName() - ".gz")
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  # unzip reads
  gunzip -d -f $lreads

  # assemble
  shasta --assemblyDirectory shasta --threads ${params.threads} ${params.shasta_additional_parameters} --input $in_reads

  # Rename contigs
  cp shasta/Assembly.fasta shasta/shasta_assembly.fasta
  """
}
