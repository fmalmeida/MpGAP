process spades_hybrid {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag { x }
  cpus params.threads

  input:
  file(lreads)
  tuple val(id), file(sread1), file(sread2)
  file(sreads)
  val(prefix)

  output:
  file "*" // Save everything
  tuple file("spades/spades_assembly.fasta"), val(lrID), val('spades') // Gets contigs file

  script:
  // Check reads
  lr   = (params.lr_type == 'nanopore') ? '--nanopore' : '--pacbio'
  lrID = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads $lr $lreads"
    x = "Performing a hybrid assembly with SPAdes, using paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads $lr $lreads"
    x = "Performing a hybrid assembly with SPAdes, using single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2 $lr $lreads"
    x = "Performing a hybrid assembly with SPAdes, using paired end reads"
  }
  """
  spades.py -o spades -t ${params.threads} \\
  ${params.spades_additional_parameters} $parameter

  # Rename
  mv spades/contigs.fasta spades/spades_assembly.fasta
  """
}
