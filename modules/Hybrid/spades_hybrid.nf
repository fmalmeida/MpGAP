process spades_hybrid {
  publishDir "${outdir}/hybrid", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus params.threads

  input:
  file lreads
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save everything
  tuple file("spades_${lreads.getSimpleName()}/contigs.fasta"), val("spades_${lreads.getSimpleName()}"), val('spades') // Gets contigs file

  script:
  lr = (params.lr_type == 'nanopore') ? '--nanopore' : '--pacbio'

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads $lr $lreads"
    x = "Executing assembly with paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads $lr $lreads"
    x = "Executing assembly with single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2 $lr $lreads"
    x = "Executing assembly with paired end reads"
  }
  """
  spades.py -o "spades_${lreads.getSimpleName()}" -t ${params.threads} \\
  ${params.spades_additional_parameters} $parameter
  """
}
