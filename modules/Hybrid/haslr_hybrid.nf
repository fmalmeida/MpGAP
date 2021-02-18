process haslr_hybrid {
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
  tuple file("haslr/*/asm.final.fa"), val(lrID), val('haslr') // Gets contigs file

  script:
  // Check reads
  lrID  = lreads.getSimpleName()
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-s $sread1 $sread2 $sreads"
    x = "Performing a hybrid assembly with haslr, using paired and single end reads"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    x = "Performing a hybrid assembly with haslr, using single end reads"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-s $sread1 $sread2"
  }

  """
  haslr.py -t ${params.threads} -o haslr -g ${params.genomeSize} \
  -l $lreads -x ${params.lr_type} ${parameter} \
  ${params.haslr_additional_parameters}
  """
}
