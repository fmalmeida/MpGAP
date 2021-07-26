process haslr_hybrid {
  publishDir "${params.outdir}/${lrID}/hybrid/strategy_1/${out_ids}", mode: 'copy'
  label 'main'
  tag { x }
  cpus params.threads

  input:
  file lreads
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "*" // Save everything
  tuple file("haslr/haslr_assembly.fa"), val(lrID), val('haslr') // Gets contigs file

  script:
  // Check reads
  lrID = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-s $sread1 $sread2 $sreads"
    x = "Performing a hybrid assembly with haslr, using paired and single end reads"
    srId = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    x = "Performing a hybrid assembly with haslr, using single end reads"
    id = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-s $sread1 $sread2"
    x = "Performing a hybrid assembly with haslr, using paired end reads"
    out_ids = "${id}"
  }

  """
  haslr.py -t ${params.threads} -o haslr -g ${params.genomeSize} \
  -l $lreads -x ${params.lr_type} ${parameter} \
  ${params.haslr_additional_parameters}

  # Rename
  cp haslr/*/asm.final.fa haslr/haslr_assembly.fa
  """
}
