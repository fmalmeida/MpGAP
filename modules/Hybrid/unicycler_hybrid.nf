process unicycler_hybrid {
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
  tuple file("unicycler/unicycler_assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  // Check reads
  lrID = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads -l $lreads"
    x = "Performing a hybrid assembly with Unicycler, using paired and single end reads"
    srId = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads -l $lreads"
    x = "Performing a hybrid assembly with Unicycler, using single end reads"
    id = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2 -l $lreads"
    x = "Performing a hybrid assembly with Unicycler, using paired end reads"
    out_ids = "${id}"
  }

  """
  unicycler ${parameter} \\
  -o unicycler -t ${params.threads} \\
  ${params.unicycler_additional_parameters}

  # Rename
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
