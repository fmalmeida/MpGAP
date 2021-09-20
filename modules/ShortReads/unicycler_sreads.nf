process unicycler_sreads_assembly {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag { x }
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2), file(sreads), val(prefix)

  output:
  file "unicycler" // Save everything
  tuple file("unicycler/unicycler_assembly.fasta"), val(out_ids), val('unicycler') // Gets contigs file

  script:
  // Check reads
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads"
    x = "Performing a illumina-only assembly with Unicycler, using paired and single end reads"
    srId = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    id = (sreads.getName() - ".gz").toString().substring(0, (sreads.getName() - ".gz").toString().lastIndexOf("."))
    x = "Performing a illumina-only assembly with Unicycler, using single end reads"
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
    x = "Performing a illumina-only assembly with Unicycler, using paired end reads"
    out_ids = "${id}"
  }

  """
  unicycler $parameter -o unicycler -t ${params.threads} ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
