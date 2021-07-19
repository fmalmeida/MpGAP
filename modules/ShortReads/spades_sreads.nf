process spades_sreads_assembly {
  publishDir "${params.outdir}/${out_ids}/shortreads_only", mode: 'copy'
  label 'main'
  tag { x }
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)
  file(sreads)

  output:
  file "spades" // Save all output
  tuple file("spades/spades_assembly.fasta"), val(out_ids), val('spades') // Gets contigs file

  script:

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sreads"
    x = "Performing a illumina-only assembly with SPAdes, using paired and single end reads"
    srId = (sreads - ".gz")[0].getBaseName()
    out_ids = "${id}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    parameter = "-s $sreads"
    id = (sreads - ".gz")[0].getBaseName()
    x = "Performing a illumina-only assembly with SPAdes, using single end reads"
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    parameter = "-1 $sread1 -2 $sread2"
    x = "Performing a illumina-only assembly with SPAdes, using paired end reads"
    out_ids = "${id}"
  }

  """
  spades.py -o spades -t ${params.threads} \\
  ${params.spades_additional_parameters} $parameter

  # Rename assembly
  mv spades/contigs.fasta spades/spades_assembly.fasta
  """
}
