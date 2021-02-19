process quast {
  publishDir "${params.outdir}/${id}/quality_assessment/${type}", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Assessing ${assembler} assembly quality"

  input:
  tuple file(contigs), val(id), val(assembler), file(reads)

  output:
  file "${assembler}/*"
  file("${assembler}")
  val(id)

  script:
  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads-only'
  } else if ((params.shortreads_paired || params.shortreads_single) && !params.longreads ) {
    type = 'shortreads-only'
  } else {
    type = 'hybrid'
  }

  // Alignment parameters
  if (params.shortreads_paired && !params.shortreads_single) {
    quast_parameter = "--pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (!params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]}"
  } else if (params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]} --pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else {
    ltype           = (params.lr_type == 'nanopore') ? "ont2d" : "pacbio"
    quast_parameter = "--${params.lr_type} ${reads}"
  }

  """
  quast.py -o ${assembler} -t ${params.threads} ${quast_parameter} \\
  --circos --glimmer --rna-finding ${contigs}
  """
}
