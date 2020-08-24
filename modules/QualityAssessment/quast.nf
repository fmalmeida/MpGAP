process quast {
  publishDir "${params.outdir}/Quality_Assessment/quast_${type}", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  tag "Assessing ${assembler} assembly quality"

  input:
  tuple file(contigs), val(id), val(assembler)
  file(reads)

  output:
  file "${assembler}/*"

  script:
  // Assembly Type - variable
  if (params.assembly_type == 'longreads-only') {
    type = 'lreadsOnly'
  } else if (params.assembly_type == 'illumina-only') {
    type = 'illuminaOnly'
  } else if (params.assembly_type == 'hybrid') {
    type = 'hybrid'
  }

  // Alignment parameters
  if (params.shortreads_paired && !params.shortreads_single) {
    quast_parameter = "--pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (!params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]}"
  } else if (params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]} --pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (params.assembly_type == 'longreads-only') {
    ltype           = (params.lr_type == 'nanopore') ? "ont2d" : "pacbio"
    quast_parameter = "--${params.lr_type} ${reads}"
  }

  """
  quast.py -o ${assembler} -t ${params.threads} ${quast_parameter} \\
  --circos --glimmer --rna-finding ${contigs}
  """
}
