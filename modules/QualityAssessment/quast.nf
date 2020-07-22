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
  if (params.shortreads_paired && !params.shortreads_single && params.assembly_type == 'illumina-only') {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${reads[1]} ${reads[2]}"
    quast_parameter = "--pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (!params.shortreads_paired && params.shortreads_single && params.assembly_type == 'illumina-only') {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${reads}"
    quast_parameter = "--single ${reads}"
  } else if (params.shortreads_paired && params.shortreads_single && params.assembly_type == 'illumina-only') {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${reads[1]} ${reads[2]}"
    quast_parameter = "--single ${reads[3]} --pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (params.assembly_type == 'longreads-only') {
    ltype           = (params.lr_type == 'nanopore') ? "ont2d" : "pacbio"
    bwa_parameter   = "-x ${ltype} -t ${params.threads} ${contigs} ${reads}"
    quast_parameter = "--${params.lr_type} ${reads}"
  }

  """
  /work/quast/quast.py -o ${assembler} -t ${params.threads} ${quast_parameter} \\
  --circos --glimmer --rna-finding ${contigs}
  """
}
