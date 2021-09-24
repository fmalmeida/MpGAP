process quast {
  publishDir "${params.outdir}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Assessing ${assembler} assembly quality for multiqc"

  input:
  tuple file(contigs), val(id), val(assembler), file(reads), val(prefix)

  output:
  file("${assembler}")

  script:
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
  --conserved-genes-finding --rna-finding --min-contig 100 \\
  ${params.quast_additional_parameters} \\
  ${contigs}
  """
}

// batch mode
process quast_batch {
  publishDir "${params.outdir}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Assessing ${assembler} assembly quality for multiqc"

  input:
  tuple file(contigs), val(id), val(assembler), val(prefix), file(sread1), file(sread2), file(single) 

  output:
  tuple file("${assembler}"), val(prefix)

  script:
  // Alignment parameters
  paired_param = (sread1 !=~ /input.?/ || sread2 !=~ /input.?/) ? "--pe1 ${sread1} --pe2 ${sread2}" : ""
  single_param = (single !=~ /input.?/) ? "--single ${single}" : ""

  """
  quast.py -o ${assembler} -t ${params.threads} ${paired_param} ${single_param} \\
  --conserved-genes-finding --rna-finding --min-contig 100 \\
  ${params.quast_additional_parameters} \\
  ${contigs}
  """
}
