process quast {
  publishDir "${params.outdir}/Quality_Assessment/quast_${type}", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Assessing ${assembler} assembly quality"

  input:
  tuple file(contigs), val(id), val(assembler)
  file sreads

  output:
  file "${assembler}/*"

  script:
  if (params.assembly_type == 'longreads-only') {
    type = 'lreadsOnly'
  } else if (params.assembly_type == 'illumina-only') {
    type = 'illuminaOnly'
  } else if (params.assembly_type == 'hybrid') {
    type = 'hybrid'
  }

  if (params.shortreads_paired && !params.shortreads_single) {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${sreads[1]} ${sreads[2]}"
    quast_parameter = "--pe1 ${sreads[1]} --pe2 ${sreads[2]}"
  } else if (!params.shortreads_paired && params.shortreads_single) {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${sreads}"
    quast_parameter = "--single ${sreads}"
  } else if (params.shortreads_paired && params.shortreads_single) {
    bwa_parameter   = "-M -t ${params.threads} ${contigs} ${sreads[1]} ${sreads[2]}"
    quast_parameter = "--single ${sreads[3]} --pe1 ${sreads[1]} --pe2 ${sreads[2]}"
  }

  """
  bwa index ${contigs} ;
  bwa mem ${bwa_parameter} > contigs_aln.sam ;
  quast.py -o ${assembler} -t ${params.threads} --sam contigs_aln.sam \\
  ${quast_parameter} --circos ${contigs}
  """
}
