process quast {
  publishDir "${params.output}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Assessing ${assembler} assembly quality for multiqc"

  input:
  tuple val(id), file(contigs), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix)

  output:
  tuple val(id), val(entrypoint), val(prefix), file("${assembler}")

  script:
  // Alignment parameters
  paired_param = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "--pe1 ${sread1} --pe2 ${sread2}" : ""
  single_param = !(single =~ /input.?/) ? "--single ${single}" : ""
  lreads_param = !(lreads =~ /input.?/) ? "--${lr_type} ${lreads}" : ""

  """
  quast.py -o ${assembler} -t ${params.threads} \\
  ${lreads_param} ${paired_param} ${single_param} \\
  --conserved-genes-finding --rna-finding --min-contig 100 \\
  ${params.quast_additional_parameters} \\
  ${contigs}
  """
}
