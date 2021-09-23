process shovill_sreads_assembly {
  publishDir "${params.outdir}/${prefix}/shovill", mode: 'copy'
  label 'main'
  tag "Performing a illumina-only assembly with shovill (${assembler}), using paired end reads"
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2), val(prefix), val(assembler)

  output:
  file "${assembler}" // Save all output
  tuple file("${assembler}/shovill_${assembler}_final.fasta"), val(id), val("shovill_${assembler}") // Gets contigs file

  when:
  ((params.shortreads_paired) && (!params.shortreads_single))

  script:
  """
  # Activate env
  source activate shovill ;

  # Run
  shovill \
    --outdir ${assembler} \
    --assembler ${assembler} \
    --R1 $sread1 --R2 $sread2 \
    --cpus ${params.threads} \
    --trim \
    ${params.shovill_additional_parameters}

  # Rename assembly
  mv ${assembler}/contigs.fa ${assembler}/shovill_${assembler}_final.fasta
  """
}
