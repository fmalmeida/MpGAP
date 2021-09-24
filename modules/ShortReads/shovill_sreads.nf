process shovill {
  publishDir "${params.outdir}/${prefix}/shovill", mode: 'copy'
  label 'main'
  tag "Shovill paired end assembly with ${assembler}"
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


// batch mode
process shovill_batch {
  publishDir "${params.outdir}/${prefix}/shovill", mode: 'copy'
  label 'main'
  tag "${id}: shovill with ${assembler}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), val(prefix), val(assembler)

  output:
  file "${assembler}" // Save all output
  tuple file("${assembler}/shovill_${assembler}_final.fasta"), val(id), val("shovill_${assembler}"), val(prefix), file(sread1), file(sread2), file(single) 

  when:
  (sread1 !=~ /input.?/ || sread2 !=~ /input.?/) && (single =~ /input.?/)

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
