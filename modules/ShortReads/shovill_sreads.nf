// batch mode
process shovill {
  publishDir "${params.output}/${prefix}/shovill", mode: 'copy'
  label 'main'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix), val(assembler)

  output:
  file "${assembler}" // Save all output
  tuple val(id), file("${assembler}/shovill_${assembler}_final.fasta"), val("shovill_${assembler}")

  when:
  !(sread1 =~ /input.*/ || sread2 =~ /input.*/) && (single =~ /input.*/) && (entrypoint == 'shortreads_only')

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
