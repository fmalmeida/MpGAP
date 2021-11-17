// batch mode
process shovill {
  publishDir "${params.output}/${prefix}/shovill", mode: 'copy'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix), val(assembler)

  output:
  file "${assembler}" // Save all output
  tuple val(id), file("${assembler}/shovill_${assembler}_final.fasta"), val("shovill_${assembler}")

  when:
  !(sread1 =~ /input.*/ || sread2 =~ /input.*/) && (single =~ /input.*/) && (entrypoint == 'shortreads_only')

  script:
  """
  # run shovill
  shovill \\
      --outdir ${assembler} \\
      --assembler ${assembler} \\
      --R1 $sread1 \\
      --R2 $sread2 \\
      --cpus ${params.threads} \\
      ${params.shovill_additional_parameters} \\
      --trim ;

  # rename results
  mv ${assembler}/contigs.fa ${assembler}/shovill_${assembler}_final.fasta
  """
}
