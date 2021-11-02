process unicycler {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), val(shasta_config), file(bams), val(prefix)

  output:
  file "unicycler/*" // Save all files
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
