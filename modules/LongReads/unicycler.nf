process unicycler {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "unicycler/*" // Save all files
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  """
  # run unicycler
  unicycler \\
      -l ${lreads} \\
      -o unicycler \\
      -t ${params.threads} \\
      ${params.unicycler_additional_parameters} \\
      --spades_path spades-3.13.0.py

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
