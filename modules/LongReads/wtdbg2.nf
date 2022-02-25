process wtdbg2 {
  publishDir "${params.output}/${prefix}/wtdbg2", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Saves all files
  tuple val(id), file("wtdbg2_assembly.fasta"), val('wtdbg2') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  fixed_id = id - ":strategy_2"
  """
  # run wtdbg2
  wtdbg2.pl \\
      -t $task.cpus \\
      -x ${wtdbg2_technology} \\
      -g ${genome_size} \\
      -o ${fixed_id} \\
      ${params.wtdbg2_additional_parameters} \\
      $lreads

  # rename results
  cp ${fixed_id}.cns.fa wtdbg2_assembly.fasta
  """
}
