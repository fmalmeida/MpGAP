process wtdbg2 {
  publishDir "${params.output}/${prefix}/wtdbg2", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Saves all files
  tuple val(id), file("wtdbg2_assembly.fasta"), val('wtdbg2') // Gets contigs file
  path('versions.yml'), emit: versions

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  fixed_id = id - ":strategy_2"
  additional_params = (params.wtdbg2_additional_parameters) ? params.wtdbg2_additional_parameters : ""
  """
  # run wtdbg2
  wtdbg2.pl \\
      -t $task.cpus \\
      -x ${wtdbg2_technology} \\
      -g ${genome_size} \\
      -o ${fixed_id} \\
      $additional_params \\
      $lreads

  # rename results
  cp ${fixed_id}.cns.fa wtdbg2_assembly.fasta

  # get version
  # --version command does not exist
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      wtdbg2: 2.5
  END_VERSIONS
  """
}
