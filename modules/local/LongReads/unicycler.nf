process unicycler {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "unicycler/*" // Save all files
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file
  path('versions.yml'), emit: versions

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  additional_params = (params.unicycler_additional_parameters) ? params.unicycler_additional_parameters : ""
  if (corrected_longreads.toBoolean() || high_quality_longreads.toBoolean()) {
    if ( !additional_params.contains('--mode') ) { additional_params = additional_params + '--mode conservative' }
  }
  """  
  # run unicycler
  unicycler \\
      -l ${lreads} \\
      -o unicycler \\
      -t $task.cpus \\
      $additional_params

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      unicycler: \$( unicycler --version | cut -f 2 -d ' ' )
  END_VERSIONS
  """
}
