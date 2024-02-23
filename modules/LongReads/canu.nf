process canu {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "canu/" // Saves all files
  tuple val(id), file("canu/canu_assembly.fasta"), val('canu') // Gets contigs file
  path('versions.yml')

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  corrected = (corrected_longreads.toBoolean() || high_quality_longreads.toBoolean()) ? '-corrected' : '' // canu does not have a specific config for high-quality, however, --corrected means, skipping canu read correction phase, which is what we want.
  fixed_id = id - ":strategy_2"
  additional_params = (params.canu_additional_parameters) ? params.canu_additional_parameters : ""
  """
  # run canu
  canu \\
      -p ${fixed_id} \\
      -d canu \\
      maxThreads=$task.cpus \\
      genomeSize=${genome_size} \\
      $additional_params \\
      $corrected \\
      $lr $lreads

  # rename results
  mv canu/${fixed_id}.contigs.fasta canu/canu_assembly.fasta

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      canu: \$( canu --version | cut -f 2 -d ' ' )
  END_VERSIONS
  """
}
