process raven {
  publishDir "${params.output}/${prefix}/raven", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "raven_assembly.*" // Saves all files
  tuple val(id), file("raven_assembly.fasta"), val('raven') // Gets contigs file
  path('versions.yml')

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  additional_params = (params.raven_additional_parameters) ? params.raven_additional_parameters : ""
  corrected = ''
  if (corrected_longreads.toBoolean()) {
    if ( additional_params.tokenize(' ').intersect( ['-k', '--kmer-len'  ] ) == 0 ) { corrected = corrected + '-k 30'}
    if ( additional_params.tokenize(' ').intersect( ['-w', '--window-len'] ) == 0 ) { corrected = corrected + '-w 10'}
  }
  if (high_quality_longreads.toBoolean()) {
    if ( additional_params.tokenize(' ').intersect( ['-k', '--kmer-len'  ] ) == 0 ) { corrected = corrected + '-k 45'}
    if ( additional_params.tokenize(' ').intersect( ['-w', '--window-len'] ) == 0 ) { corrected = corrected + '-w 15'}
  }
  """
  # run raven
  raven \\
      --threads $task.cpus \\
      --graphical-fragment-assembly raven_assembly.gfa \\
      $additional_params \\
      $corrected \\
      $lreads > raven_assembly.fasta ;
  
  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      raven: \$( raven --version )
  END_VERSIONS
  """
}
