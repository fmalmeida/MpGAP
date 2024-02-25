process spades_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Save everything
  tuple val(id), file("spades/spades_assembly.fasta"), val('spades') // Gets contigs file
  path('versions.yml'), emit: versions

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/) && (entrypoint == 'hybrid_strategy_1')

  script:
  // Check reads
  lr   = (lr_type == 'nanopore') ? '--nanopore' : '--pacbio'
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "-1 $sread1 -2 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "-s $single" : ""
  additional_params = (params.spades_additional_parameters) ? params.spades_additional_parameters : ""
  """
  # run spades
  spades.py \\
      -o spades \\
      -t $task.cpus \\
      $additional_params \\
      ${paired_reads} \\
      ${single_reads} \\
      ${lr} ${lreads}

  # rename results
  mv spades/contigs.fasta spades/spades_assembly.fasta

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      spades: \$( spades.py --version | cut -f 4 -d ' ' )
  END_VERSIONS
  """
}
