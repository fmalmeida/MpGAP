process unicycler_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "unicycler" // Save everything
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file
  path('versions.yml')

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/) && (entrypoint == 'hybrid_strategy_1')

  script:
  // Check reads
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "-1 $sread1 -2 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "-s $single" : ""
  additional_params = (params.unicycler_additional_parameters) ? params.unicycler_additional_parameters : ""
  """
  # run unicycler
  unicycler \\
      ${paired_reads} \\
      ${single_reads} \\
      -l ${lreads} \\
      -o unicycler \\
      -t $task.cpus \\
      $additional_params

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      spades: \$( spades.py --version | cut -f 4 -d ' ' )
      unicycler: \$( unicycler --version | cut -f 2 -d ' ' )
  END_VERSIONS
  """
}
