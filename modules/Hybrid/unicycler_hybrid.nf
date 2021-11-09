process unicycler_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Save everything
  tuple val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler') // Gets contigs file

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/) && (entrypoint == 'hybrid_strategy_1')

  script:
  // Check reads
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "-1 $sread1 -2 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "-s $single" : ""
  """
  # run unicycler
  unicycler ${paired_reads} ${single_reads} -l ${lreads} -o unicycler -t ${params.threads} ${params.unicycler_additional_parameters}

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
