process raven {
  publishDir "${params.output}/${prefix}/raven", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "raven_assembly.*" // Saves all files
  tuple val(id), file("raven_assembly.fasta"), val('raven') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  corrected = (corrected_long_reads == 'true') ? '-k 30 -w 10' : ''
  additional_params = (params.raven_additional_parameters) ? params.raven_additional_parameters : ""
  """
  # run raven
  raven \\
      --threads $task.cpus \\
      --graphical-fragment-assembly raven_assembly.gfa \\
      $additional_params \\
      $corrected \\
      $lreads > raven_assembly.fasta ;
  """
}
