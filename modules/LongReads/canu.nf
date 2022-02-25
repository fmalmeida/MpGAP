process canu {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "canu/" // Saves all files
  tuple val(id), file("canu/canu_assembly.fasta"), val('canu') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  corrected = (corrected_long_reads == 'true') ? '-corrected' : ''
  fixed_id = id - ":strategy_2"
  """
  # run canu
  canu \\
      -p ${fixed_id} \\
      -d canu \\
      maxThreads=$task.cpus \\
      genomeSize=${genome_size} \\
      ${params.canu_additional_parameters} \\
      $corrected \\
      $lr $lreads

  # rename results
  mv canu/${fixed_id}.contigs.fasta canu/canu_assembly.fasta
  """
}
