process flye {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "flye" // Saves all files
  tuple val(id), file("flye/flye_assembly.fasta"), val('flye') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '--nano' : '--pacbio'
  if (corrected_longreads.toBoolean())    { lrparam = lr + '-corr' }
  else if (high_quality_longreads.toBoolean()) { lrparam = lr + (lr_type == 'nanopore') ? '-hq' : '-hifi' }
  else { lrparam = lr + '-raw' }
  gsize     = (genome_size) ? "--genome-size ${genome_size}" : ""
  additional_params = (params.flye_additional_parameters) ? params.flye_additional_parameters : ""
  """
  # run flye
  flye \\
      ${lrparam} $lreads \\
      ${gsize} \\
      --out-dir flye \\
      $additional_params \\
      --threads $task.cpus &> flye.log ;

  # rename results
  mv flye/assembly.fasta flye/flye_assembly.fasta
  """
}
