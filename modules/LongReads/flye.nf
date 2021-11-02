process flye {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), val(shasta_config), file(bams), val(prefix)

  output:
  file "flye" // Saves all files
  tuple val(id), file("flye/flye_assembly.fasta"), val('flye') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '--nano' : '--pacbio'
  corrected = (corrected_lreads == 'true') ? '-corr' : '-raw'
  lrparam   = lr + corrected
  gsize     = (genomeSize) ? "--genome-size ${genomeSize}" : ""
  """
  source activate flye ;
  flye ${lrparam} $lreads ${gsize} --out-dir flye \
  --threads ${params.threads} ${params.flye_additional_parameters} &> flye.log ;
  mv flye/assembly.fasta flye/flye_assembly.fasta
  """
}
