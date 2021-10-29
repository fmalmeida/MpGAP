process raven {
  publishDir "${params.outdir}/${prefix}/raven", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix)

  output:
  file "raven_assembly.*" // Saves all files
  tuple val(id), file("raven_assembly.fasta"), val('raven') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  corrected = (corrected_lreads == 'true') ? '--weaken' : ''
  """
  # Activate env
  source activate RAVEN;

  # Run
  raven --threads ${params.threads} --graphical-fragment-assembly raven_assembly.gfa \
  ${params.raven_additional_parameters} $corrected $lreads > raven_assembly.fasta ;
  """
}
