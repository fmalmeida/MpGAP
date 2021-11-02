process wtdbg2 {
  publishDir "${params.output}/${prefix}/wtdbg2", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix)

  output:
  file "*" // Saves all files
  tuple val(id), file("wtdbg2_assembly.fasta"), val('wtdbg2') // Gets contigs file

  when:
  (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr = (lr_type == 'nanopore') ? 'ont' : wtdbg2_technology
  """
  wtdbg2.pl -t ${params.threads} -x $lr -g ${genomeSize} -o ${id} ${params.wtdbg2_additional_parameters} $lreads

  # Rename contigs
  cp ${id}.cns.fa wtdbg2_assembly.fasta
  """
}
