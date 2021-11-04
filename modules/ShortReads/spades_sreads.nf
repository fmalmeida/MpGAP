// batch mode
process spades {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "spades" // Save all output
  tuple val(id), file("spades/spades_assembly.fasta"), val('spades')

  when:
  (entrypoint == 'shortreads_only')

  script:
  param_paired = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "-1 $sread1 -2 $sread2" : ""
  param_single = !(single =~ /input.*/) ? "-s $single" : ""
  """
  # run spades
  spades.py -o spades -t ${params.threads} ${params.spades_additional_parameters} $param_paired $param_single

  # Rename assembly
  mv spades/contigs.fasta spades/spades_assembly.fasta
  """
}
