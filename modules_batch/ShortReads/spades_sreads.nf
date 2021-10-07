// batch mode
process spades {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}: spades assembly"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix)

  output:
  file "spades" // Save all output
  tuple val(id), file("spades/spades_assembly.fasta"), val('spades')

  when:
  (entrypoint == 'sr-only')

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
