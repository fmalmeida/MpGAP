// batch mode
process unicycler {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "unicycler" // Save everything
  tuple  val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler')

  when:
  (entrypoint == 'shortreads_only')

  script:
  param_paired = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "-1 $sread1 -2 $sread2" : ""
  param_single = !(single =~ /input.*/) ? "-s $single" : ""
  """
  # run unicycler
  unicycler $param_paired $param_single -o unicycler -t ${params.threads} ${params.unicycler_additional_parameters}

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}