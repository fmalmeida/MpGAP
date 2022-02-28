// batch mode
process unicycler {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

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
  additional_params = (params.unicycler_additional_parameters) ? params.unicycler_additional_parameters : ""
  """
  # copy spades 3.13 to dir
  src_dir=\$(which shasta | sed 's/shasta//g')
  spades_path="\${src_dir}/spades-3.13.tar.gz"
  cp \${spades_path} .
  tar zxvf spades-3.13.tar.gz && rm spades-3.13.tar.gz

  # run unicycler
  unicycler \\
      ${param_paired} \\
      ${param_single} \\
      -o unicycler \\
      -t $task.cpus \\
      $additional_params \\
      --spades_path SPAdes-3.13.0-Linux/bin/spades.py

  # rename results
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}