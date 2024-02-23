// batch mode
process megahit {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_longreads), val(high_quality_longreads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "megahit" // Save everything
  tuple val(id), file("megahit/megahit_assembly.fasta"), val('megahit')
  path('versions.yml')

  when:
  (entrypoint == 'shortreads_only')

  script:
  param_paired = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "-1 $sread1 -2 $sread2" : ""
  param_single = !(single =~ /input.*/) ? "-r $single" : ""
  additional_params = (params.megahit_additional_parameters) ? params.megahit_additional_parameters : ""
  memory = "$task.memory" - " GB" + "e9"
  """
  # run megahit
  megahit \\
      ${param_paired} \\
      ${param_single} \\
      -o megahit \\
      -t $task.cpus \\
      -m $memory \\
      $additional_params

  # rename results
  mv megahit/final.contigs.fa megahit/megahit_assembly.fasta

  # get version
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      megahit: \$( megahit --version | cut -f 2 -d ' ' )
  END_VERSIONS
  """
}