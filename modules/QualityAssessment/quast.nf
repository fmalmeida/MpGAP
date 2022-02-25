process quast {
  publishDir "${params.output}/${prefix}/00_quality_assessment", mode: 'copy'
  tag "${id}"
  label 'process_low'

  input:
  tuple val(id), file(contigs), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  tuple val(id), val(entrypoint), val(prefix), file("${assembler}")

  script:
  // Alignment parameters
  paired_param = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "--pe1 ${sread1} --pe2 ${sread2}" : ""
  single_param = !(single =~ /input.?/) ? "--single ${single}" : ""
  lreads_param = !(lreads =~ /input.?/) ? "--${lr_type} ${lreads}" : ""
  additional_params = (params.quast_additional_parameters) ? params.quast_additional_parameters : ""

  if (params.selected_profile == "docker" || params.selected_profile == "conda")
  """
  # run quast
  quast.py \\
      -o ${assembler} \\
      -t $task.cpus \\
      ${lreads_param} \\
      ${paired_param} \\
      ${single_param} \\
      --conserved-genes-finding \\
      --rna-finding \\
      --min-contig 100 \\
      $additional_params \\
      ${contigs}
  """

  else if (params.selected_profile == "singularity")
  """
  # fix busco usage in singularity
  mkdir -p ~/.quast/busco
  cp -R /opt/conda/envs/mpgap-*/lib/python3.6/site-packages/quast_libs/busco ~/.quast

  # run quast
  quast.py \\
      -o ${assembler} \\
      -t $task.cpus \\
      ${lreads_param} \\
      ${paired_param} \\
      ${single_param} \\
      --conserved-genes-finding \\
      --rna-finding \\
      --min-contig 100 \\
      $additional_params \\
      ${contigs}
  """
}
