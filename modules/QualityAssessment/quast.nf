process quast {
  publishDir "${params.output}", mode: 'copy', saveAs: { filename ->
    if ( filename.tokenize('/').contains('input_assembly') ) "final_assemblies/${asm_copy_prefix}_${filename.tokenize('/')[1]}"
    else "${prefix}/00_quality_assessment/${filename}"
  }
  tag "${id}"

  input:
  tuple val(id), file(contigs), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  tuple val(id), val(entrypoint), val(prefix), file("${assembler}"), emit: results
  file("input_assembly/*")

  script:

  // Alignment parameters
  paired_param = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "--pe1 ${sread1} --pe2 ${sread2}" : ""
  single_param = !(single =~ /input.?/) ? "--single ${single}" : ""
  lreads_param = !(lreads =~ /input.?/) ? "--${lr_type} ${lreads}" : ""
  additional_params = (params.quast_additional_parameters) ? params.quast_additional_parameters : ""

  // define prefix for final assemblies
  asm_copy_prefix = id.replaceAll(':', '_') // fixes hybrid prefixes that has a ':'

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
  
  # save assembly
  mkdir -p input_assembly
  cp ${contigs} input_assembly/${contigs}
  """

  else if (params.selected_profile == "singularity")
  """
  # fix busco usage in singularity
  mkdir -p ~/.quast/busco
  cp -R /opt/conda/envs/mpgap-*/lib/python3.8/site-packages/quast_libs/busco ~/.quast

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
  
  # save assembly
  mkdir -p input_assembly
  cp ${contigs} input_assembly/${contigs}
  """
}
