process shasta {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "shasta/" // Saves all files
  tuple val(id), file("shasta/shasta_assembly.fasta"), val('shasta') // Gets contigs file

  when:
  (lr_type == 'nanopore') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  in_reads  = (lreads.getName() - ".gz")
  additional_params = (params.shasta_additional_parameters) ? params.shasta_additional_parameters : ""
  """
  # unzip reads
  gunzip -dcf $lreads > uncompressed_${in_reads} ;

  # assemble
  shasta \\
      --assemblyDirectory shasta \\
      --threads $task.cpus \\
      $additional_params \\
      --input uncompressed_${in_reads} \\
      --config ${shasta_config} ;

  # rename contigs
  cp shasta/Assembly.fasta shasta/shasta_assembly.fasta ;
  """
}
