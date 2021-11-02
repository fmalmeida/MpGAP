process shasta {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), val(shasta_config), file(bams), val(prefix)

  output:
  file "shasta/" // Saves all files
  tuple val(id), file("shasta/shasta_assembly.fasta"), val('shasta') // Gets contigs file

  when:
  (lr_type == 'nanopore') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  lr        = (lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  in_reads  = (lreads.getName() - ".gz")

  """
  # unzip reads
  gunzip -dcf $lreads > uncompressed_${in_reads} ;

  # assemble
  shasta \
    --assemblyDirectory shasta \
    --threads ${params.threads} \
    ${params.shasta_additional_parameters} \
    --input uncompressed_${in_reads} \
    --config ${shasta_config};

  # Rename contigs
  cp shasta/Assembly.fasta shasta/shasta_assembly.fasta ;
  """
}
