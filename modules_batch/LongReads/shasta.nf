process shasta {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}: shasta assembly"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model),file(fast5), val(fast5_dir), file(bams), val(nBams), val(prefix)

  output:
  file "shasta/" // Saves all files
  tuple val(id), file("shasta/shasta_assembly.fasta"), val('shasta') // Gets contigs file

  when:
  (lr_type == 'nanopore')

  script:
  lr        = (lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  in_reads  = (lreads.getName() - ".gz")

  """
  # unzip reads
  gunzip -d -f $lreads

  # assemble
  shasta --assemblyDirectory shasta --threads ${params.threads} ${params.shasta_additional_parameters} --input $in_reads

  # Rename contigs
  cp shasta/Assembly.fasta shasta/shasta_assembly.fasta
  """
}
