process shasta {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with shasta"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "shasta/" // Saves all files
  tuple file("shasta/shasta_assembly.fasta"), val(lrID), val('shasta') // Gets contigs file

  script:
  lr        = (params.lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  in_reads  = (lreads.getName() - ".gz")
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  """
  # unzip reads
  gunzip -dcf $lreads > uncompressed_${in_reads}

  # assemble
  shasta --assemblyDirectory shasta --threads ${params.threads} ${params.shasta_additional_parameters} --input uncompressed_${in_reads}

  # Rename contigs
  cp shasta/Assembly.fasta shasta/shasta_assembly.fasta
  """
}
