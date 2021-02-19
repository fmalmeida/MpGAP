process shovill_sreads_assembly {
  publishDir "${params.outdir}/${id}/shortreads_only", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Performing a illumina-only assembly with shovill, using paired end reads"
  cpus params.threads

  input:
  tuple val(id), file(sread1), file(sread2)

  output:
  file "shovill" // Save all output
  tuple file("shovill/shovill_contigs.fa"), val(id), val('shovill') // Gets contigs file

  when:
  ((params.shortreads_paired) && (!params.shortreads_single))

  script:
  """
  # Activate env
  source activate shovill ;

  # Run
  shovill --outdir shovill --R1 $sread1 --R2 $sread2 \
  --cpus ${params.threads} --trim ${params.shovill_additional_parameters}

  # Rename assembly
  mv shovill/contigs.fa shovill/shovill_contigs.fa
  """
}
