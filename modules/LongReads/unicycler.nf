process unicycler {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with Unicycler"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "unicycler/*" // Save all files
  tuple file("unicycler/unicycler_assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  lrID = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
