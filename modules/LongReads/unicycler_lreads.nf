process unicycler_lreads {
  publishDir "${params.outdir}/${lrID}/longreads_only", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with Unicycler"

  input:
  file lreads

  output:
  file "*" // Save all files
  tuple file("unicycler/unicycler_assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  lrID = lreads.getSimpleName()
  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters} &> unicycler.log

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
