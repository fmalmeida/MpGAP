process unicycler_lreads {
  publishDir "${params.outdir}/${lrID}/longreads_only", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Performing a longreads only assembly with Unicycler"

  input:
  file lreads

  output:
  file "unicycler/" // Save all files
  tuple file("unicycler/assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  lrID = lreads.getSimpleName()
  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters} &> unicycler.log
  """
}
