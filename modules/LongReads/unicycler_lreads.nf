process unicycler_lreads {
  publishDir "${params.outdir}/longreads-only", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  file lreads

  output:
  file "unicycler_${lrID}/" // Save all files
  tuple file("unicycler_${lrID}/assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  lrID = lreads.getSimpleName()
  """
  unicycler -l $lreads \
  -o unicycler_${lrID} -t ${params.threads} \
  ${params.unicycler_additional_parameters} &> unicycler.log
  """
}
