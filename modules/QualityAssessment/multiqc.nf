process multiqc {
  publishDir "${params.outdir}/${id}/quality_assessment/${dir}", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Collecting Quast quality reports"

  input:
  file(quast_dirs)
  val(id)
  val(dir)
  val(nfRun)

  output:
  file "multiqc_report_${nfRun}.html"

  script:
  """
  # Run
  multiqc */report.tsv

  # Rename to have nf run name
  mv multiqc_report.html multiqc_report_${nfRun}.html
  """
}
