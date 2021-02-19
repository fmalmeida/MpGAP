process multiqc {
  publishDir "${params.outdir}/${id}/quality_assessment/${dir}", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Collecting Quast quality reports"

  input:
  file(quast_dirs)
  val(id)
  val(dir)

  output:
  file "multiqc_report.html"

  script:
  """
  multiqc */report.tsv
  """
}
