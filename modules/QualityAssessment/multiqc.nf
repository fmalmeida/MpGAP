process multiqc {
  publishDir "${params.outdir}/${id}/${dir}/00_quality_assessment", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Collecting Quast quality reports"

  input:
  file(quast_dirs)
  val(id)
  val(dir)
  val(nfRun)

  output:
  file "multiqc_report_${nfRun}.html"
  file "multiqc_data" optional true

  script:
  """
  # Run
  multiqc */report.tsv */busco_stats/short_summary_* ;

  # Rename to have nf run name
  mv multiqc_report.html multiqc_report_${nfRun}.html ;
  """
}
