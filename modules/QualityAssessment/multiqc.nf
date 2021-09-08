process multiqc {
  publishDir "${params.outdir}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Collecting Quast quality reports"

  input:
  file(quast_dirs)
  val prefix
  val nfRun

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
