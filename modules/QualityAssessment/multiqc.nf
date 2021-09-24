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
  file "multiqc_data"
  file "ASSEMBLY_SUMMARY.txt"

  script:
  """
  # Run
  multiqc */report.tsv */busco_stats/short_summary_* ;

  # Rename to have nf run name
  mv multiqc_report.html multiqc_report_${nfRun}.html ;

  # Create the markdown file resuming the main statistics
  echo \"# A summary of the main assembly statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main QUAST statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk cut -t -f 1,14,15,16,17,18,22,25,26,32 multiqc_data/multiqc_quast.txt | csvtk -t pretty >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main BUSCO statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk -t pretty multiqc_data/multiqc_busco.txt >> ASSEMBLY_SUMMARY.txt
  """
}


// batch mode
process multiqc_batch {
  publishDir "${params.outdir}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Collecting Quast quality reports"

  input:
  tuple file(quast_dirs), val(prefix)
  val nfRun

  output:
  file "multiqc_report_${nfRun}.html"
  file "multiqc_data"
  file "ASSEMBLY_SUMMARY.txt"

  script:
  """
  # Run
  multiqc */report.tsv */busco_stats/short_summary_* ;

  # Rename to have nf run name
  mv multiqc_report.html multiqc_report_${nfRun}.html ;

  # Create the markdown file resuming the main statistics
  echo \"# A summary of the main assembly statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main QUAST statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk cut -t -f 1,14,15,16,17,18,22,25,26,32 multiqc_data/multiqc_quast.txt | csvtk -t pretty >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main BUSCO statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk -t pretty multiqc_data/multiqc_busco.txt >> ASSEMBLY_SUMMARY.txt
  """
}
