process multiqc {
  publishDir "${params.output}/${prefix}/00_quality_assessment", mode: 'copy'
  label 'process_low'

  input:
  tuple val(id), val(entrypoint), val(prefix), file(quast_dirs)
  val nfRun
  path versions
  path config

  output:
  file "multiqc_report_${nfRun}.html"
  file "multiqc_data"
  file "ASSEMBLY_SUMMARY.txt"

  script:
  """
  # Run
  multiqc . \\
    --ignore "*.sam" \\
    --ignore "*.bam" \\
    --ignore "*.err" \\
    --ignore "*.stat" \\
    --config $config

  # Rename to have nf run name
  mv multiqc_report.html multiqc_report_${nfRun}.html ;

  # Create the markdown file resuming the main statistics
  echo \"# A summary of the main assembly statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main QUAST statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk cut -t -f 1,14,15,16,17,18,22,27,28,31 multiqc_data/multiqc_quast.txt | csvtk -t pretty >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  echo \"## Main BUSCO statistics\" >> ASSEMBLY_SUMMARY.txt
  echo \"\" >> ASSEMBLY_SUMMARY.txt
  csvtk -t pretty multiqc_data/multiqc_busco.txt >> ASSEMBLY_SUMMARY.txt
  """
}
