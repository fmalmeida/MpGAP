process quast {
  publishDir "${params.outdir}/${id}/${out_dir}/00_quality_assessment", mode: 'copy'
  label 'main'
  tag "Assessing ${assembler} assembly quality"

  input:
  tuple file(contigs), val(id), val(assembler), file(reads)

  output:
  file("${assembler}")
  val(id)
  val(out_dir)

  script:
  // Check available reads
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    srId = (reads[3].getName() - ".gz").toString().substring(0, (reads[3].getName() - ".gz").toString().lastIndexOf("."))
    prId = reads[0]
    out_ids = "${prId}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    srId = (reads[3].getName() - ".gz").toString().substring(0, (reads[3].getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${srId}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    prId = reads[0]
    out_ids = "${prId}"
  }

  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
    out_dir = "${type}"
  } else if ((params.shortreads_paired || params.shortreads_single) && !params.longreads ) {
    type = 'shortreads_only'
    out_dir = "${type}"
  } else {
    if (params.strategy_2) {
    type = 'hybrid/strategy_2'
    } else {
    type = 'hybrid/strategy_1'
    }
    out_dir = "${type}/${out_ids}"
  }

  // Alignment parameters
  if (params.shortreads_paired && !params.shortreads_single) {
    quast_parameter = "--pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else if (!params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]}"
  } else if (params.shortreads_paired && params.shortreads_single) {
    quast_parameter = "--single ${reads[3]} --pe1 ${reads[1]} --pe2 ${reads[2]}"
  } else {
    ltype           = (params.lr_type == 'nanopore') ? "ont2d" : "pacbio"
    quast_parameter = "--${params.lr_type} ${reads}"
  }

  """
  quast.py -o ${assembler} -t ${params.threads} ${quast_parameter} \\
  --circos --conserved-genes-finding --rna-finding ${contigs} --min-contig 100
  """
}
