process medaka {
  publishDir "${params.outdir}/${lrID}/${type}/medaka_polished_contigs", mode: 'copy'
  container 'fmalmeida/mpgap'
  tag "Polishing assembly with Medaka"

  input:
  tuple file(draft), val(lrID), val(assembler), file(reads)

  output:
  file "${assembler}" // Save everything
  tuple file("${assembler}/${assembler}_medaka_consensus.fa"), val(lrID), val("${assembler}_medaka") // Save medaka contigs

  script:
  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  source activate MEDAKA ;
  medaka_consensus -i $reads -d $draft -o ${assembler} -t ${params.threads} -m ${params.medaka_sequencing_model} ;
  mv ${assembler}/consensus.fasta ${assembler}/${assembler}_medaka_consensus.fa
  """
}
