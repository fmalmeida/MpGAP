process medaka {
  publishDir "${params.outdir}/longreads-only/medaka_polished_contigs", mode: 'copy', overwrite: true
  container 'fmalmeida/mpgap'
  tag "Polishing assembly with Medaka"

  input:
  tuple file(draft), val(lrID), val(assembler)
  file reads

  output:
  file "${assembler}_${lrID}" // Save everything
  tuple file("${assembler}_${lrID}/${assembler}_${lrID}_consensus.fa"), val(lrID), val("${assembler}-medaka") // Save medaka contigs

  script:
  """
  source activate MEDAKA ;
  medaka_consensus -i $reads -d $draft -o ${assembler}_${lrID} -t ${params.threads} -m ${params.medaka_sequencing_model} ;
  mv ${assembler}_${lrID}/consensus.fasta ${assembler}_${lrID}/${assembler}_${lrID}_consensus.fa
  """
}
