process variantCaller {
  publishDir "${params.outdir}/${lrID}/longreads_only/arrow_polished_contigs", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  tuple file(draft), val(lrID), val(assembler)
  file bams
  val nBams

  output:
  file "${assembler}_pbvariants.gff" // Save gff
  tuple file("${assembler}_pbconsensus.fasta"), val("${lrID}"), val("${assembler}_arrow") // Save contigs

  script:
  id = "${bams}" - ".bam"

  // Guide: https://sr-c.github.io/2018/05/29/polish-pacbio-assembly/
  """
  # Activate env
  source activate pacbio;

  # Single bam
  if [ $nBams -eq 1 ];
  then
    pbalign --nproc ${params.threads} ${bams} ${draft} ${id}_pbaligned.bam
    samtools index ${id}_pbaligned.bam;
    pbindex ${id}_pbaligned.bam;
    samtools faidx ${draft};
    arrow -j ${params.threads} --referenceFilename ${draft} -o ${assembler}_pbconsensus.fasta \
    -o ${assembler}_pbvariants.gff ${id}_pbaligned.bam

  # Multiple bams
  elif [ $nBams -gt 1 ];
  then
    for BAM in ${bams.join(" ")} ; do pbalign --nproc ${params.threads}  \
    \$BAM ${draft} \${BAM%%.bam}_pbaligned.bam; done;
    for BAM in *_pbaligned.bam ; do samtools sort -@ ${params.threads} \
    -o \${BAM%%.bam}_sorted.bam \$BAM; done;
    samtools merge pacbio_merged.bam *_sorted.bam;
    samtools index pacbio_merged.bam;
    pbindex pacbio_merged.bam;
    samtools faidx ${draft};
    arrow -j ${params.threads} --referenceFilename ${draft} -o ${assembler}_pbconsensus.fasta \
    -o ${assembler}_pbvariants.gff pacbio_merged.bam
  fi
  """
}
