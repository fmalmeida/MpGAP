process variantCaller {
  publishDir "${params.outdir}/longreads-only/arrowPolished_contigs", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  tuple file(draft), val(lrID), val(assembler)
  file bams
  val nBams

  output:
  file "${params.prefix}_${assembler}_pbvariants.gff" // Save gff
  file "${params.prefix}_${assembler}_pbconsensus.fasta" // Save contigs

  script:
  id = "${bams}" - ".bam"
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
    arrow -j ${params.threads} --referenceFilename ${draft} -o ${params.prefix}_${assembler}_pbconsensus.fasta \
    -o ${params.prefix}_${assembler}_pbvariants.gff ${id}_pbaligned.bam

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
    arrow -j ${params.threads} --referenceFilename ${draft} -o ${params.prefix}_${assembler}_pbconsensus.fasta \
    -o ${params.prefix}_${assembler}_pbvariants.gff pacbio_merged.bam
  fi
  """
}
