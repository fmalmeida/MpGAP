process gcpp {
  publishDir "${params.outdir}/${lrID}/${type}/gcpp_polished_contigs", mode: 'copy'
  label 'main'
  tag "Computing pacbio assembly consensus with gcpp"
  cpus params.threads

  input:
  tuple file(draft), val(lrID), val(assembler), file(bams), val(nBams)

  output:
  file "${assembler}_pbvariants.gff" // Save gff
  tuple file("${assembler}_pbconsensus.fasta"), val("${lrID}"), val("${assembler}_gcpp") // Save contigs

  script:
  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  # Activate env
  source activate pacbio;

  # generate genome index
  pbmm2 index -j ${params.threads} ${draft} draft.mmi

  # Single bam
  if [ $nBams -eq 1 ];
  then
    pbmm2 align -j ${params.threads} --sort draft.mmi ${bams.join(" ")} final_pbaligned.bam ;

  # Multiple bams
  elif [ $nBams -gt 1 ];
  then
    for file in ${bams.join(" ")} ; do pbmm2 align -j ${params.threads} --sort draft.mmi \$file \${file%%.bam}_pbaligned.bam ; done
    samtools merge --threads ${params.threads} pacbio_merged.bam *_pbaligned.bam ;
    samtools sort -@ ${params.threads} -o final_pbaligned.bam pacbio_merged.bam ;
  
  fi

  # run polisher
  samtools index final_pbaligned.bam ;
  samtools faidx ${draft} ;
  gcpp -r ${draft} -o ${assembler}_pbconsensus.fasta,${assembler}_pbvariants.gff -j ${params.threads} final_pbaligned.bam ;
  """
}
