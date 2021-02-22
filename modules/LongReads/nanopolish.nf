process nanopolish {
  publishDir "${params.outdir}/${lrID}/${type}/nanopolished_contigs/${assembler}", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  tuple file(draft), val(lrID), val(assembler), file(reads), file(fast5), val(fast5_dir)

  output:
  tuple file("${assembler}_nanopolished.fa"), val(lrID), val("${assembler}_nanopolish") // Save nanopolished contigs
  file "${assembler}_nanopolished.complete.vcf" // Save VCF

  script:
  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2'
  }
  """
  source activate NANOPOLISH ;
  seqtk seq -A ${reads} > reads.fa ;
  nanopolish index -d "${fast5_dir}" reads.fa ;
  minimap2 -d draft.mmi ${draft} ;
  minimap2 -ax map-ont -t ${params.threads} ${draft} reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp ;
  samtools index reads.sorted.bam ;
  python /miniconda/envs/NANOPOLISH/bin/nanopolish_makerange.py ${draft} | parallel --results nanopolish.results -P ${params.cpus} \
  nanopolish variants --consensus -o polished.{1}.vcf \
    -w {1} \
    -r reads.fa \
    -b reads.sorted.bam \
    -g ${draft} \
    --max-haplotypes ${params.nanopolish_max_haplotypes} \
    --min-candidate-frequency 0.1;
  nanopolish vcf2fasta --skip-checks -g ${draft} polished.*.vcf > ${assembler}_nanopolished.fa ;
  cat polished.*.vcf >> ${assembler}_nanopolished.complete.vcf
  """
}
