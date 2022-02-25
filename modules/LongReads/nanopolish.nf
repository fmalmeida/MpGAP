process nanopolish {
  publishDir "${params.output}/${prefix}/nanopolished_contigs/${assembler}", mode: 'copy'
  tag "${id}"

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  tuple val(id), file("${assembler}_nanopolish_consensus.fa"), val("${assembler}_nanopolish") // Save nanopolished contigs
  file "${assembler}_nanopolish_consensus.complete.vcf" // Save VCF

  when:
  !(fast5 =~ /input.*/) && (lr_type == 'nanopore') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  fast5_dir = fast5.toString()
  """
  # save reads as fasta
  seqtk seq -A ${lreads} > reads.fa ;

  # filter contigs for polish
  # nanopolish accepts sequences greater than 40bp
  # will be used only to create window ranges
  seqtk seq -L 40 ${draft} > filtered_assembly.fa ;

  # index fast5 files
  nanopolish index -d "${fast5_dir}" reads.fa ;

  # index assembly
  minimap2 \\
      -d draft.mmi \\
      ${draft} ;
  
  # map reads to assembly
  minimap2 \\
      -ax map-ont \\
      -t $task.cpus \\
      ${draft} \\
      reads.fa | \\
      samtools \\
          sort \\
          -o reads.sorted.bam \\
          -T reads.tmp ;
  
  # index bam
  samtools index reads.sorted.bam ;

  # run nanopolish
  nanopolish_makerange.py \\
      filtered_assembly.fa | \\
      parallel --results nanopolish.results -P 1 \\
      nanopolish variants --consensus -o polished.{1}.vcf \\
          -w {1} \\
          -r reads.fa \\
          -b reads.sorted.bam \\
          -g ${draft} \\
          -t $task.cpus \\
          --max-haplotypes ${nanopolish_max_haplotypes} ;
  
  # call polished fasta from vcf
  nanopolish \\
      vcf2fasta \\
      --skip-checks \\
      -g ${draft} \\
      polished.*.vcf > ${assembler}_nanopolish_consensus.fa ;
  
  # rename contigs
  cat polished.*.vcf >> ${assembler}_nanopolish_consensus.complete.vcf
  """
}
