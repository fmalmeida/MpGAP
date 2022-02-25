process gcpp {
  publishDir "${params.output}/${prefix}/gcpp_polished_contigs", mode: 'copy'
  tag "${id}"
  label 'process_assembly'

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "${assembler}_gcpp_variants.gff" // Save gff
  tuple val(id), file("${assembler}_gcpp_consensus.fasta"), val("${assembler}_gcpp") // Save contigs

 when:
 !(bams =~ /input.*/) && (lr_type == 'pacbio') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  """
  # generate genome index
  pbmm2 \\
      index \\
      -j $task.cpus \\
      ${draft} \\
      draft.mmi ;

  # align bam
  pbmm2 \\
      align \\
      -j $task.cpus \\
      --sort \\
      draft.mmi \\
      ${bams} \\
      final_pbaligned.bam ;

  # index bam and fasta
  samtools index final_pbaligned.bam ;
  samtools faidx ${draft} ;

  # run polisher
  gcpp \\
      -r ${draft} \\
      -o ${assembler}_gcpp_consensus.fasta,${assembler}_gcpp_variants.gff \\
      -j $task.cpus \\
      final_pbaligned.bam ;
  """
}
