process gcpp {
  publishDir "${params.output}/${prefix}/gcpp_polished_contigs", mode: 'copy'
  label 'main'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file "${assembler}_pbvariants.gff" // Save gff
  tuple val(id), file("${assembler}_pbconsensus.fasta"), val("${assembler}_gcpp") // Save contigs

 when:
 !(bams =~ /input.*/) && (lr_type == 'pacbio') && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2')

  script:
  """
  # Activate env
  source activate pacbio;

  # count bams
  nBams=\$(ls *.bam | wc -l) ;

  # generate genome index
  pbmm2 index -j ${params.threads} ${draft} draft.mmi ;

  # Align bam
  pbmm2 align -j ${params.threads} --sort draft.mmi ${bams} final_pbaligned.bam ;

  # run polisher
  samtools index final_pbaligned.bam ;
  samtools faidx ${draft} ;
  gcpp -r ${draft} -o ${assembler}_pbconsensus.fasta,${assembler}_pbvariants.gff -j ${params.threads} final_pbaligned.bam ;
  """
}
