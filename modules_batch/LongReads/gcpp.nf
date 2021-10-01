process gcpp {
  publishDir "${params.outdir}/${prefix}/gcpp_polished_contigs", mode: 'copy'
  label 'main'
  tag "${id}: gcpp consensus"
  cpus params.threads

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model),file(fast5), val(fast5_dir), file(bams), val(nBams), val(prefix)

  output:
  file "${assembler}_pbvariants.gff" // Save gff
  tuple val(id), file("${assembler}_pbconsensus.fasta"), val("${assembler}_gcpp") // Save contigs

  when:
  !(bams =~ /input.*/)

  script:
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
