process pilon_polish {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}"

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file("pilon_polished_${assembler}/*") // Get everything
  tuple val(id), file("pilon_polished_${assembler}/${assembler}_final_pilon_polish.fasta"), val("${assembler}_pilon_polished")

  script:
  has_paired = (sread1 =~ /input.*/) ? false : true
  has_single = (single =~ /input.*/) ? false : true
  // has paired reads but doesn't have unpaired reads
  if(has_paired && !has_single)
      """
      # Create the results dir
      mkdir pilon_polished_${assembler};

      # Execute Unicycler polishing pilon wrapper
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon*/pilon*jar \
      -a $draft -1 ${sread1} -2 ${sread2} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_polished_${assembler};
      mv pilon_polished_${assembler}/*_final_polish.fasta pilon_polished_${assembler}/${assembler}_final_pilon_polish.fasta ;
      """
  // doesn't have paired reads but has unpaired reads
  else if(!has_paired && has_single)
      """
      # Create the results dir
      mkdir pilon_polished_${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t ${params.threads} ${draft} ${single} > ${id}_${assembler}_aln.sam ;
      samtools view -bS ${id}_${assembler}_aln.sam | samtools sort > ${id}_${assembler}_aln.bam ;
      samtools index ${id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon*/pilon*jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output ${assembler}_final_pilon_polish \
      --outdir pilon_polished_${assembler} &> pilon.log

      # save bam file in the desired directory
      mv ${id}_${assembler}_aln.bam pilon.log pilon_polished_${assembler};
      """
  // has paired reads and has unpaired reads
  else if(has_paired && has_single)
      """
      # Create the results dir
      mkdir pilon_polished_${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t ${params.threads} ${draft} ${single} > ${id}_${assembler}_aln.sam ;
      samtools view -bS ${id}_${assembler}_aln.sam | samtools sort > ${id}_${assembler}_aln.bam ;
      samtools index ${id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon*/pilon*jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output first_polish \
      --outdir . &> pilon.log

      # Execute Unicycler polishing pilon wrapper (for paired reads)
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon*/pilon*jar \
      -a first_polish.fasta -1 ${sread1} -2 ${sread2} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_polished_${assembler};
      mv pilon_polished_${assembler}/*_final_polish.fasta pilon_polished_${assembler}/${assembler}_final_pilon_polish.fasta ;
      """
}
