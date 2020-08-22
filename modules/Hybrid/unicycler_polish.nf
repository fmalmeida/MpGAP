process pilon_polish {
  publishDir "${params.outdir}/polished-with-shortreads", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads
  tag "Polishing a longreads-only assembly with shortreads (through Pilon)"

  input:
  tuple file(draft), val(id), val(assembler)
  file(reads)

  output:
  file("pilon_results_${assembler}") // Get everything
  tuple file("pilon_results_${assembler}/${assembler}_final_polish.fasta"), val(id), val("${assembler}-pilon-polished")

  script:
  if(params.shortreads_paired != '' && params.shortreads_single == '')
      """
      # Create the results dir
      mkdir pilon_results_${assembler};

      # Execute Unicycler polishing pilon wrapper
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
      -a $draft -1 ${reads[1]} -2 ${reads[2]} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_results_${assembler};
      mv pilon_results_${assembler}/*_final_polish.fasta pilon_results_${assembler}/${assembler}_final_polish.fasta ;
      """

  else if(params.shortreads_paired == '' && params.shortreads_single != '')
      """
      # Create the results dir
      mkdir pilon_results_${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t ${params.threads} ${draft} ${reads[3]} > ${id}_${assembler}_aln.sam ;
      samtools view -bS ${id}_${assembler}_aln.sam | samtools sort > ${id}_${assembler}_aln.bam ;
      samtools index ${id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output ${assembler}_final_polish \
      --outdir pilon_results_${assembler} &> pilon.log

      # save bam file in the desired directory
      mv ${id}_${assembler}_aln.bam pilon.log pilon_results_${assembler};
      """
  else if((params.shortreads_paired != '' && params.shortreads_single != ''))
      """
      # Create the results dir
      mkdir pilon_results_${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t ${params.threads} ${draft} ${reads[3]} > ${id}_${assembler}_aln.sam ;
      samtools view -bS ${id}_${assembler}_aln.sam | samtools sort > ${id}_${assembler}_aln.bam ;
      samtools index ${id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output first_polish \
      --outdir . &> pilon.log

      # Execute Unicycler polishing pilon wrapper (for paired reads)
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
      -a first_polish.fasta -1 ${reads[1]} -2 ${reads[2]} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_results_${assembler};
      mv pilon_results_${assembler}/*_final_polish.fasta pilon_results_${assembler}/${assembler}_final_polish.fasta ;
      """
}
