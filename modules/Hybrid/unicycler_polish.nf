process pilon_polish {
  publishDir "${params.outdir}/${id}/hybrid/strategy_2/${out_ids}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Polishing a longreads-only assembly with shortreads (through Pilon)"

  input:
  tuple file(draft), val(id), val(assembler), file(reads)

  output:
  file("pilon_results_${assembler}/*") // Get everything
  tuple file("pilon_results_${assembler}/${assembler}_final_pilon_polish.fasta"), val(id), val("${assembler}_pilon_polished")

  script:
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    srId = (reads[3].getName() - ".gz").toString().substring(0, (reads[3].getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${reads[1]}_and_${srId}"
  } else if ((params.shortreads_single) && (!params.shortreads_paired)) {
    id = (reads[3].getName() - ".gz").toString().substring(0, (reads[3].getName() - ".gz").toString().lastIndexOf("."))
    out_ids = "${id}"
  } else if ((params.shortreads_paired) && (!params.shortreads_single)) {
    out_ids = "${reads[0]}"
  }

  if(params.shortreads_paired != '' && params.shortreads_single == '')
      """
      # Create the results dir
      mkdir pilon_results_${assembler};

      # Execute Unicycler polishing pilon wrapper
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon*/pilon*jar \
      -a $draft -1 ${reads[1]} -2 ${reads[2]} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_results_${assembler};
      mv pilon_results_${assembler}/*_final_polish.fasta pilon_results_${assembler}/${assembler}_final_pilon_polish.fasta ;
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
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon*/pilon*jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output ${assembler}_final_pilon_polish \
      --outdir pilon_results_${assembler} &> pilon.log

      # save bam file in the desired directory
      mv ${id}_${assembler}_aln.bam pilon.log pilon_results_${assembler};
      """
  else if(params.shortreads_paired != '' && params.shortreads_single != '')
      """
      # Create the results dir
      mkdir pilon_results_${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t ${params.threads} ${draft} ${reads[3]} > ${id}_${assembler}_aln.sam ;
      samtools view -bS ${id}_${assembler}_aln.sam | samtools sort > ${id}_${assembler}_aln.bam ;
      samtools index ${id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java -Xmx${params.pilon_memory_limit}G -jar /miniconda/share/pilon*/pilon*jar \
      --genome ${draft} --bam ${id}_${assembler}_aln.bam --output first_polish \
      --outdir . &> pilon.log

      # Execute Unicycler polishing pilon wrapper (for paired reads)
      unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon*/pilon*jar \
      -a first_polish.fasta -1 ${reads[1]} -2 ${reads[2]} --threads ${params.threads} &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log pilon_results_${assembler};
      mv pilon_results_${assembler}/*_final_polish.fasta pilon_results_${assembler}/${assembler}_final_pilon_polish.fasta ;
      """
}
