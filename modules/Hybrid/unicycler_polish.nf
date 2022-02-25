process pilon_polish {
  publishDir "${params.output}/${prefix}/pilon_polished_contigs", mode: 'copy'
  tag "${id}"

  input:
  tuple val(id), file(draft), val(assembler), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genome_size), val(corrected_long_reads), val(medaka_model), file(fast5), val(nanopolish_max_haplotypes), val(shasta_config), file(bams), val(prefix)

  output:
  file("${assembler}/*") // Get everything
  tuple val(id), file("${assembler}/${assembler}_pilon_consensus.fasta"), val("${assembler}_pilon_polished")

  script:
  has_paired = (sread1 =~ /input.*/) ? false : true
  has_single = (single =~ /input.*/) ? false : true
  fixed_id   = id - ":strategy_2"
  // has paired reads but doesn't have unpaired reads
  if(has_paired && !has_single)
      """
      # get tools path
      miniasm_path=\$(find \$CONDA_PREFIX -name "miniasm" | grep "mpgap" | head -n 1)
      pilonjar_path=\$(find \$CONDA_PREFIX -name "pilon.jar" | grep "mpgap" | head -n 1)

      # Create the results dir
      mkdir ${assembler};

      # Execute Unicycler polishing pilon wrapper
      unicycler_polish \\
          --minimap \$miniasm_path \\
          --pilon \$pilonjar_path \\
          -a $draft \\
          -1 ${sread1} \\
          -2 ${sread2} \\
          --threads $task.cpus &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log ${assembler};
      mv ${assembler}/*_final_polish.fasta ${assembler}/${assembler}_pilon_consensus.fasta ;
      """
  // doesn't have paired reads but has unpaired reads
  else if(!has_paired && has_single)
      """
      # get tools path
      miniasm_path=\$(find \$CONDA_PREFIX -name "miniasm" | grep "mpgap" | head -n 1)
      pilonjar_path=\$(find \$CONDA_PREFIX -name "pilon.jar" | grep "mpgap" | head -n 1)

      # Create the results dir
      mkdir ${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t $task.cpus ${draft} ${single} > ${fixed_id}_${assembler}_aln.sam ;
      samtools view -bS ${fixed_id}_${assembler}_aln.sam | samtools sort > ${fixed_id}_${assembler}_aln.bam ;
      samtools index ${fixed_id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java \\
          -Xmx${params.pilon_memory_limit}G \\
          -jar \$pilonjar_path \\
          --genome ${draft} \\
          --bam ${fixed_id}_${assembler}_aln.bam \\
          --output ${assembler}_pilon_consensus \\
          --outdir ${assembler} &> pilon.log

      # save bam file in the desired directory
      mv ${fixed_id}_${assembler}_aln.bam pilon.log ${assembler};
      """
  // has paired reads and has unpaired reads
  else if(has_paired && has_single)
      """
      # get tools path
      miniasm_path=\$(find \$CONDA_PREFIX -name "miniasm" | grep "mpgap" | head -n 1)
      pilonjar_path=\$(find \$CONDA_PREFIX -name "pilon.jar" | grep "mpgap" | head -n 1)

      # Create the results dir
      mkdir ${assembler};

      # Index and align reads with bwa
      bwa index ${draft} ;
      bwa mem -M -t $task.cpus ${draft} ${single} > ${fixed_id}_${assembler}_aln.sam ;
      samtools view -bS ${fixed_id}_${assembler}_aln.sam | samtools sort > ${fixed_id}_${assembler}_aln.bam ;
      samtools index ${fixed_id}_${assembler}_aln.bam ;

      # Execute pilon a single time (for single end reads)
      java \\
          -Xmx${params.pilon_memory_limit}G \\
          -jar \$pilonjar_path \\
          --genome ${draft} \\
          --bam ${fixed_id}_${assembler}_aln.bam \\
          --output first_polish \\
          --outdir . &> pilon.log

      # Execute Unicycler polishing pilon wrapper (for paired reads)
      unicycler_polish \\
          --minimap \$miniasm_path \\
          --pilon \$pilonjar_path \\
          -a first_polish.fasta \\
          -1 ${sread1} \\
          -2 ${sread2} \\
          --threads $task.cpus &> polish.log ;

      # Save files in the desired directory
      mv 0* polish.log ${assembler};
      mv ${assembler}/*_final_polish.fasta ${assembler}/${assembler}_pilon_consensus.fasta ;
      """
}
