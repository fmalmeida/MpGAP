process pilon_polish {
  publishDir "${params.outdir}/polished-with-shortreads", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus params.threads

  input:
  tuple file(draft), val(id), val(assembler)
  file(reads)

  output:
  file("pilon_results_${assembler}") // Get everything
  tuple file("pilon_results_${assembler}/${assembler}_final_polish.fasta"), val(id), val("${assembler}-pilon-polished")

  script:
  """
  mkdir pilon_results_${assembler};
  unicycler_polish --minimap /miniconda/bin/miniasm --pilon /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
  -a $draft -1 ${reads[1]} -2 ${reads[2]} --threads ${params.threads} &> polish.log ;
  mv 0* polish.log pilon_results_${assembler};
  mv pilon_results_${assembler}/*_final_polish.fasta pilon_results_${assembler}/${assembler}_final_polish.fasta;
  """
}
