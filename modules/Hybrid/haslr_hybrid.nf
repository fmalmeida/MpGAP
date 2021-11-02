process haslr_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), val(shasta_config), file(bams), val(prefix)

  output:
  file "*" // Save everything
  tuple val(id), file("haslr/haslr_assembly.fa"), val('haslr') // Gets contigs file

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/)

  script:
  // Check reads
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "$sread1 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "$single" : ""
  """
  haslr.py -t ${params.threads} -o haslr -g ${genomeSize} \
  -l $lreads -x ${lr_type} -s ${paired_reads} ${single_reads} \
  ${params.haslr_additional_parameters}

  # Rename
  cp haslr/*/asm.final.fa haslr/haslr_assembly.fa
  """
}
