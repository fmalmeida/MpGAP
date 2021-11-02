process spades_hybrid {
  publishDir "${params.output}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model), file(fast5), file(bams), val(prefix)

  output:
  file "*" // Save everything
  tuple val(id), file("spades/spades_assembly.fasta"), val('spades') // Gets contigs file

  when:
  ((!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) || !(single =~ /input.*/)) && !(lreads =~ /input.*/)

  script:
  // Check reads
  lr   = (lr_type == 'nanopore') ? '--nanopore' : '--pacbio'
  paired_reads = (!(sread1 =~ /input.*/) && !(sread2 =~ /input.*/)) ? "-1 $sread1 -2 $sread2" : ""
  single_reads = !(single =~ /input.*/) ? "-s $single" : ""
  """
  spades.py -o spades -t ${params.threads} \\
  ${params.spades_additional_parameters} ${paired_reads} ${single_reads} ${lr} ${lreads}

  # Rename
  mv spades/contigs.fasta spades/spades_assembly.fasta
  """
}
