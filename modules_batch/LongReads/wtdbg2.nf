process wtdbg2 {
  publishDir "${params.outdir}/${prefix}/wtdbg2", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "${id}: wtdgb2 assembly"

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), val(wtdbg2_technology), val(genomeSize), val(corrected_lreads), val(medaka_model),file(fast5), val(fast5_dir), file(bams), val(nBams), val(prefix)

  output:
  file "*" // Saves all files
  tuple val(id), file("wtdbg2_assembly.fasta"), val('wtdbg2') // Gets contigs file

  script:
  lr = (lr_type == 'nanopore') ? 'ont' : wtdbg2_technology

  """
  wtdbg2.pl -t ${params.threads} -x $lr -g ${genomeSize} -o ${id} ${params.wtdbg2_additional_parameters} $lreads

  # Rename contigs
  cp ${id}.cns.fa wtdbg2_assembly.fasta
  """
}
