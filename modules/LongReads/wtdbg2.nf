process wtdbg2 {
  publishDir "${params.outdir}/${prefix}/wtdbg2", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with wtdbg2"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "*" // Saves all files
  tuple file("wtdbg2_assembly.fasta"), val(lrID), val('wtdbg2') // Gets contigs file

  script:
  lr                = (params.lr_type == 'nanopore') ? 'ont' : params.wtdbg2_technology
  lrID              = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  """
  wtdbg2.pl -t ${params.threads} -x $lr -g ${params.genomeSize} -o ${lrID} ${params.wtdbg2_additional_parameters} $lreads

  # Rename contigs
  cp ${lrID}.cns.fa wtdbg2_assembly.fasta
  """
}
