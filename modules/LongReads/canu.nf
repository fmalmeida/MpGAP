process canu_assembly {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with Canu"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "canu/" // Saves all files
  tuple file("canu/canu_assembly.fasta"), val(lrID), val('canu') // Gets contigs file

  script:
  lr        = (params.lr_type == 'nanopore') ? '-nanopore' : '-pacbio'
  corrected = (params.corrected_lreads) ? '-corrected' : ''
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))

  """
  canu -p ${lrID} -d canu maxThreads=${params.threads} genomeSize=${params.genomeSize} \
  ${params.canu_additional_parameters} $corrected $lr $lreads

  # Rename contigs
  mv canu/${lrID}.contigs.fasta canu/canu_assembly.fasta
  """
}
