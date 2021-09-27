process raven {
  publishDir "${params.outdir}/${prefix}/raven", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with raven"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "raven_assembly.*" // Saves all files
  tuple file("raven_assembly.fasta"), val(lrID), val('raven') // Gets contigs file

  script:
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  corrected = (params.corrected_lreads) ? '--weaken' : ''

  """
  # Activate env
  source activate RAVEN;

  # Run
  raven --threads ${params.threads} --graphical-fragment-assembly raven_assembly.gfa \
  ${params.raven_additional_parameters} $corrected $lreads > raven_assembly.fasta ;
  """
}
