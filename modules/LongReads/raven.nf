process raven_assembly {
  publishDir "${params.outdir}/${lrID}/${type}/raven", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with raven"

  input:
  file lreads

  output:
  file "raven_assembly.*" // Saves all files
  tuple file("raven_assembly.fa"), val(lrID), val('raven') // Gets contigs file

  script:
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  corrected = (params.corrected_lreads) ? '--weaken' : ''

  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  # Activate env
  source activate RAVEN;

  # Run
  raven --threads ${params.threads} --graphical-fragment-assembly raven_assembly.gfa \
  ${params.raven_additional_parameters} $corrected $lreads > raven_assembly.fa ;
  """
}
