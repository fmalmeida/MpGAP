/*
 * Define log message
 */
def logMessage() {
  log.info "===================================================================="
  log.info " Container-based, fmalmeida/mpgap, generic genome assembly pipeline "
  log.info "===================================================================="
  def summary = [:]
  // Generic parameters
  summary['Output directory      '] = params.outdir
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
  summary['Assembly method       '] = 'longreads-only'
  } else if ((params.shortreads_paired || params.shortreads_single) && !params.longreads ) {
  summary['Assembly method       '] = 'shortreads-only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads ) {
  summary['Assembly method       '] = 'hybrid'
  }
  summary['Number of threads     '] = params.threads
  // Long reads?
  if (params.longreads) {
  summary['Longreads             '] = params.longreads
  }
  if (params.longreads) {
  summary['Longread technology   '] = params.lr_type
  }
  if (params.nanopolish_fast5 && params.lr_type =='nanopore') {
  summary['Fast5 files dir       '] = params.nanopolish_fast5
  }
  if (params.medaka_sequencing_model && params.lr_type =='nanopore') {
  summary['Medaka model          '] = params.medaka_sequencing_model
  }
  if (params.pacbio_bam && params.lr_type =='pacbio') {
  summary['Pacbio subreads BAM   '] = params.pacbio_bam
  }
  // Short reads?
  if (params.shortreads_single) {
  summary['Short single end reads'] = params.shortreads_single
  }
  if (params.shortreads_paired) {
  summary['Short paired end reads'] = params.shortreads_paired
  }
  // Workflow information
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  //summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  summary['Command used']   = "$workflow.commandLine"
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "========================================="
}
