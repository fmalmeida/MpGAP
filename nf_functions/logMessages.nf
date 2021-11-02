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
  summary['Number of threads     '] = params.threads
  // Workflow information
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  //summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  summary['Command used']   = "$workflow.commandLine"
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "========================================="
}
