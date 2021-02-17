/*
 * Define log message
 */
def logMessage() {
  log.info "================================================================="
  log.info " Docker-based, fmalmeida/mpgap, generic genome assembly pipeline "
  log.info "================================================================="
  def summary = [:]
  // Generic parameters
  summary['Output directory']    = params.outdir
  summary['Assembly method']     = params.assembly_type
  summary['Number of threads']   = params.threads
  // Long reads?
  if (params.longreads) { summary['Longreads']   = params.longreads }
  if (params.longreads) { summary['Longread technology'] = params.lr_type }
  if (params.nanopolish_fast5Path && params.lr_type == 'nanopore') { summary['Fast5 files dir']   = params.nanopolish_fast5Path }
  if (params.medaka_sequencing_model && params.lr_type == 'nanopore') { summary['Medaka model']   = params.medaka_sequencing_model }
  if (params.pacbio_all_bam_path && params.lr_type == 'pacbio') { summary['Pacbio subreads BAM']   = params.pacbio_all_bam_path }
  // Short reads?
  if (params.shortreads_single) { summary['Short single end reads']   = params.shortreads_single }
  if (params.shortreads_paired) { summary['Short paired end reads']   = params.shortreads_paired }
  // Workflow information
  if(workflow.revision) summary['Pipeline Release'] = workflow.revision
  summary['Current home']   = "$HOME"
  summary['Current user']   = "$USER"
  summary['Current path']   = "$PWD"
  summary['Command used']   = "$workflow.commandLine"
  log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
  log.info "========================================="
}
