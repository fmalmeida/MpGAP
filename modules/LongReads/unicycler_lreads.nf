process unicycler_lreads_assembly {
  publishDir "${params.outdir}/${lrID}/${type}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with Unicycler"

  input:
  file lreads

  output:
  file "unicycler/*" // Save all files
  tuple file("unicycler/unicycler_assembly.fasta"), val(lrID), val('unicycler') // Gets contigs file

  script:
  lrID = (lreads - ".gz")[0].getBaseName()

  // Check available reads
  if (!params.shortreads_paired && !params.shortreads_single && params.longreads && params.lr_type) {
    type = 'longreads_only'
  } else if ((params.shortreads_paired || params.shortreads_single) && params.longreads && params.lr_type) {
    type = 'hybrid/strategy_2/longreads_only'
  }
  """
  unicycler -l $lreads \
  -o unicycler -t ${params.threads} \
  ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}
