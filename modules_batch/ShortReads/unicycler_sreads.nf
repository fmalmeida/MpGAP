// batch mode
process unicycler {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  tag "${id}: unicycler assembly"
  cpus params.threads

  input:
  tuple val(id), val(entrypoint), file(sread1), file(sread2), file(single), file(lreads), val(lr_type), file(fast5), val(fast5_dir), file(bams), val(nBams), val(prefix)

  output:
  file "unicycler" // Save everything
  tuple  val(id), file("unicycler/unicycler_assembly.fasta"), val('unicycler')

  script:
  param_paired = !(sread1 =~ /input.*/ || sread2 =~ /input.*/) ? "-1 $sread1 -2 $sread2" : ""
  param_single = !(single =~ /input.*/) ? "-s $single" : ""
  """
  # run unicycler
  unicycler $param_paired $param_single -o unicycler -t ${params.threads} ${params.unicycler_additional_parameters}

  # Rename assembly
  mv unicycler/assembly.fasta unicycler/unicycler_assembly.fasta
  """
}