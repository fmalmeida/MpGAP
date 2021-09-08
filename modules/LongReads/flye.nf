process flye_assembly {
  publishDir "${params.outdir}/${prefix}", mode: 'copy'
  label 'main'
  cpus params.threads
  tag "Performing a longreads only assembly with Flye"

  input:
  tuple file(lreads), val(prefix)

  output:
  file "flye" // Saves all files
  tuple file("flye/flye_assembly.fasta"), val(lrID), val('flye') // Gets contigs file

  script:
  lr        = (params.lr_type == 'nanopore') ? '--nano' : '--pacbio'
  corrected = (params.corrected_lreads) ? '-corr' : '-raw'
  lrparam   = lr + corrected
  lrID      = (lreads.getName() - ".gz").toString().substring(0, (lreads.getName() - ".gz").toString().lastIndexOf("."))
  gsize     = (params.genomeSize) ? "--genome-size ${params.genomeSize}" : ""

  """
  source activate flye ;
  flye ${lrparam} $lreads ${gsize} --out-dir flye \
  --threads ${params.threads} ${params.flye_additional_parameters} &> flye.log ;
  mv flye/assembly.fasta flye/flye_assembly.fasta
  """
}
