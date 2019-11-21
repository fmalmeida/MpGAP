#!/usr/bin/env nextflow

/*
                      A docker-based pipeline for generic hybrid, illumina-only
                      or long reads only assembly. It accepts Illumina, ONT and
                      Pacbio data.

                      It uses Unicycler, Flye, Canu or Spades to assemble reads.
                      And uses Nanopolish, VariantCaller or Pilon to polish assemblies.

*/

def helpMessage() {
   log.info """
   Usage:
   nextflow run fmalmeida/MpGAP [--help] [ -c nextflow.config ] [OPTIONS] [-with-report] [-with-trace] [-with-timeline]

   Comments:
   This pipeline contains a massive amount of configuration variables and its usage as CLI parameters would
   cause the command to be huge. Therefore, it is extremely recommended to use the nextflow.config configuration file in order to make
   parameterization easier and more readable.

   Creating a configuration file:
   nextflow run fmalmeida/MpGAP [--get_hybrid_config] [--get_lreads_config] [--get_sreads_config] [--get_yaml]

   Show command line examples:
   nextflow run fmalmeida/MpGAP --show

   Execution Reports:
   nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-report
   nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-trace
   nextflow run fmalmeida/MpGAP [ -c nextflow.config ] -with-timeline

   OBS: These reports can also be enabled through the configuration file.

   OPTIONS:
            General Parameters - Mandatory

    --outDir <string>                      Output directory name
    --prefix <string>                      Set prefix for output files
    --threads <int>                        Number of threads to use
    --yaml <string>                        Sets path to yaml file containing additional parameters to assemblers.
    --assembly_type <string>               Selects assembly mode: hybrid, illumina-only or longreads-only
    --try_canu                             Execute assembly with Canu. Multiple assemblers can be chosen.
    --try_unicycler                        Execute assembly with Unicycler. Multiple assemblers can be chosen.
    --try_flye                             Execute assembly with Flye. Multiple assemblers can be chosen.
    --try_spades                           Execute assembly with Spades. Multiple assemblers can be chosen.


            Parameters for illumina-only mode. Can be executed by SPAdes and Unicycler assemblers.
            Users can use paired or single end reads. If both types are given at once, assemblers
            will be executed with a mix of both.

    --shortreads_paired <string>           Path to Illumina paired end reads.
    --shortreads_single <string>           Path to Illumina single end reads.
    --ref_genome <string>                  Path to reference genome for guided assembly. Used only by SPAdes.

            Parameters for hybrid mode. Can be executed by spades and unicycler assemblers.

    --shortreads_paired <string>           Path to Illumina paired end reads.
    --shortreads_single <string>           Path to Illumina single end reads.
    --ref_genome <string>                  Path to reference genome for guided assembly. Used only by SPAdes.
    --longreads <string>                   Path to longreads in FASTA or FASTQ formats.
    --lr_type <string>                     Sets wich type of long reads are being used: pacbio or nanopore

            Parameters for longreads-only mode. Can be executed by canu, flye and unicycler assemblers.
            In the end, long reads only assemblies can be polished with illumina reads through pilon.

    --longreads <string>                   Path to longreads in FASTA or FASTQ formats.
    --fast5Path <string>                   Path to directory containing FAST5 files for given reads.
                                           Whenever set, the pipeline will execute a polishing step
                                           with Nanopolish. This makes the pipeline extremely SLOW!!
    --pacbio_all_baxh5_path <string>       Path to all bax.h5 files for given reads. Whenever set, the pipeline
                                           will execute a polishing step with VarianCaller.
    --pacbio_all_bam_path <string>         Path to all subreads bam files for given reads. Whenever set, the pipeline
                                           will execute a polishing step with VarianCaller.
    --genomeSize                           Canu and Flye require an estimative of genome size in order
                                           to be executed. Examples: 5.6m; 1.2g
    --lr_type <string>                     Sets wich type of long reads are being used: pacbio or nanopore
    --illumina_polish_longreads_contigs    This tells the pipeline to polish long reads only assemblies
                                           with Illumina reads through Pilon. This is another hybrid methodology.
                                           For that, users have to set path to Illumina reads through
                                           --shortreads_paired or --shortreads_single.



   """.stripIndent()
}

def exampleMessage() {
   log.info """
   """.stripIndent()
}

/*
          Display Help Message
*/
params.help = false
 // Show help emssage
 if (params.help){
   helpMessage()
   //file('work').deleteDir()
   exit 0
}

/*
          Display CLI examples
*/
params.show = false
 // Show help emssage
 if (params.show){
   exampleMessage()
   exit 0
}

/*
          Download configuration file, if necessary.
*/
params.get_hybrid_config = false
if (params.get_hybrid_config) {
  new File("hybrid.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/hybrid.config").getText()
  println ""
  println "hybrid.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/NGS-preprocess -c ./hybrid.config"
  println "Nice code!\n"

  exit 0
}
params.get_lreads_config = false
if (params.get_lreads_config) {
  new File("lreads-only.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/lreads.config").getText()
  println ""
  println "lreads.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/NGS-preprocess -c ./lreads.config"
  println "Nice code!\n"

  exit 0
}
params.get_sreads_config = false
if (params.get_sreads_config) {
  new File("sreads-only.config") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/configuration_example/sreads.config").getText()
  println ""
  println "sreads.config file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/NGS-preprocess -c ./sreads.config"
  println "Nice code!\n"

  exit 0
}

params.get_yaml = false
if ( params.get_yaml ) {
  new File("additional_parameters.yaml") << new URL ("https://github.com/fmalmeida/MpGAP/raw/master/additional_parameters.yaml").getText()
  println ""
  println "additional_parameters.yaml file saved in working directory"
  println "After configuration, run:"
  println "nextflow run fmalmeida/NGS-preprocess -c ./*.config"
  println "Nice code!\n"

  exit 0
}

/*
                    Setting Default Parameters.
                    Do not change any of this values
                    directly in main.nf.
                    Use the config file instead.

*/

params.longreads = ''
params.fast5Path = ''
params.pacbio_all_baxh5_path = ''
params.pacbio_all_bam_path = ''
params.lr_type = ''
params.shortreads_paired = ''
params.shortreads_single = ''
params.ref_genome = ''
params.assembly_type = ''
params.illumina_polish_longreads_contigs = false
params.pilon_memmory_limit = 50
params.try_canu = false
params.try_unicycler = false
params.try_flye = false
params.try_spades = false
params.genomeSize = ''
params.outDir = 'output'
params.prefix = 'out'
params.threads = 3
params.cpus = 2
params.yaml = ""

/*
                    Loading Parameters properly
                    set through config file.

*/

prefix = params.prefix
outdir = params.outDir
threads = params.threads
genomeSize = params.genomeSize
assembly_type = params.assembly_type.toLowerCase()
ref_genome = (params.ref_genome) ? file(params.ref_genome) : ''

/*

                    PARSING YAML FILE

*/

import org.yaml.snakeyaml.Yaml

//Def method for addtional parameters
class MyClass {
def getAdditional(String file, String value) {
  def yaml = new Yaml().load(new FileReader("$file"))
  def output = ""
  if ( "$value" == "canu" ) {
    yaml."$value".each {
  	   def (k, v) = "${it}".split( '=' )
  	    if ((v ==~ /null/ ) || (v == "")) {} else {
  	       output = output + " " + "${it}"
  	}}
    return output
  } else {
  yaml."$value".each {
    def (k, v) = "${it}".split( '=' )
    if ( v ==~ /true/ ) {
      output = output + " --" + k
      } else if ( v ==~ /false/ ) {}
        else if ((v ==~ /null/ ) || (v == "")) {} else {
          if ( k ==~ /k/ ) { output = output + " -" + k + " " + v }
          else { output = output + " --" + k + " " + v }
    }}
    return output
  }}}

def additionalParameters = [:]

if (params.yaml) {

//Creating map for additional parameters
additionalParameters['Spades'] = new MyClass().getAdditional(params.yaml, 'spades')
additionalParameters['Unicycler'] = new MyClass().getAdditional(params.yaml, 'unicycler')
additionalParameters['Canu'] = new MyClass().getAdditional(params.yaml, 'canu')
additionalParameters['Pilon'] = new MyClass().getAdditional(params.yaml, 'pilon')
additionalParameters['Flye'] = new MyClass().getAdditional(params.yaml, 'flye')

} else {
// Empty Map
additionalParameters['Spades'] = ""
additionalParameters['Unicycler'] = ""
additionalParameters['Canu'] = ""
additionalParameters['Pilon'] = ""
additionalParameters['Flye'] = ""
}

/*
 * Header log info
 */
log.info "========================================="
log.info "     Docker-based assembly Pipeline      "
log.info "========================================="
def summary = [:]
summary['Long Reads']   = params.longreads
summary['Fast5 files dir']   = params.fast5Path
summary['Long Reads']   = params.longreads
summary['Short single end reads']   = params.shortreads_single
summary['Short paired end reads']   = params.shortreads_paired
summary['Fasta Ref']    = params.ref_genome
summary['Output dir']   = params.outDir
summary['Assembly assembly_type chosen'] = params.assembly_type
summary['Long read sequencing technology'] = params.lr_type
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Command used']   = "$workflow.commandLine"
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*

                  PIPELINE BEGINS - Long Reads Only Assembly

*/

// Loading long reads files
canu_lreads = (params.longreads && (params.try_canu) && params.assembly_type == 'longreads-only') ?
              file(params.longreads) : Channel.empty()
unicycler_lreads = (params.longreads && (params.try_unicycler) && params.assembly_type == 'longreads-only') ?
                   file(params.longreads) : Channel.empty()
flye_lreads = (params.longreads && (params.try_flye) && params.assembly_type == 'longreads-only') ?
                   file(params.longreads) : Channel.empty()

// Check if fast5 dir path is set
if (params.fast5Path && params.assembly_type == 'longreads-only') {
  fast5 = Channel.fromPath( params.fast5Path )
  nanopolish_lreads = file(params.longreads)
  fast5_dir = Channel.fromPath( params.fast5Path, type: 'dir' )
} else { Channel.empty().into{fast5; fast5_dir; nanopolish_lreads} }

// Checking Long Reads technology
lreads_outdir = (params.lr_type == 'nanopore') ? 'ONT' : 'Pacbio'

/*

          Canu Assembler - Long reads only

*/

process canu_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  file lreads from canu_lreads

  output:
  file "${lreads_outdir}/canu_${lrID}/"
  file("${lreads_outdir}/canu_${lrID}/*.contigs.fasta") into canu_contigs

  when:
  (params.try_canu) && assembly_type == 'longreads-only'

  script:
  lr = (params.lr_type == 'nanopore') ? '-nanopore-raw' : '-pacbio-raw'
  lrID = lreads.getSimpleName()
  """
  canu -p ${prefix} -d ${lreads_outdir}/canu_${lrID} maxThreads=${params.threads}\
  genomeSize=${genomeSize} ${additionalParameters['Canu']} $lr $lreads
  """
}

/*

        Unicycler Assembler - Long reads only

*/

process unicycler_longReads {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  file lreads from unicycler_lreads

  output:
  file "${lreads_outdir}/unicycler_${lrID}/"
  file("${lreads_outdir}/unicycler_${lrID}/assembly.fasta") into unicycler_longreads_contigs

  when:
  (params.try_unicycler) && assembly_type == 'longreads-only'

  script:
  lrID = lreads.getSimpleName()
  """
  unicycler -l $lreads \
  -o ${lreads_outdir}/unicycler_${lrID} -t ${params.threads} \
  ${additionalParameters['Unicycler']} &> unicycler.log
  """
}

/*

            Flye Assembler - Long reads only

*/

process flye_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  file lreads from flye_lreads

  output:
  file "${lreads_outdir}/flye_${lrID}"
  file("${lreads_outdir}/flye_${lrID}/scaffolds.fasta") optional true
  file("${lreads_outdir}/flye_${lrID}/assembly_flye.fasta") into flye_contigs

  when:
  (params.try_flye) && assembly_type == 'longreads-only'

  script:
  lr = (params.lr_type == 'nanopore') ? '--nano-raw' : '--pacbio-raw'
  lrID = lreads.getSimpleName()
  """
  source activate flye ;
  mkdir ${lreads_outdir}/ ;
  flye ${lr} $lreads --genome-size ${genomeSize} --out-dir ${lreads_outdir}/flye_${lrID} \
  --threads $threads ${additionalParameters['Flye']} &> flye.log ;
  mv ${lreads_outdir}/flye_${lrID}/assembly.fasta ${lreads_outdir}/flye_${lrID}/assembly_flye.fasta
  """
}

/*
                                Managing Channels for Nanopolish and VariantCaller
*/

// You have fast5 and will execute nanopolish
if (params.fast5Path && params.lr_type == 'nanopore') {
    longread_assembly_nanopolish = Channel.empty().mix(flye_contigs, canu_contigs, unicycler_longreads_contigs)
    longread_assemblies_variantCaller = Channel.empty()
// You have h5 files and will execute VariantCaller
} else if (params.lr_type == 'pacbio' && params.pacbio_all_baxh5_path) {
  longread_assembly_nanopolish = Channel.empty()
  longread_assemblies_variantCaller = Channel.empty().mix(flye_contigs, canu_contigs, unicycler_longreads_contigs)
// You won't execute any of those
} else {
  longread_assembly_nanopolish = Channel.empty()
  longread_assemblies_variantCaller = Channel.empty()
}

/*

            NANOPOLISH - A tool to polish nanopore only assemblies

*/

process nanopolish {
  publishDir "${outdir}/nanopolish_output", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  each file(draft) from longread_assembly_nanopolish
  file(reads) from nanopolish_lreads
  file fast5
  val fast5_dir from fast5_dir

  output:
  file("${prefix}_${assembler}_nanopolished.fa") into nanopolished_contigs

  when:
  assembly_type == 'longreads-only' && (params.fast5Path)

  script:
  // Check assembler used
  if (draft.getName()  == 'assembly.fasta' || draft.getName() =~ /unicycler/) {
    assembler = 'unicycler'
  } else if (draft.getName()  == 'assembly_flye.fasta' || draft.getName() =~ /flye/) {
    assembler = 'flye'
  } else { assembler = 'canu' }
  """
  source activate NANOPOLISH ;
  zcat -f ${reads} > reads ;
  if [ \$(grep -c "^@" reads) -gt 0 ] ; then sed -n '1~4s/^@/>/p;2~4p' reads > reads.fa ; else mv reads reads.fa ; fi ;
  nanopolish index -d "${fast5_dir}" reads.fa ;
  minimap2 -d draft.mmi ${draft} ;
  minimap2 -ax map-ont -t ${params.threads} ${draft} reads.fa | samtools sort -o reads.sorted.bam -T reads.tmp ;
  samtools index reads.sorted.bam ;
  python /miniconda/envs/NANOPOLISH/bin/nanopolish_makerange.py ${draft} | parallel --results nanopolish.results -P ${params.cpus} \
  nanopolish variants --consensus -o polished.{1}.fa \
    -w {1} \
    -r reads.fa \
    -b reads.sorted.bam \
    -g ${draft} \
    --min-candidate-frequency 0.1;
  python /miniconda/envs/NANOPOLISH/bin/nanopolish_merge.py polished.*.fa > ${prefix}_${assembler}_nanopolished.fa
  """
}

/*

            VariantCaller - A pacbio only polishing step

*/
// Loading pacbio raw files
baxh5 = (params.pacbio_all_baxh5_path) ? Channel.fromPath(params.pacbio_all_baxh5_path).buffer( size: 3 ) : Channel.empty()
bams =  (params.pacbio_all_bam_path) ? Channel.fromPath(params.pacbio_all_bam_path).buffer( size: 3 ) : Channel.empty()

process bax2bam {
  publishDir "${outdir}/subreads", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  file(bax) from baxh5

  output:
  file "*.subreads.bam" into pacbio_bams

  when:
  params.lr_type == 'pacbio' && params.pacbio_all_baxh5_path != ''

  script:
  """
  source activate pacbio ;
  bax2bam ${bax.join(" ")} --subread  \
  --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag;
  """
}

// Get pacbio bams together
variantCaller_bams = Channel.empty().mix(pacbio_bams,bams).collect()

process variantCaller {
  publishDir "${outdir}/lreadsOnly_pacbio_consensus", mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  each file(draft) from longread_assemblies_variantCaller
  file bams from variantCaller_bams

  output:
  file "${prefix}_${assembler}_pbvariants.gff"
  file "${prefix}_${assembler}_pbconsensus.fasta" into variant_caller_contigs

  when:
  params.lr_type == 'pacbio' && ( params.pacbio_all_baxh5_path != '' || params.pacbio_all_bam_path != '' )

  script:
  // Check assembler used
  if (draft.getName()  == 'assembly.fasta' || draft.getName() =~ /unicycler/) {
    assembler = 'unicycler'
  } else if (draft.getName()  == 'assembly_flye.fasta' || draft.getName() =~ /flye/) {
    assembler = 'flye'
  } else { assembler = 'canu' }
  """
  source activate pacbio;
  for BAM in ${bams.join(" ")} ; do pbalign --nproc ${params.threads}  \
  \$BAM ${draft} \${BAM%%.bam}_pbaligned.bam; done;
  for BAM in *_pbaligned.bam ; do samtools sort -@ ${params.threads} \
  -o \${BAM%%.bam}_sorted.bam \$BAM; done;
  samtools merge pacbio_merged.bam *_sorted.bam;
  samtools index pacbio_merged.bam;
  pbindex pacbio_merged.bam;
  samtools faidx ${draft};
  arrow -j ${params.threads} --referenceFilename ${draft} -o ${prefix}_${assembler}_pbconsensus.fasta \
  -o ${prefix}_${assembler}_pbvariants.gff pacbio_merged.bam
  """

}

/*

                      HYBRID ASSEMBLY WITH Unicycler and Spades

*/

/*
                      SPAdes Assembler - Hybrid mode

*/
// Loading paired end short reads
short_reads_spades_hybrid_paired = (params.shortreads_paired && params.assembly_type == 'hybrid' \
                                    && (params.try_spades)) ?
                                    Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', ''])
// Loading single end short reads
short_reads_spades_hybrid_single = (params.shortreads_single && params.assembly_type == 'hybrid' \
                                    && (params.try_spades)) ?
                                    Channel.fromPath(params.shortreads_single) : ''
// Long reads
spades_hybrid_lreads = (params.longreads && params.assembly_type == 'hybrid' && (params.try_spades)) ?
                        file(params.longreads) : ''

// Assembly begin
process spades_hybrid_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus threads

  input:
  file lreads from spades_hybrid_lreads
  set val(id), file(sread1), file(sread2) from short_reads_spades_hybrid_paired
  file(sread) from short_reads_spades_hybrid_single
  file ref_genome from ref_genome

  output:
  file("spades_hybrid_results_${rid}/contigs.fasta") into spades_hybrid_contigs
  file "*"

  when:
  assembly_type == 'hybrid' && (params.try_spades)

  script:
  lr = (params.lr_type == 'nanopore') ? '--nanopore' : '--pacbio'
  spades_opt = (params.ref_genome) ? "--trusted-contigs $ref_genome" : ''

  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sread $lr $lreads"
    rid = sread.getSimpleName() + "_and_" + sread1.getSimpleName()
    x = "Executing assembly with paired and single end reads"
  } else if ((params.shortreads_single) && (params.shortreads_paired == '')) {
    parameter = "-s $sread $lr $lreads"
    rid = sread.getSimpleName()
    x = "Executing assembly with single end reads"
  } else if ((params.shortreads_paired) && (params.shortreads_single == '')) {
    parameter = "-1 $sread1 -2 $sread2 $lr $lreads"
    rid = sread1.getSimpleName()
    x = "Executing assembly with paired end reads"
  }
  """
  spades.py -o "spades_hybrid_results_${rid}" -t ${params.threads} ${additionalParameters['Spades']} \\
  $parameter ${spades_opt}
  """
}

/*

                          Unicycler Assembler - Hybrid mode

*/
// Loading paired end short reads
short_reads_unicycler_hybrid_paired = (params.shortreads_paired && params.assembly_type == 'hybrid' \
                                       && (params.try_unicycler)) ?
                                       Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', ''])
// Loading single end short reads
short_reads_unicycler_hybrid_single = (params.shortreads_single && params.assembly_type == 'hybrid' \
                                       && (params.try_unicycler)) ?
                                       Channel.fromPath(params.shortreads_single) : ''
// Long reads
unicycler_hybrid_lreads = (params.longreads && params.assembly_type == 'hybrid' && (params.try_unicycler)) ?
                          file(params.longreads) : ''

// Assembly begin
process unicycler_hybrid_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus threads

  input:
  set val(id), file(sread1), file(sread2) from short_reads_unicycler_hybrid_paired
  file(sread) from short_reads_unicycler_hybrid_single
  file lreads from unicycler_hybrid_lreads

  output:
  file "*"
  file("unicycler_hybrid_results_${rid}/assembly.fasta") into unicycler_hybrid_contigs

  when:
  assembly_type == 'hybrid' && (params.try_unicycler)

  script:
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sread -l $lreads"
    rid = sread.getSimpleName() + "_and_" + sread1.getSimpleName()
    x = "Executing assembly with paired and single end reads"
  } else if ((params.shortreads_single) && (params.shortreads_paired == '')) {
    parameter = "-s $sread -l $lreads"
    rid = sread.getSimpleName()
    x = "Executing assembly with single end reads"
  } else if ((params.shortreads_paired) && (params.shortreads_single == '')) {
    parameter = "-1 $sread1 -2 $sread2 -l $lreads"
    rid = sread1.getSimpleName()
    x = "Executing assembly with paired end reads"
  }
  """
  unicycler $parameter \\
  -o unicycler_hybrid_results_${rid} -t ${params.threads} \\
  ${additionalParameters['Unicycler']} &>unicycler.log
  """
}

/*

                          ILLUMINA-ONLY ASSEMBLY WITH Unicycler and Spades

*/
/*

                          SPAdes Assembler - Illumina only mode

*/
// Loading short reads
short_reads_spades_illumina_paired = (params.shortreads_paired && params.assembly_type == 'illumina-only' \
                                      && (params.try_spades)) ?
                                      Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', ''])
// Loading short reads
short_reads_spades_illumina_single = (params.shortreads_single && params.assembly_type == 'illumina-only' \
                                      && (params.try_spades)) ?
                                      Channel.fromPath(params.shortreads_single) : ''
// Assembly begin
process spades_illumina_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus threads

  input:
  set val(id), file(sread1), file(sread2) from short_reads_spades_illumina_paired
  file(sread) from short_reads_spades_illumina_single
  file ref_genome from ref_genome

  output:
  file("spades_illuminaOnly_results_${rid}/contigs.fasta") into spades_illumina_contigs
  file "*"

  when:
  assembly_type == 'illumina-only' && (params.try_spades)

  script:
  spades_opt = (params.ref_genome) ? "--trusted-contigs $ref_genome" : ''
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sread"
    rid = sread.getSimpleName() + "_and_" + sread1.getSimpleName()
    x = "Executing assembly with paired and single end reads"
  } else if ((params.shortreads_single) && (params.shortreads_paired == '')) {
    parameter = "-s $sread"
    rid = sread.getSimpleName()
    x = "Executing assembly with single end reads"
  } else if ((params.shortreads_paired) && (params.shortreads_single == '')) {
    parameter = "-1 $sread1 -2 $sread2"
    rid = sread1.getSimpleName()
    x = "Executing assembly with paired end reads"
  }
  """
  spades.py -o "spades_illuminaOnly_results_${rid}" -t ${params.threads} ${additionalParameters['Spades']} \\
  $parameter ${spades_opt}
  """
}

/*

                          Unicycler Assembler - Illumina only mode

*/
// Loading short reads
short_reads_unicycler_illumina_single = (params.shortreads_single && params.assembly_type == 'illumina-only' \
                                         && (params.try_unicycler)) ?
                                         Channel.fromPath(params.shortreads_single) : ''
short_reads_unicycler_illumina_paired = (params.shortreads_paired && params.assembly_type == 'illumina-only' \
                                         && (params.try_unicycler)) ?
                                         Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 ) : Channel.value(['', '', ''])
// Assembly begin
process unicycler_illumina_assembly {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  tag { x }
  cpus threads

  input:
  file(sread) from short_reads_unicycler_illumina_single
  set val(id), file(sread1), file(sread2) from short_reads_unicycler_illumina_paired

  output:
  file "*"
  file("unicycler_illuminaOnly_results_${rid}/assembly.fasta") into unicycler_illumina_contigs

  when:
  assembly_type == 'illumina-only' && (params.try_unicycler)

  script:
  if ((params.shortreads_single) && (params.shortreads_paired)) {
    parameter = "-1 $sread1 -2 $sread2 -s $sread"
    rid = sread.getSimpleName() + "_and_" + sread1.getSimpleName()
    x = "Executing assembly with paired and single end reads"
  } else if ((params.shortreads_single) && (params.shortreads_paired == '')) {
    parameter = "-s $sread"
    rid = sread.getSimpleName()
    x = "Executing assembly with single end reads"
  } else if ((params.shortreads_paired) && (params.shortreads_single == '')) {
    parameter = "-1 $sread1 -2 $sread2"
    rid = sread1.getSimpleName()
    x = "Executing assembly with paired end reads"
  }
  """
  unicycler $parameter \\
  -o unicycler_illuminaOnly_results_${rid} -t ${params.threads} \\
  ${additionalParameters['Unicycler']} &>unicycler.log
  """
}

/*
 * STEP 2 - ASSEMBLY POLISHING
 */

/*

                          Whenever the user have paired end short reads, this pipeline will execute
                          the polishing step with Unicycler polish pipeline.

                                                Unicycler Polishing Pipeline

*/

// You have already polished your genome with VariantCaller
if (params.pacbio_all_baxh5_path && (params.shortreads_paired) && params.illumina_polish_longreads_contigs == true) {
 Channel.empty().mix(variant_caller_contigs).set { unicycler_polish }
// You have already polished your genome with Nanopolish
} else if (params.fast5Path && (params.shortreads_paired) && params.illumina_polish_longreads_contigs == true) {
 Channel.empty().mix(nanopolished_contigs).set { unicycler_polish }
// You have not polished your genome whatsoever but want to polish it with short reads
} else if (params.pacbio_all_baxh5_path == '' && params.fast5Path == '' && (params.shortreads_paired) && params.illumina_polish_longreads_contigs == true) {
 Channel.empty().mix(flye_contigs, canu_contigs, unicycler_longreads_contigs).set { unicycler_polish }
// You have not polished your genome whatsoever and don't want to polish it with short reads
} else if (params.pacbio_all_baxh5_path == '' && params.fast5Path == '' && (params.shortreads_paired) && params.illumina_polish_longreads_contigs == false) {
 Channel.empty().set { unicycler_polish }
 Channel.empty().mix(flye_contigs, canu_contigs, unicycler_longreads_contigs).set { lreads_contigs }
} else {
 Channel.empty().set { unicycler_polish }
}

//Loading reads for polishing
short_reads_lreads_polish = (params.shortreads_paired) ? Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 )
                                                       : Channel.value(['', '', ''])
process illumina_polish_longreads_contigs {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  each file(draft) from unicycler_polish.collect()
  set val(id), file(sread1), file(sread2) from short_reads_lreads_polish

  output:
  file("${assembler}_lreadsOnly_exhaustive_polished")
  file("${assembler}_lreadsOnly_exhaustive_polished/${assembler}_final_polish.fasta") into unicycler_polished_contigs

  when:
  (assembly_type == 'longreads-only' && (params.illumina_polish_longreads_contigs) && (params.shortreads_paired))

  script:
  if (draft.getName()  == 'assembly.fasta' || draft.getName() =~ /unicycler/) {
    assembler = 'unicycler'
  } else if (draft.getName()  == 'assembly_flye.fasta' || draft.getName() =~ /flye/) {
    assembler = 'flye'
  } else { assembler = 'canu' }
  """
  mkdir ${assembler}_lreadsOnly_exhaustive_polished;
  unicycler_polish --ale /work/ALE/src/ALE --samtools /miniconda/bin/samtools --pilon /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
  -a $draft -1 $sread1 -2 $sread2 --threads $threads &> polish.log ;
  mv 0* polish.log ${assembler}_lreadsOnly_exhaustive_polished;
  mv ${assembler}_lreadsOnly_exhaustive_polished/*_final_polish.fasta ${assembler}_lreadsOnly_exhaustive_polished/${assembler}_final_polish.fasta;
  """
}

/*

                                    Whenever the user have unpaired short reads, this pipeline will execute
                                            the polishing step with a single Pilon round pipeline.

                                                              Pilon Polishing Pipeline

*/

// You have already polished your genome with VarianCaller
if (params.pacbio_all_baxh5_path && (params.shortreads_single) && params.illumina_polish_longreads_contigs == true) {
Channel.empty().mix(variant_caller_contigs).set { pilon_polish }
// You have already polished your genome with Nanopolish
} else if (params.fast5Path && (params.shortreads_single) && params.illumina_polish_longreads_contigs == true) {
Channel.empty().mix(nanopolished_contigs).set { pilon_polish }
// You have not polished your genome whatsoever but want to polish with unpaired short reads
} else if (params.pacbio_all_baxh5_path == '' && params.fast5Path == '' && (params.shortreads_single) && params.illumina_polish_longreads_contigs == true) {
Channel.empty().mix(flye_contigs, canu_contigs, unicycler_longreads_contigs).set { pilon_polish }
} else { Channel.empty().set { pilon_polish } }

//Load reads
short_reads_pilon_single = (params.shortreads_single) ?
                     Channel.fromPath(params.shortreads_single) : ''

process pilon_polish {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'
  cpus threads

  input:
  each file(draft) from pilon_polish.collect()
  file(sread) from short_reads_pilon_single

  output:
  file "pilon_results_${assembler}/pilon*"
  file("pilon_results_${assembler}/pilon*.fasta") into pilon_polished_contigs

  when:
  (assembly_type == 'longreads-only' && (params.illumina_polish_longreads_contigs) && (params.shortreads_single))

  script:
  parameter = "$sread"
  rid = sread.getSimpleName()
  x = "Polishing assembly with single end reads"

  if (draft.getName()  == 'assembly.fasta' || draft.getName() =~ /unicycler/) {
    assembler = 'unicycler'
  } else if (draft.getName()  == 'assembly_flye.fasta' || draft.getName() =~ /flye/) {
    assembler = 'flye'
  } else { assembler = 'canu' }
  """
  bwa index ${draft} ;
  bwa mem -M -t ${params.threads} ${draft} $parameter > ${rid}_${assembler}_aln.sam ;
  samtools view -bS ${rid}_${assembler}_aln.sam | samtools sort > ${rid}_${assembler}_aln.bam ;
  samtools index ${rid}_${assembler}_aln.bam ;
  java -Xmx${params.pilon_memmory_limit}G -jar /miniconda/share/pilon-1.23-2/pilon-1.23.jar \
  --genome ${draft} --bam ${rid}_${assembler}_aln.bam --output pilon_${assembler}_${rid} \
  --outdir pilon_results_${assembler} ${additionalParameters['Pilon']} &>pilon.log
  """
}

/*

                                                STEP 3 -  Assembly quality assesment with QUAST

*/
// If you have polished your genome with short reads
if (params.illumina_polish_longreads_contigs) {
  Channel.empty().mix(unicycler_polished_contigs, pilon_polished_contigs).set { final_assembly }
// You have polished your genome with VarianCaller but not with short reads
} else if (params.pacbio_all_baxh5_path && params.illumina_polish_longreads_contigs == false ) {
  Channel.empty().mix(variant_caller_contigs).set { final_assembly }
// You have polished your genome with Nanopolish but not with short reads
} else if (params.fast5Path && params.illumina_polish_longreads_contigs == false ) {
  Channel.empty().mix(nanopolished_contigs).set { final_assembly }
// You have not repolished your genome. Is a draft from any of the methods
} else {
  Channel.empty().mix(lreads_contigs, spades_hybrid_contigs, unicycler_hybrid_contigs,
    unicycler_illumina_contigs, spades_illumina_contigs).set { final_assembly }
}
//Loading reads for quast
short_reads_quast_single = (params.shortreads_single) ? Channel.fromPath(params.shortreads_single) : ''
short_reads_quast_paired = (params.shortreads_paired) ? Channel.fromFilePairs( params.shortreads_paired, flat: true, size: 2 )
                                                      : Channel.value(['', '', ''])
long_reads_quast = (params.longreads) ? Channel.fromPath(params.longreads) : ''

process quast {
  publishDir outdir, mode: 'copy'
  container 'fmalmeida/mpgap'

  input:
  each file(contigs) from final_assembly
  file 'reference_genome' from ref_genome
  file('sread') from short_reads_quast_single
  file('lreads') from long_reads_quast
  set val(id), file('pread1'), file('pread2') from short_reads_quast_paired

  output:
  file "quast_${type}_outputs_${assembler}/*"

  script:
  if ((params.shortreads_single) && (params.shortreads_paired) && assembly_type != 'longreads-only') {
    ref_parameter = "-M -t ${params.threads} reference_genome sread pread1 pread2"
    parameter = "-M -t ${params.threads} ${contigs} pread1 pread2"
    x = "Assessing assembly with paired and single end reads"
    sreads_parameter = "--single sread"
    preads_parameter = "--pe1 pread1 --pe2 pread2"
    lreads_parameter = ""
  } else if ((params.shortreads_single) && (params.shortreads_paired == '') && assembly_type != 'longreads-only') {
    ref_parameter = "-M -t ${params.threads} reference_genome sread"
    parameter = "-M -t ${params.threads} ${contigs} sread"
    x = "Assessing assembly with single end reads"
    sreads_parameter = "--single sread"
    preads_parameter = ""
    lreads_parameter = ""
  } else if ((params.shortreads_paired) && (params.shortreads_single == '') && assembly_type != 'longreads-only') {
    ref_parameter = "-M -t ${params.threads} reference_genome pread1 pread2"
    parameter = "-M -t ${params.threads} ${contigs} pread1 pread2"
    x = "Assessing assembly with paired end reads"
    sreads_parameter = ""
    preads_parameter = "--pe1 pread1 --pe2 pread2"
    lreads_parameter = ""
  } else if (assembly_type == 'longreads-only') {
    ltype = (params.lr_type == 'nanopore') ? "ont2d" : "pacbio"
    parameter = "-x ${ltype} -t ${params.threads} ${contigs} lreads"
    ref_parameter = "-x ${ltype} -t ${params.threads} reference_genome lreads"
    x = "Assessing assembly with long reads"
    sreads_parameter = ""
    preads_parameter = ""
    lreads_parameter = "--${params.lr_type} lreads"
  }
  if (contigs.getName()  == 'assembly.fasta' || contigs.getName() =~ /unicycler/) {
    assembler = 'unicycler'
  } else if (contigs.getName()  == 'contigs.fasta' || contigs.getName() =~ /spades/) {
    assembler = 'spades'
  } else if (contigs.getName()  == 'assembly_flye.fasta' || contigs.getName() =~ /flye/) {
    assembler = 'flye'
  } else { assembler = 'canu' }

  if (assembly_type == 'longreads-only') {
    type = 'lreadsOnly'
  } else if (assembly_type == 'illumina-only') {
    type = 'illuminaOnly'
  } else if (assembly_type == 'hybrid') {
    type = 'hybrid'
  }
  if (params.ref_genome != '')
  """
  source activate QUAST ;
  bwa index reference_genome ;
  bwa index ${contigs} ;
  bwa mem $parameter > contigs_aln.sam ;
  bwa mem $ref_parameter > reference_aln.sam ;
  quast.py -o quast_${type}_outputs_${assembler} -t ${params.threads} --ref-sam reference_aln.sam --sam contigs_aln.sam \\
  $sreads_parameter $preads_parameter $lreads_parameter -r reference_genome --circos ${contigs}
  """
  else
  """
  source activate QUAST ;
  bwa index ${contigs} ;
  bwa mem $parameter > contigs_aln.sam ;
  quast.py -o quast_${type}_outputs_${assembler} -t ${params.threads} --sam contigs_aln.sam \\
  $sreads_parameter $preads_parameter $lreads_parameter --circos ${contigs}
  """
}


// Completition message
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Execution duration: $workflow.duration"
}
