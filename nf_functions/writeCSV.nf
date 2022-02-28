def write_csv(in_list) {

  return in_list.collectFile( name:"samples.txt", sort: { it[0] }, newLine: true ) {


                    // CHEKING FILE INPUTS

    /*
     * Checking the illumina input key
     */
    // defaults
    paired   = "missing_paired"
    fwd_pair = "missing_pairFWD"
    rev_pair = "missing_pairREV"
    single   = "missing_single"

    // illumina key is used for this sample?
    if (it.illumina) {
      if (it.illumina.size() == 1) {           // only unpaired reads
        single   = it.illumina[0]
      } else if (it.illumina.size() == 2) {    // only paired reads
        fwd_pair = it.illumina[0]
        rev_pair = it.illumina[1]
        paired   = 'has_paired'
      } else if (it.illumina.size() == 3) {    // both paired and unpaired reads
        fwd_pair = it.illumina[0]
        rev_pair = it.illumina[1]
        paired   = 'has_paired'
        single   = it.illumina[2]
      }
    }

    /*
     * Checking the long reads (ont or pacbio) input keys
     */
    // defaults
    lr_type = "missing_lr_type"
    lreads  = "missing_lreads"

    // ont or pacbio key are used for this sample?
    if (it.nanopore) {
      lr_type = "nanopore"
      lreads  = it.nanopore
    } else if (it.pacbio) {
      lr_type = "pacbio"
      lreads  = it.pacbio
    } else if (it.nanopore && it.pacbio) {
      println """
      ERROR!
      A major error has occurred with sample: ${it.id}
        ==> The pipeline pipeline is not yet capable of assembling both nanopore and pacbio reads together.
        Please, use only one of them for the sample. The only possibility of hybrid assemblies is between
        nanopore + illumina, or pacbio + illumina.
      Cheers.
      """.stripIndent()
      exit 1
    }

    /*
     * Checking for nanopolish_fast5 directory
     */
    nanopolish_fast5 = (it.nanopolish_fast5) ? it.nanopolish_fast5 : "missing_nanopolish_fast5"

    /*
     * Check for pacbio bams
     */
    pacbio_bam = (it.pacbio_bam) ? it.pacbio_bam : "missing_pacbio_bam"


                  // CHECKING FOR ASSEMBLY CONFIGURATIONS

     
    /*
     * Checking if long reads are corrected
     */
    // defaults
    corrected_long_reads = (params.corrected_long_reads) ? true : false

    // corrected_long_reads input key is used for the sample?
    if (it.corrected_long_reads) {
      corrected_long_reads = it.corrected_long_reads
      if (corrected_long_reads.toString().toLowerCase() != "true" && corrected_long_reads.toString().toLowerCase() != "false") {
        println """
        ERROR!
        A minor error has occurred!
          ==> In the YAML, the 'corrected_long_reads:' must be either true or false.
        Please the re-check the parameters. Problem in sample: ${it.id}.
        Cheers.
        """.stripIndent()
        exit 1
      }
    }

    /*
     * Checking for nanopolish max haplotypes
     */
    // defaults -- are params given outside YAML
    nanopolish_max_haplotypes = params.nanopolish_max_haplotypes

    // nanopolish_max_haplotypes input key is used for the sample? -- it overwrites params for that sample
    if (it.nanopolish_max_haplotypes) { nanopolish_max_haplotypes = it.nanopolish_max_haplotypes }

    /*
     * Check for medaka model
     */
    // defaults -- are params given outside YAML
    medaka_model = params.medaka_model

    // medaka_model input key is used for the sample? -- it overwrites params for that sample
    if (it.medaka_model) { medaka_model = it.medaka_model }

    /*
     * Check for shasta config
     */
    // defaults -- are params given outside YAML
    shasta_config = params.shasta_config

    // shasta_config input key is used for the sample? -- it overwrites params for that sample
    if (it.shasta_config) { shasta_config = it.shasta_config }

    /*
     * Check for sample genome_size
     */
    // defaults -- are params given outside YAML
    genome_size = params.genome_size

    // genome_size input key is used for the sample? -- it overwrites params for that sample
    if (it.genome_size) { genome_size = it.genome_size }

    /*
     * Check for wtdbg2 technology
     */
    // wtdbg2 default configuration -- creating one automatically
    if (lr_type == "nanopore") {
        wtdbg2_technology = "ont"
    } else if (lr_type == "pacbio") {
        wtdbg2_technology = "sq"
    } else {
        wtdbg2_technology = "not_used"
    }

    // overwrite defaults with params given outside YAML
    if (params.wtdbg2_technology) { wtdbg2_technology = params.wtdbg2_technology }

    // wtdbg2_technology input key is used for the sample? -- it overwrites params for that sample
    if (it.wtdbg2_technology) { wtdbg2_technology = it.wtdbg2_technology }

    /*
     * Check for desired hybrid strategy model
     */
    // defaults -- are params given outside YAML
    hybrid_strategy = params.hybrid_strategy.toString().toLowerCase()

    // hybrid_strategy input key is used for the sample?
    if (it.hybrid_strategy) { hybrid_strategy = it.hybrid_strategy.toString().toLowerCase() }

    /*
     * Create entrypoints
     */
    // any short reads are given?
    has_short_reads = (paired != "missing_paired" || single != "missing_single") ? true : false

    // any long reads are given?
    has_long_reads = (lr_type != "missing_lr_type" && lreads != "missing_lreads") ? true : false

    // has only short reads
    if (has_short_reads && !has_long_reads) { entrypoint = 'shortreads_only' }

    // has both short reads and long reads -- hybrid
    else if (has_short_reads && has_long_reads) {
      entrypoint = 'hybrid_strategy_' + hybrid_strategy.toString().toLowerCase()

      // check if the hybrid strategy is valid
      if (entrypoint != 'hybrid_strategy_1' && entrypoint != 'hybrid_strategy_2' && entrypoint != 'hybrid_strategy_both') {
        println """
        ERROR!
        A major error has occurred!
          ==> The 'hybrid_strategy:' key in the YAML, or the --hybrid_strategy command line, must be either '1', '2' or 'both'.
        Please the re-check the parameters. Problem in sample: ${it.id}.
        Cheers.
        """.stripIndent()
        exit 1
      }
    }
       
    // has only long reads
    else if (!has_short_reads && has_long_reads) { entrypoint = 'longreads_only' } 
       
    // has nothing
    else {
      println """
      ERROR!
      A major error has occurred!
        ==> At least one read type must be set: illumina, nanopore or pacbio!
      Please, re-check sample: ${it.id}.
      Cheers.
      """.stripIndent()
      exit 1
    }

    /*
     * Paramaters check
     * check whether the input parameters for a sample are enough for analysis
     */
    // genome size is given if using canu?
    if (!params.skip_canu && !genome_size && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2' || entrypoint == 'hybrid_strategy_both')) {
      println """
      ERROR!
      A minor error has occurred
        ==> The pipeline will try to run canu assembler on sample ${it.id}, but user forgot to tell the expected genome size.
      Please either set it with --genome_size or inside the YAML (genome_size:) for that sample, or turn off the canu assembly with --skip_canu.
      Cheers.
      """.stripIndent()
      exit 1
    }

    // genome size is given if using haslr?
    if (!params.skip_haslr && !genome_size && (entrypoint == 'hybrid_strategy_1' || entrypoint == 'hybrid_strategy_both')) {
      println """
      ERROR!
      A minor error has occurred
        ==> The pipeline will try to run haslr assembler on sample ${it.id}, but user forgot to tell the expected genome size.
      Please either set it with --genome_size or inside the YAML (genome_size:) for that sample, or turn off the haslr assembly with --skip_haslr
      Cheers.
      """.stripIndent()
      exit 1
    }

    // genome size is given if using wtdbg2?
    if (!params.skip_wtdbg2 && !genome_size && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2' || entrypoint == 'hybrid_strategy_both')) {
      println """
      ERROR!
      A minor error has occurred
        ==> The pipeline will try to run wtdbg2 assembler on sample ${it.id}, but user forgot to tell the expected genome size.
      Please either set it with --genome_size or inside the YAML (genome_size:) for that sample, or turn off the haslr assembly with --skip_wtdbg2
      Cheers.
      """.stripIndent()
      exit 1
    }

    // Checking if shovill parameter contain "--assembler"
    if (params.shovill_additional_parameters =~ /--assembler/) {
      println """
      ERROR!
      A minor error has occurred
        ==> The pipeline already executes shovill with both spades, skesa and megahit. Therefore, you can't pass the 
      "--assembler" parameter as an additional parameter for shovill inside "--shovill_additional_parameters".
      Please, remove the shovill additional parameter "--assembler" and re-run the pipeline.
      Cheers.
      """.stripIndent()
      exit 1
    }



    /*
     * Output samplesheet as CSV
     */
    if (entrypoint == 'hybrid_strategy_both') { // creates two lines for the sample, for both hybrid strategies
      "${it.id}:strategy_1,hybrid_strategy_1,${fwd_pair},${rev_pair},${single},${lreads},${lr_type},${wtdbg2_technology},${genome_size},${corrected_long_reads},${medaka_model},${nanopolish_fast5},${nanopolish_max_haplotypes},${shasta_config},${pacbio_bam}\n${it.id}:strategy_2,hybrid_strategy_2,${fwd_pair},${rev_pair},${single},${lreads},${lr_type},${wtdbg2_technology},${genome_size},${corrected_long_reads},${medaka_model},${nanopolish_fast5},${nanopolish_max_haplotypes},${shasta_config},${pacbio_bam}"
    } else {
      "${it.id},${entrypoint},${fwd_pair},${rev_pair},${single},${lreads},${lr_type},${wtdbg2_technology},${genome_size},${corrected_long_reads},${medaka_model},${nanopolish_fast5},${nanopolish_max_haplotypes},${shasta_config},${pacbio_bam}"
    }
  
  } // end of collectFile function

} // end of write_csv function
