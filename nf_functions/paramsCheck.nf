
def paramsCheck() {

  // genome size is given if using canu?
  if (!params.skip_canu && (genomeSize == "missing_genomeSize") && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2' || entrypoint == 'hybrid_strategy_both')) {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to run canu assembler on sample ${it.id}, but user forgot to tell the expected genome size.
    Please either set it with --genomeSize or inside the YAML (genomeSize:) for that sample, or turn off the canu assembly with --skip_canu.
    Cheers.
    """.stripIndent()
    exit 1
  }

  // genome size is given if using haslr?
  if (!params.skip_haslr && genomeSize == "missing_genomeSize" && (entrypoint == 'hybrid_strategy_1' || entrypoint == 'hybrid_strategy_both')) {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to run haslr assembler on sample ${it.id}, but user forgot to tell the expected genome size.
    Please either set it with --genomeSize or inside the YAML (genomeSize:) for that sample, or turn off the haslr assembly with --skip_haslr
    Cheers.
    """.stripIndent()
    exit 1
  }

  // genome size is given if using wtdbg2?
  if (!params.skip_wtdbg2 && genomeSize == "missing_genomeSize" && (entrypoint == 'longreads_only' || entrypoint == 'hybrid_strategy_2' || entrypoint == 'hybrid_strategy_both')) {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to run wtdbg2 assembler on sample ${it.id}, but user forgot to tell the expected genome size.
    Please either set it with --genomeSize or inside the YAML (genomeSize:) for that sample, or turn off the haslr assembly with --skip_wtdbg2
    Cheers.
    """.stripIndent()
    exit 1
  }

  // wtdgb2 is used with pacbio, and the technology is not properly set?
  if (!params.skip_wtdbg2 && lr_type == 'pacbio' && wtdbg2_technology == "ont") {
    println """
    ERROR!
    A minor error has occurred
      ==> The pipeline will try to assemble pacbio reads from sample ${it.id} with wtdbg2, but user forgot to tell the correct techonology. Options: "rs" for PacBio RSII, "sq" for PacBio Sequel, "ccs" for PacBio CCS reads.
    Please either set it with --wtdbg2_technology or inside the YAML (wtdbg2_technology:) for that sample, or turn off the wtdbg2 assembly with --skip_wtdbg2.
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

} // end of paramsCheck function
