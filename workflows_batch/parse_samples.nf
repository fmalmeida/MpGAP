include { write_csv } from '../nf_functions/writeCSV.nf'
workflow parse_samplesheet {

  take:
    data
    
  main:

    // iterate over input list
    custom_csv = write_csv(Channel.fromList(data))
  
    // now we parse the csv created
    parsed_csv = custom_csv.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'wtdbg2_technology','genomeSize', 'corrected_lreads', 'medaka_model','fast5', 'pacbio_bams']).map{ row ->

    if (row.entrypoint == 'sr-only') { prefix = "${row.name}/shortreads_only" }
    if (row.entrypoint == 'hybrid-strategy-1') { prefix = "${row.name}/hybrid_strategy_1" }
    if (row.entrypoint == 'hybrid-strategy-2') { prefix = "${row.name}/hybrid_strategy_2" }
    if (row.entrypoint == 'lr-only') { prefix = "${row.name}/longreads_only" }
    
    tuple(
        row.name, // sample name
        row.entrypoint, // assembly type
        (row.fwd == "missing_pairFWD") ? row.fwd : file(row.fwd), // short reads pair 1
        (row.rev == "missing_pairREV") ? row.rev : file(row.rev), // short reads pair 2
        (row.single == "missing_single") ? row.single : file(row.single), // short reads unpaired
        (row.lreads == "missing_lreads") ? row.lreads : file(row.lreads), // long reads
        row.lr_type, // long reads type
        row.wtdbg2_technology, // which wtdbg2 tech to use?
        row.genomeSize, // expected genome size
        row.corrected_lreads.toLowerCase(), // are reads corrected?
        row.medaka_model, // change medaka model?
        (row.fast5 == "missing_fast5") ? row.fast5 : file(row.fast5), // nanopolish fast5 as file
        (row.pacbio_bams == "missing_pacbio_bams") ? row.pacbio_bams : file(row.pacbio_bams).collect(),
        prefix, // ouput prefix
      )
    }

    // now we create the filtered channels
    parsed_csv.branch{
      lronly: it[1] == 'lr-only'
      sronly: it[1] == 'sr-only'
      hybrid: it[1] =~ /hybrid/
    }.set { results }

  emit:
  results.sronly
  results.lronly
  results.hybrid

}
