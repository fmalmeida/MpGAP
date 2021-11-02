include { write_csv } from '../nf_functions/writeCSV.nf'
workflow parse_samplesheet {

  take:
    data
    
  main:

    // iterate over input list
    custom_csv = write_csv(Channel.fromList(data))
  
    // now we parse the csv created
    parsed_csv = custom_csv.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'wtdbg2_technology','genomeSize', 'corrected_lreads', 'medaka_model','fast5', 'shasta_config','pacbio_bam']).map{ row ->

    if (row.entrypoint == 'shortreads_only') { 
      prefix = "${row.name}/shortreads_only" 
    }
    if (row.entrypoint == 'hybrid_strategy_1') {
      fixed_name = row.name - ":strategy_1"
      prefix = "${fixed_name}/hybrid_strategy_1" 
    }
    if (row.entrypoint == 'hybrid_strategy_2') { 
      fixed_name = row.name - ":strategy_2"
      prefix = "${fixed_name}/hybrid_strategy_2" 
    }
    if (row.entrypoint == 'longreads_only') { 
      prefix = "${row.name}/longreads_only" 
    }
    if (row.corrected_lreads && row.corrected_lreads.toLowerCase() != 'false') {
      corrected_lreads = 'true'
    } else {
      corrected_lreads = 'false'
    }

    // create input tuple   
    tuple(
        row.name, // sample name
        row.entrypoint, // assembly type
        (row.fwd == "missing_pairFWD")   ? row.fwd : file(row.fwd), // short reads pair 1
        (row.rev == "missing_pairREV")   ? row.rev : file(row.rev), // short reads pair 2
        (row.single == "missing_single") ? row.single : file(row.single), // short reads unpaired
        (row.lreads == "missing_lreads") ? row.lreads : file(row.lreads), // long reads
        row.lr_type, // long reads type
        row.wtdbg2_technology, // which wtdbg2 tech to use?
        row.genomeSize, // expected genome size
        row.corrected_lreads.toLowerCase(), // are reads corrected?
        row.medaka_model, // change medaka model?
        (row.fast5 == "missing_fast5") ? row.fast5 : file(row.fast5), // nanopolish fast5 as file
        row.shasta_config, // shasta config
        (row.pacbio_bam == "missing_pacbio_bam") ? row.pacbio_bam : file(row.pacbio_bam),
        prefix, // output prefix
      )
    }

    // now we create the filtered channels
    parsed_csv.branch{
      lronly: it[1] == 'longreads_only'
      sronly: it[1] == 'shortreads_only'
      hybrid: it[1] =~ /hybrid/
    }.set { results }

  emit:
  results.sronly
  results.lronly
  results.hybrid

}
