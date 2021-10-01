def parse_csv(in_ch) {
  
  return in_ch.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'wtdbg2_technology','genomeSize', 'corrected_lreads', 'medaka_model','fast5', 'pacbio_bams']).map{ row ->

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
        (row.fast5 == "missing_fast5") ? row.fast5 : Channel.fromPath(row.fast5), // nanopolish fast5 as file
        (row.fast5 == "missing_fast5") ? row.fast5 : Channel.fromPath(row.fast5, type: 'dir'), // nanopolish fast5 as dir value
        (row.pacbio_bams == "missing_pacbio_bams") ? row.pacbio_bams : Channel.fromPath(row.pacbio_bams).collect(), // pacbio bams as file
        (row.pacbio_bams == "missing_pacbio_bams") ? row.pacbio_bams : Channel.fromPath(row.pacbio_bams).count().subscribe { println it }, // pacbio bams as count
        prefix, // ouput prefix
    )
  }

}

def filter_ch(in_ch, entrypoint) {
  parse_csv(in_ch.map { it.text }) | filter { it[1] == "${entrypoint}" }
}
