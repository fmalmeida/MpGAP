def parse_csv(in_ch) {

  return in_ch.splitCsv(header: ['name', 'entrypoint', 'fwd', 'rev', 'single', 'lreads', 'lr_type', 'fast5']).map{ row ->

    if (row.entrypoint == "sr-only") {
      tuple(
        row.name, // prefix
        row.entrypoint, // assembly type
        (row.fwd == "missing_pairFWD") ? row.fwd : file(row.fwd), // short reads pair 1
        (row.rev == "missing_pairREV") ? row.rev : file(row.rev), // short reads pair 2
        (row.single == "missing_single") ? row.single : file(row.single), // short reads unpaired
        "shortreads_only/${row.name}",
      )
    }
  }

}

def filter_ch(in_ch, entrypoint) {
  parse_csv(in_ch.map { it.text }) | filter { it[1] == "${entrypoint}" }
}
