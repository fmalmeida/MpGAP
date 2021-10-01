workflow parse_samplesheet {

  take:
    data
  main:

    // Parse input to check for missing entries
    parsed = []

    data.each {

      /*
       * Check for illumina reads
       */
      if (it.illumina) {
        if (it.illumina.size() == 1) {
          it['paired'] = "missing_paired"
          it['single'] = it.illumina[0]
        } else if (it.illumina.size() == 2) {
          it['paired'] = [it.illumina[0], it.illumina[1]]
          it['single'] = "missing_single"
        } else if (it.illumina.size() == 3) {
          it['paired'] = [it.illumina[0], it.illumina[1]]
          it['single'] = it.illumina[2]
        }
      } else {
        it['paired'] = "missing_paired"
        it['single'] = "missing_single"
      }

      /*
       * Check for long reads
       */
      if (it.nanopore) {
        it['lr_type'] = "nanopore"
        it['lreads']  = it.nanopore
      } else if (it.pacbio) {
        it['lr_type'] = "pacbio"
        it['lreads']  = it.pacbio
      } else {
        it['lr_type'] = "missing_lr_type"
        it['lreads']  = "missing_lreads"
      }

      /*
       * Check if long reads are corrected
       */
      if (params.corrected_lreads || it.corrected_lreads) {
        it['corrected_lreads'] = 'true'
      } else {
        it['corrected_lreads'] = 'false'
      }

      /*
       * Check for fast5 directory
       */
      it['fast5'] = (it.fast5)   ? it.fast5  : "missing_fast5"

      /*
       * Check for medaka model
       */
      if (it.medaka_model) {
        it['medaka_model'] = it.medaka_model
      } else {
        it['medaka_model'] = params.medaka_sequencing_model
      }

      /*
       * Check for sample genomeSize
       */
      if (params.genomeSize) {
        it['genomeSize'] = params.genomeSize
      } else if (!params.genomeSize && it.genomeSize) {
        it['genomeSize'] = it.genomeSize
      } else {
        it['genomeSize'] = 'missing_genomeSize'
      }

      /*
       * Check for pacbio bams
       */
      it['pacbio_bams'] = (it.pacbio_bams) ? it.pacbio_bams : "missing_pacbio_bams"

      /*
       * Check for wtdbg2 technology
       */
      if (it.wtdbg2_technology) {
        it['wtdbg2_technology'] = it.wtdbg2_technology
      } else {
        it['wtdbg2_technology'] = params.wtdbg2_technology
      }

      // Save
      parsed.add(it)
    }

    emit:
      Channel.from(parsed)

}
