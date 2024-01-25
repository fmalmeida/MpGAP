
Advanced assembler customization options
""""""""""""""""""""""""""""""""""""""""

.. note::

  Additional parameters must be set inside double quotes separated by blank spaces.

.. list-table::
   :widths: 30 10 60
   :header-rows: 1

   * - Arguments
     - Default value
     - Description

   * - ``--quast_additional_parameters``
     - NA
     - | Give additional parameters to Quast while assessing assembly metrics. Must be given as shown in Quast manual. E.g. ``" --large --eukaryote "``.
   
   * - ``--skip_raw_assemblies_polishing``
     - false
     - | This will make the pipeline not polish raw assemblies on hybrid strategy 2.
       | For example, if a sample is assembled with flye and polished with medaka, by default, both assemblies will be passed to pilon so you can compare them.
       | If you don't need this comparison and don't want to polish the raw assembly, use this parameter.

   * - ``--skip_canu``
     - false
     - Skip the execution of Canu

   * - ``--canu_additional_parameters``
     - NA
     - | Passes additional parameters for Canu assembler. E.g. ``" correctedErrorRate=0.075 corOutCoverage=200 "``. Must be given as shown in Canu's manual.

   * - ``--skip_flye``
     - false
     - Skip the execution of Flye

   * - ``--flye_additional_parameters``
     - NA
     - | Passes additional parameters for Flye assembler. E.g. ``" --meta --iterations 4 "``. Must be given as shown in Flye's manual.

   * - ``--skip_raven``
     - false
     - Skip the execution of Raven

   * - ``--raven_additional_parameters``
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``" --polishing-rounds 4 "``. Must be given as shown in Raven's manual.
   
   * - ``--skip_shasta``
     - false
     - Skip the execution of Shasta

   * - ``--shasta_additional_parameters``
     - NA
     - | Passes additional parameters for Raven assembler. E.g. ``" --Assembly.detangleMethod 1 "``. Must be given as shown in Shasta's manual.
   
   * - ``--skip_wtdbg2``
     - false
     - Skip the execution of Raven

   * - ``--wtdbg2_additional_parameters``
     - NA
     - | Passes additional parameters for wtdbg2 assembler. E.g. ``" -k 250 "``. Must be given as shown in wtdbg2's manual. Remember, the script called for wtdbg2 is ``wtdbg2.pl`` thus you must give the parameters used by it.

   * - ``--skip_unicycler``
     - false
     - Skip the execution of Unicycler

   * - ``--unicycler_additional_parameters``
     - NA
     - | Passes additional parameters for Unicycler assembler. E.g. ``" --mode conservative --no_correct "``. Must be given as shown in Unicycler's manual.

   * - ``--skip_spades``
     - false
     - Skip the execution of SPAdes

   * - ``--spades_additional_parameters``
     - NA
     - | Passes additional parameters for SPAdes assembler. E.g. ``" --meta --plasmids "``. Must be given as shown in Spades' manual.

   * - ``--skip_haslr``
     - false
     - Skip the execution of Haslr

   * - ``--haslr_additional_parameters``
     - NA
     - | Passes additional parameters for Haslr assembler. E.g. ``" --cov-lr 30 "``. Must be given as shown in Haslr' manual.

   * - ``--skip_shovill``
     - false
     - Skip the execution of Shovill

   * - ``--shovill_additional_parameters``
     - NA
     - | Passes additional parameters for Shovill assembler. E.g. ``" --depth 15 "``. Must be given as shown in Shovill's manual.
       | The pipeline already executes shovill with spades, skesa and megahit, so please, do not use it with shovill's ``--assembler`` parameter.
   
   * - ``--skip_megahit``
     - false
     - Skip the execution of Megahit

   * - ``--megahit_additional_parameters``
     - NA
     - | Passes additional parameters for Megahit assembler. E.g. ``" --presets meta-large "``. Must be given as shown in Megahit's manual.
