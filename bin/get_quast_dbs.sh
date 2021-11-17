#!/bin/bash -ue

mkdir -p $CONDA_PREFIX/lib/python3.6/site-packages/quast_libs/busco/
### bacteria db
BACTERIA_DB="$CONDA_PREFIX/lib/python3.6/site-packages/quast_libs/busco/bacteria.tar.gz"
[ ! -f "$BACTERIA_DB" ] && wget -O "$BACTERIA_DB" https://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
### eukaryota db
EUKARYOTA_DB="$CONDA_PREFIX/lib/python3.6/site-packages/quast_libs/busco/eukaryota.tar.gz"
[ ! -f "$EUKARYOTA_DB" ] && wget -O "$EUKARYOTA_DB" https://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
### fungi db
FUNGI_DB="$CONDA_PREFIX/lib/python3.6/site-packages/quast_libs/busco/fungi.tar.gz"
[ ! -f "$FUNGI_DB" ] &&  wget -O "$FUNGI_DB" https://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
