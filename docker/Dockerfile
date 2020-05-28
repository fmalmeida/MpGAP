FROM ubuntu:18.04

MAINTAINER Felipe Marques de Almeida <180039822@aluno.unb.br>
SHELL ["/bin/bash", "-c"]

# SET DIR
WORKDIR /work

# Dependencies
RUN apt-get update && apt-get upgrade -y && \
	apt-get install -y git libc6-dev build-essential default-jdk gnuplot python libc6-dev gcc g++ make dpkg-dev python-pip python-pip-whl bowtie2 samtools ncbi-blast+ cmake wget python3 python3-pip

# Set miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
	bash /miniconda.sh -b -p /miniconda

ENV PATH="/miniconda/bin:$PATH"

# Set ncbi-blast+
RUN apt-get install -y build-essential wget git hmmer python perl bedtools libidn11 && \
        wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
        tar -zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz ; cd ncbi-blast-2.9.0+/bin/ ; ln -rs $(ls) /usr/local/bin/

# INSTALL Racon
RUN git clone --recursive https://github.com/isovic/racon.git racon && \
	cd racon && mkdir build && cd build && \
	cmake -DCMAKE_BUILD_TYPE=Release .. && \
	make && make install

# INSTALL Medaka
RUN conda create -y -n MEDAKA && \
	conda install -y -n MEDAKA -c conda-forge -c bioconda medaka

# Set Pilon
RUN conda install -y -c bioconda pilon

# Set Canu assembler
RUN wget https://github.com/marbl/canu/releases/download/v1.9/canu-1.9.Linux-amd64.tar.xz && \
tar -xJf canu-1.9.*.tar.xz && \
ln -rs canu-1.9/*/bin/* /usr/local/bin

# Set SPAdes assembler
RUN conda install -y -c bioconda spades

# Set Unicycler
RUN apt-get install -y zlib1g zlib1g-dev && \
		conda install -y -c bioconda unicycler

## Set Unicycler complementary tools
## Only Illumina polish
RUN git clone https://github.com/sc932/ALE.git && \
		cd ALE/src && make && ln -rs ALE /usr/local/bin && \
		ALE --help
RUN conda install -y -c bioconda minimap2 && \
		chmod a+x /miniconda/bin/*

# Set Nanopolish
RUN apt-get install -y libhdf5-dev parallel && \
       conda config --add channels defaults && \
       conda config --add channels conda-forge && \
       conda config --add channels bioconda && \
       conda create -y -n NANOPOLISH && \
       conda install -y -n NANOPOLISH nanopolish

# Set PACBIO GENOMICONSENSUS
RUN conda create --name pacbio python=2 && \
	conda install -y -c bioconda -n pacbio bax2bam bam2fastx lima genomicconsensus pbalign blasr bax2bam pbcore pbcommand && \
	chmod a+x /miniconda/bin/*

# Set Flye
RUN conda create -y -n flye && \
	conda install -y -n flye -c bioconda flye

# Set Quast
RUN conda create -y -n QUAST && \
	conda install -y -n QUAST -c bioconda quast && \
	apt-get update && apt-get install -y bwa

# Check envs
RUN conda info -e

# Set WorkDir
WORKDIR /work
