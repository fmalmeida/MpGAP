FROM ubuntu:16.04

MAINTAINER Felipe Marques de Almeida <180039822@aluno.unb.br>
SHELL ["/bin/bash", "-c"]

# SET DIR
WORKDIR /work

# Dependencies
RUN apt-get update && apt-get upgrade -y && \
	apt-get install -y git libc6-dev build-essential default-jdk gnuplot python libc6-dev gcc g++ make dpkg-dev python-pip python-pip-whl bowtie2 samtools ncbi-blast+ cmake wget python3 python3-pip

# Set miniconda
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
	bash /miniconda.sh -b -p /miniconda && \
	wget http://releases.llvm.org/6.0.0/clang+llvm-6.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz && \
	tar xf clang+llvm-6.0.0-x86_64-linux-gnu-ubuntu-16.04.tar.xz && \
	ln -rs clang+llvm-6.0.0-x86_64-linux-gnu-ubuntu-16.04/bin/* /usr/local/bin/
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

# Set Pilon
RUN conda install -y -c bioconda pilon

# Set Canu assembler
RUN wget https://github.com/marbl/canu/releases/download/v1.9/canu-1.9.Linux-amd64.tar.xz && \
tar -xJf canu-1.9.*.tar.xz && \
ln -rs canu-1.9/*/bin/* /usr/local/bin

# Set SPAdes assembler
RUN conda install -y -c bioconda spades

# Set Unicycler
RUN pip install --upgrade pip && \
	pip3 install --upgrade pip && \
	hash -d pip && \
	hash -d pip3 && \
	pip install setuptools && \
	pip3 install setuptools && \
	apt install zlib1g-dev && git clone https://github.com/rrwick/Unicycler.git && cd Unicycler && python3 setup.py install

## Set Unicycler complementary tools
RUN apt-get install -y mummer && \
	conda install -y -c bioconda samtools && \
	conda install -y -c bioconda bowtie2 && \
	conda install -y -c bioconda freebayes && \
	conda create -y -n legacy_python python=2.7 && conda install -y -c bioconda -n legacy_python pbalign && \
	conda install -y -n legacy_python -c bioconda bax2bam && \
	conda install -y -c conda-forge arrow && \
	git clone https://github.com/sc932/ALE.git && cd ALE/src && make && ln -rs ALE /usr/local/bin && \
	conda install -y -c bioconda minimap && \
	curl -L https://github.com/lh3/minimap2/releases/download/v2.10/minimap2-2.10_x64-linux.tar.bz2 | tar -jxvf - ; ln -rs minimap2-2.10_x64-linux/minimap2 /usr/local/bin && \
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
	conda install -y -c bioconda -n pacbio genomicconsensus pbalign blasr bax2bam pbcore pbcommand && \
	chmod a+x /miniconda/bin/*

# Set Flye
RUN conda create -y -n flye && \
	conda install -y -n flye -c bioconda flye

# Set Quast
RUN rm -rf /var/lib/apt/lists/*  && \
	apt-get clean && apt-get update && apt-get upgrade -y && \
	DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y hmmer unzip zlib1g-dev pkg-config libfreetype6-dev libpng-dev python-matplotlib && \
	got clone https://github.com/ablab/quast.git && \
	cd quast && ./setup.py install_full && \
	apt-get update && apt-get install -y bwa

# Set WorkDir
WORKDIR /work
