# Pull base image
FROM phusion/baseimage:0.11

# Who maintains this image
LABEL maintainer Jason Smith "jasonsmith@virginia.edu"

# Version info
LABEL version 0.7

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes \
    autoconf \
    automake \
    autotools-dev \
    curl \
    default-jre \
    default-jdk \
    git \
    gsl-bin \
    libgsl-dbg \
    libgsl-dev \
    libcommons-math3-java \
    libcurl4-gnutls-dev \ 
    libjbzip2-java \
    libpng-dev \
    libssl-dev \
    libtbb2 \
    libtbb-dev \
    libtool \
    openssl \
    pigz \
    python \
    python-pip python-dev build-essential \
    python3 \
    python3-pip \
    wget

# Install MySQL server
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes mysql-server \
    mysql-client \
    libmysqlclient-dev
    
# Install python tools
RUN pip install --upgrade pip
RUN pip install cutadapt && \
    pip install numpy && \
    pip install https://github.com/pepkit/looper/zipball/master && \
    pip install pararead && \
    pip install pandas && \
    pip install piper && \
    pip install refgenie
    

RUN pip3 install --upgrade pip
RUN pip3 install cutadapt && \
    pip3 install numpy && \
    pip3 install https://github.com/pepkit/looper/zipball/master && \
    pip3 install pararead && \
    pip3 install pandas && \
    pip3 install piper && \
    pip3 install refgenie

# Install R
RUN echo '# Ubuntu 18.04 (Bionic) - R' >> /etc/apt/sources.list && \
    echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get --assume-yes install r-base r-base-dev && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org/'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('devtools')" && \
    Rscript -e "devtools::install_github('pepkit/pepr')" && \
    Rscript -e "install.packages('BiocManager')" && \
    Rscript -e "BiocManager::install(); BiocManager::install('GenomicRanges')" && \
    Rscript -e "devtools::install_github('databio/GenomicDistributions')"
ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true
RUN Rscript -e "devtools::install_github('databio/peppro', subdir='PEPPROr')"

# Install htslib
WORKDIR /home/src/
RUN wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2 && \
    tar xf htslib-1.7.tar.bz2 && \
    cd /home/src/htslib-1.7 && \
    ./configure --prefix /home/tools/ && \
    make && \
    make install

# Install samtools
WORKDIR /home/src/
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 && \
    tar xf samtools-1.7.tar.bz2 && \
    cd /home/src/samtools-1.7 && \
    ./configure && \
    make && \
    make install && \
    ln -s /home/src/samtools-1.7/samtools /usr/bin/

# Install bedtools
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes \
    ant \
    bedtools

# Install bowtie2
WORKDIR /home/src/
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.5.1/bowtie2-2.3.5.1-source.zip && \
    unzip bowtie2-2.3.5.1-source.zip && \
    cd /home/src/bowtie2-2.3.5.1 && \
    make && \
    make install && \
    ln -s /home/src/bowtie2-2.3.5.1/bowtie2 /usr/bin/

# Install seqkit
WORKDIR /home/src/
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.10.1/seqkit_linux_amd64.tar.gz && \
    tar -zxvf seqkit_linux_amd64.tar.gz && \
    ln -s /home/src/seqkit /usr/bin/

# Install fastp
WORKDIR /home/src/
RUN git clone https://github.com/OpenGene/fastp.git && \
    cd fastp && \
    make && \
    make install && \
    ln -s /home/src/fastp/fastp /usr/bin/

# Install seqtk
WORKDIR /home/src/
RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make && \
    ln -s /home/src/seqtk/seqtk /usr/bin/

# Install preseq
WORKDIR /home/src/
RUN wget http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2 && \
    tar -jxvf preseq_linux_v2.0.tar.bz2 && \
    rm preseq_linux_v2.0.tar.bz2 && \
    cd preseq_v2.0/ && \
    ln /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgsl.so.0 && \
    ln -s /home/src/preseq_v2.0/preseq /usr/bin/ && \
    ln -s /home/src/preseq_v2.0/bam2mr /usr/bin/

# Install UCSC tools
WORKDIR /home/tools/
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat && \
    chmod +x /home/tools/wigToBigWig && \
    chmod +x /home/tools/bigWigCat && \
    ln -s /home/tools/wigToBigWig /usr/bin/ && \
    ln -s /home/tools/bigWigCat /usr/bin/


# OPTIONAL REQUIREMENTS
# Install fastqc
WORKDIR /home/tools/
RUN wget https://github.com/s-andrews/FastQC/archive/v0.11.8.zip && \
    unzip v0.11.8.zip && \
    cd /home/tools/FastQC-0.11.8 && \
    chmod 755 fastqc && \ 
    ln -s /home/tools/FastQC-0.11.8/fastqc /usr/bin/

# Install cargo
WORKDIR /home/tools/
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
RUN echo 'source $HOME/.cargo/env' >> $HOME/.bashrc

# Install fqdedup
WORKDIR /home/tools/
RUN git clone https://github.com/guertinlab/fqdedup.git && \
    cd fqdedup && \
    bash -c 'source $HOME/.cargo/env; cargo build --release' && \
    ln -s /home/tools/fqdedup/target/release/fqdedup /usr/bin/

# Install fastx_toolkit
WORKDIR /home/tools/
RUN git clone https://github.com/agordon/libgtextutils.git && \
    cd libgtextutils && \
    ./reconf && \
    ./configure && \
    make && \
    make install
WORKDIR /home/tools/
RUN git clone https://github.com/agordon/fastx_toolkit && \
    cd fastx_toolkit && \
    ./reconf && \
    ./configure && \
    find . -type f -exec sed -ie 's/-Werror//g' {} \; && \
    make && \
    make install

# Install seqOutBias
WORKDIR /home/tools/
RUN git clone https://github.com/guertinlab/seqOutBias.git && \
    cd seqOutBias && \
    bash -c 'source $HOME/.cargo/env; cargo build --release' && \
    ln -s /home/tools/seqOutBias/target/release/seqOutBias /usr/bin/


# Set environment variables
ENV PATH=/home/tools/bin:/home/tools/:/home/src/bowtie2-2.3.5.1:/home/src/samtools-1.7:/home/src/htslib-1.7:$PATH \
    R_LIBS_USER=/usr/local/lib/R/site-library/ \
    PYTHONPATH=/usr/local/lib/python3.6/dist-packages:$PYTHONPATH

# Define default command
WORKDIR /home/
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*