# Pull base image
FROM phusion/baseimage:latest

# Who maintains this image
LABEL maintainer Jason Smith "jasonsmith@virginia.edu"

# Version info
LABEL version 0.7

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install dependencies
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes \    
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
    openssl \
    pigz \
    python \
    python-pip python-dev build-essential \
    wget

# Install MySQL server
RUN DEBIAN_FRONTEND=noninteractive apt-get install --assume-yes mysql-server \
    mysql-client \
    libmysqlclient-dev
    
# Install python tools
RUN pip install --upgrade pip
RUN pip install numpy && \
    pip install https://github.com/pepkit/looper/zipball/master && \
    pip install pararead && \
    pip install pandas && \
    pip install piper && \
    pip install cutadapt

# Install R
RUN DEBIAN_FRONTEND=noninteractive apt-get --assume-yes install r-base r-base-dev && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('devtools')" && \
    Rscript -e "devtools::install_github('databio/peppro', subdir='PEPPROr')"

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
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && \
    unzip fastqc_v0.11.7.zip && \
    cd /home/tools/FastQC && \
    chmod 755 fastqc && \ 
    ln -s /home/tools/FastQC/fastqc /usr/bin/

# Install fqdedup
WORKDIR /home/tools/
RUN git clone https://github.com/guertinlab/fqdedup.git && \
    cd fqdedup && \
    cargo build --release && \
    ln -s /home/tools/fqdedup/target/release/fqdedup /usr/bin/

# Install fastx_toolkit
WORKDIR /home/tools/
RUN git clone https://github.com/agordon/libgtextutils.git && \
    cd libgtextutils && \
    ./reconf && \
    ./configure && \
    make && \
    make install && \
    cd /home/tools/ && \
    git clone https://github.com/agordon/fastx_toolkit && \
    cd fastx_toolkit && \
    ./reconf && \
    ./configure && \
    make && \
    make install

# Install seqOutBias
WORKDIR /home/tools/
RUN curl https://sh.rustup.rs -sSf | sh && \
    source $HOME/.cargo/env && \
    git clone https://github.com/guertinlab/seqOutBias.git && \
    cd seqOutBias && \
    cargo build --release && \
    ln -s /home/tools/seqOutBias/target/release/seqOutBias /usr/bin/

# Install pigz
WORKDIR /home/tools/
RUN wget https://zlib.net/pigz/pigz-2.4.tar.gz && \
    tar -zxvf pigz-2.4.tar.gz && \
    ln -s /home/tools/pigz-2.4/pigz /usr/bin/


# Set environment variables
ENV PATH=/home/tools/bin:/home/tools/:/home/src/bowtie2-2.3.5.1:/home/src/samtools-1.7:/home/src/htslib-1.7:$PATH \
    R_LIBS_USER=/usr/local/lib/R/site-library/ \
    PYTHONPATH=/usr/local/lib/python3.6/dist-packages:$PYTHONPATH

# Define default command
WORKDIR /home/
CMD ["/bin/bash"]

# Clean up APT when done.
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*