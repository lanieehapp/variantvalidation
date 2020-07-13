# Base Image
FROM r-base:3.5.2

# Maintainer
MAINTAINER DaveLab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    curl \
    vim \
    devscripts \
    bcftools \ 
    less \
    wget \
    unzip \
    cmake \
    python3 \
    gawk \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    git \
    autoconf \
    bsdmainutils \
    bedtools \
    gcc-8-base \
    libmpx2 \
    libgcc-8-dev \
    bedops \
    tabix


WORKDIR /usr/local/bin

# install R required dependencies
RUN R --vanilla -e 'install.packages(c("vcfR", "stringr", "plyr"), repos="http://cran.us.r-project.org")'

# clone variantvalidation
ADD https://api.github.com/repos/lanieehapp/variantvalidation/git/refs/heads/ version.json
RUN git clone https://github.com/lanieehapp/variantvalidation.git

# add variantvalidation repo to SYSPATH
ENV PATH variantvalidation:$PATH

# change the permission of the repo
RUN chmod 777 -R variantvalidation
WORKDIR /usr/local/bin/variantvalidation



