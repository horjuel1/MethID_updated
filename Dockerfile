FROM rocker/r-ver:4.4.0

ENV DEBIAN_FRONTEND=noninteractive

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    tar \
    ca-certificates \
    pkg-config \
    cmake \
    samtools \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libgit2-dev \
    libudunits2-dev \
    libuv1-dev \
    && rm -rf /var/lib/apt/lists/*

# Install modkit v0.6.1
WORKDIR /tmp
RUN wget -O modkit.tar.gz https://github.com/nanoporetech/modkit/releases/download/v0.6.1/modkit_v0.6.1_u16_x86_64.tar.gz && \
    tar -xzf modkit.tar.gz && \
    mv dist_modkit_v0.6.1_481e3c9/modkit /usr/local/bin/modkit && \
    chmod +x /usr/local/bin/modkit && \
    rm -f modkit.tar.gz && \
    rm -rf dist_modkit_v0.6.1_481e3c9

WORKDIR /app

# Install R dependencies first for Docker layer caching
COPY bin/requirements.r /app/requirements.r
RUN Rscript /app/requirements.r

# Copy the rest of the pipeline
COPY . /app

# Validate tools and packages
RUN modkit --help >/dev/null && \
    samtools --version && \
    R -e "suppressWarnings(suppressMessages(library(GenomicRanges))); \
          suppressWarnings(suppressMessages(library(tidyverse))); \
          suppressWarnings(suppressMessages(library(data.table))); \
          suppressWarnings(suppressMessages(library(gridExtra))); \
          suppressWarnings(suppressMessages(library(limSolve))); \
          suppressWarnings(suppressMessages(library(argparser))); \
          suppressWarnings(suppressMessages(library(jsonlite))); \
          suppressWarnings(suppressMessages(library(R.utils))); \
          cat('All packages loaded successfully\n')"

CMD ["R"]
