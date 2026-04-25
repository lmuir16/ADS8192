# corinet — R package for gene co-expression networks from bulk RNA-seq
# Base: Ubuntu 24.04 LTS

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=UTC

# Minimal OS deps: fetch Miniconda only (R and libs come from conda)
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    wget \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_DIR=/opt/conda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p $CONDA_DIR \
    && rm /tmp/miniconda.sh
ENV PATH=$CONDA_DIR/bin:$PATH

# Use only conda-forge + bioconda (no defaults) to avoid Anaconda pkgs/* TOS prompts in Docker.
RUN conda config --system --remove channels defaults 2>/dev/null || true \
    && conda config --system --prepend channels conda-forge \
    && conda config --system --append channels bioconda \
    && conda config --system --set channel_priority strict

# R + Bioconductor + vignette/test stack from environment.yml
COPY environment.yml /tmp/environment.yml
RUN sed -i '/^prefix:/d' /tmp/environment.yml \
    && conda env create -y -p /opt/conda/envs/ads8192 -f /tmp/environment.yml \
    && rm /tmp/environment.yml \
    && conda clean -afy

ENV PATH=/opt/conda/envs/ads8192/bin:$PATH \
    R_HOME=/opt/conda/envs/ads8192/lib/R

WORKDIR /ads8192

# Package sources (exclude dev-only paths via .dockerignore when using COPY .)
COPY . .

# Rapp CLIs to ~/.local/bin; explicit CRAN mirror avoids interactive mirror selection
RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); install.packages("Rapp"); library(Rapp); Rapp::install_pkg_cli_apps("Rapp");'

RUN R CMD INSTALL /ads8192

RUN Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); Rapp::install_pkg_cli_apps("corinet")'

ENV PATH=/root/.local/bin:$PATH

# Default: run package tests (WORKDIR is package root)
CMD ["R", "--no-save", "-q", "-e", "devtools::test()"]
