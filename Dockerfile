# UMAP / LocusZoom / Z-score toolkit: R plotting + Python extract_z_lz.py
#
# Build (use amd64 if you are on Apple Silicon and need x86_64 plink2):
#   docker build --platform linux/amd64 -t umap-lz-plot .
#
# Run R plot (mount your data under /work):
#   docker run --rm -v "$PWD":/work -w /work umap-lz-plot \
#     Rscript plot_umap_simplified_multimodules_one_side.R \
#       --lz_files /work/data/lz.csv --out /work/out/fig --png
#
# Run Python extract (example):
#   docker run --rm -v "$PWD":/work -w /work umap-lz-plot \
#     python3 extract_z_lz.py --qtl_module_adata /work/a.h5ad ... --out /work/out/base
#
# Interactive shell:
#   docker run --rm -it -v "$PWD":/work -w /work umap-lz-plot bash

FROM rocker/tidyverse:4.4.3

ARG PLINK2_VER=v2.0.0-a.7.1

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/venv/bin:${PATH}"
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    wget \
    unzip \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# PLINK 2 (optional but recommended for extract_z_lz.py --safeld).
RUN set -eux; \
    wget -qO /tmp/plink2.zip \
      "https://github.com/chrchang/plink-ng/releases/download/${PLINK2_VER}/plink2_linux_x86_64.zip"; \
    unzip -oq /tmp/plink2.zip -d /tmp/plink2_install; \
    bin="$(find /tmp/plink2_install -type f \( -name plink2 -o -name plink2.exe \) | head -n 1)"; \
    test -n "$bin"; \
    install -m 755 "$bin" /usr/local/bin/plink2; \
    rm -rf /tmp/plink2.zip /tmp/plink2_install; \
    plink2 --version

# R packages used by plot_umap_simplified_multimodules_one_side.R (tidyverse image
# already provides ggplot2, dplyr, readr, stringr, etc.)
RUN install2.r --error --deps TRUE \
    data.table \
    optparse \
    patchwork \
    ggrepel \
    ragg \
    scattermore

WORKDIR /app

COPY requirements.txt /app/requirements.txt
RUN python3 -m venv /opt/venv \
    && /opt/venv/bin/pip install --upgrade pip \
    && /opt/venv/bin/pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

# Default: shell so you can run Rscript or python3 with your own args
CMD ["/bin/bash"]
