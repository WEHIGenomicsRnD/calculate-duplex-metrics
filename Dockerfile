# syntax=docker/dockerfile:1
ARG R_VER=4.4.1
FROM rocker/r-ver:${R_VER}

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    python3 python3-pip python3-venv python3-distutils \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

# Hint to findpython (required by argparse)
ENV PYTHON=/usr/bin/python3
ENV PYTHON3=/usr/bin/python3

WORKDIR /app

# Install R package dependencies first (improves layer caching)
COPY DESCRIPTION .
RUN R -q -e "install.packages(c('remotes','BiocManager'), repos='https://cloud.r-project.org')" \
 && R -q -e "remotes::install_deps('.', dependencies=TRUE, repos=BiocManager::repositories(), upgrade='never')"

# Copy full source and install the package
COPY . .
RUN R CMD INSTALL . \
 && ln -s \
    "$(Rscript --vanilla -e 'cat(system.file("exec", "calc-duplex-metrics", package="CalcDuplexMetrics"))')" \
    /usr/local/bin/calc-duplex-metrics

CMD ["/bin/bash"]
