# --- Builder Stage for InterProScan ---
FROM ubuntu:24.04 AS builder

ARG IPS_VERSION=5.75-106.0
RUN apt-get update && apt-get install -y --no-install-recommends wget && rm -rf /var/lib/apt/lists/*
RUN wget --no-check-certificate -qO ipr.tgz https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPS_VERSION}/interproscan-${IPS_VERSION}-64-bit.tar.gz && \
    mkdir -p /opt/interproscan && \
    tar -xzf ipr.tgz -C /opt/interproscan --strip-components=1 && \
    rm ipr.tgz

# --- Final Stage: runtime image with conda-managed tools ---
# This is the actual base image that will be created.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install minimal base dependencies in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget ca-certificates bzip2 unzip libssl-dev openjdk-11-jre r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Miniconda
RUN wget -qO ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# 3. Copy the environment file
COPY environment.yml /tmp/environment.yml

# 4. Create the main conda environment.
RUN conda install -n base -c conda-forge --override-channels mamba -y && \
    mamba env create --override-channels -f /tmp/environment.yml && \
    conda clean -afy

# 5. Pre-build the smaller StringTie default environments
RUN conda create -y -n stringtie_short_211 -c bioconda -c conda-forge --override-channels "stringtie=2.1.1" && \
    conda create -y -n stringtie_mix_221 -c bioconda -c conda-forge --override-channels "stringtie=2.2.1" && \
    conda clean -afy

# Create a non-root user for improved security
RUN groupadd --gid 1001 smedanno && \
    useradd --uid 1001 --gid 1001 --shell /bin/bash --create-home smedanno

# Copy the InterProScan files from the 'builder' stage into this final stage
COPY --from=builder /opt/interproscan /opt/interproscan

# Switch to the non-root user
USER smedanno

# 7. Activate the main conda environment by default for all subsequent commands
SHELL ["/opt/conda/bin/conda", "run", "--live-stream", "-n", "smedanno", "/bin/bash", "-lc"]
