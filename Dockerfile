# FINAL Dockerfile
FROM smedanno-base
WORKDIR /smedanno

# become root for system and conda writes 
USER root
SHELL ["/bin/bash", "-c"]

# optional apt bits that conda doesn't give you
RUN apt-get update && apt-get install -y --no-install-recommends perl \
 && rm -rf /var/lib/apt/lists/*

# add bioconda packages to the existing  env
RUN /opt/conda/bin/mamba install -y -n smedanno -c bioconda bedtools regtools

# revert to the normal conda shell
USER smedanno
SHELL ["/opt/conda/bin/conda", "run", "--live-stream", "-n", "smedanno", "/bin/bash", "-lc"]

# pipeline code and models
COPY --chmod=755 smedanno.sh functional_annotation.R tagXSstrandedData.awk \
                 deepsplice.py ./
COPY --chmod=755 DSmodels/model_*.pt /smedanno/DSmodels/

ENTRYPOINT ["conda","run","--live-stream","-n","smedanno","/smedanno/smedanno.sh"]
