# =========================================================
# Reproducible environment for the NGS variant calling pipeline.
# Built on the bioconda/miniforge base so all tool versions are
# pinned via environment.yml rather than "whatever apt has today".
# =========================================================
FROM condaforge/miniforge3:24.3.0-0

LABEL maintainer="Agnidipa Sett"
LABEL description="WGS/WES variant calling & clinical annotation pipeline (GIAB HG002, chr20)"

WORKDIR /workspace

# Docker CLI is needed inside the container so the deepvariant rule
# can shell out to `docker run google/deepvariant:...` (docker-outside-of-docker
# pattern). If you'd rather run everything conda-native, swap this rule
# for the bioconda `deepvariant` package instead -- see README notes.
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
        docker.io \
        && rm -rf /var/lib/apt/lists/*

COPY environment.yml /workspace/environment.yml
RUN conda env create -f /workspace/environment.yml && conda clean -afy

# Activate the pipeline env by default for interactive shells and RUN steps
SHELL ["conda", "run", "-n", "ngs-pipeline", "/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/ngs-pipeline/bin:$PATH

COPY . /workspace

ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "ngs-pipeline"]
CMD ["snakemake", "--cores", "4", "-p"]
