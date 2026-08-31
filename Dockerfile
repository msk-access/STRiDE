# syntax=docker/dockerfile:1
#
# Dockerfile for STRiDE — MSI prediction pipeline for MSK-ACCESS
#
# Build: docker build -t stride .
# Run:   docker run stride --help
#
# Nextflow requirements:
#   - bash at /bin/bash (for script execution)
#   - procps (ps command for task resource metrics)
#   - entrypoint overridden via --entrypoint="" in Nextflow Docker profile

FROM python:3.11-slim-bookworm

# OCI Labels for GitHub Container Registry
LABEL org.opencontainers.image.title="STRiDE"
LABEL org.opencontainers.image.description="MSI prediction pipeline for MSK-ACCESS"
LABEL org.opencontainers.image.url="https://github.com/msk-access/STRiDE"
LABEL org.opencontainers.image.source="https://github.com/msk-access/STRiDE"
LABEL org.opencontainers.image.vendor="Memorial Sloan Kettering Cancer Center"
LABEL org.opencontainers.image.licenses="AGPL-3.0"
LABEL org.opencontainers.image.authors="Karmelina Charlambous <charlk@mskcc.org>, Ronak Shah <shahr2@mskcc.org>"

# Version set dynamically via --build-arg in CI
ARG VERSION=dev
LABEL org.opencontainers.image.version="${VERSION}"

# procps: required by Nextflow for task resource metrics (ps command)
# bash:   required by Nextflow for script execution (/bin/bash)
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps bash && rm -rf /var/lib/apt/lists/*

# Install STRiDE with QC extras (plotly for interactive HTML reports)
COPY . /build
RUN pip install --no-cache-dir "/build[qc]" && rm -rf /build

# Verify installation
RUN stride --version

ENTRYPOINT ["stride"]
CMD ["--help"]
