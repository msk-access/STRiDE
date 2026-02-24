# syntax=docker/dockerfile:1
FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.title="STRiDE"
LABEL org.opencontainers.image.description="MSI prediction pipeline for MSK-ACCESS"
LABEL org.opencontainers.image.url="https://github.com/msk-access/STRiDE"
LABEL org.opencontainers.image.source="https://github.com/msk-access/STRiDE"
LABEL org.opencontainers.image.vendor="Memorial Sloan Kettering Cancer Center"
LABEL org.opencontainers.image.licenses="MIT"

# Version and creation timestamp set dynamically via --build-arg in CI
ARG VERSION=dev
LABEL org.opencontainers.image.version="${VERSION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps bash && rm -rf /var/lib/apt/lists/*

COPY . /build
RUN pip install --no-cache-dir /build && rm -rf /build

ENTRYPOINT ["stride"]
CMD ["--help"]
