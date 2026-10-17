# PhageTermVirome 4.3 (Pasteur) — no public biocontainer exists, so the pipeline builds this image locally.
# Build:  docker build -t phage-pipeline/phagetermvirome:4.3 pipeline/containers/phageterm
FROM python:3.11-slim
# the package declares a dependency on the deprecated "sklearn" PyPI stub; this variable lets pip install it
ENV SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True \
    MPLBACKEND=Agg \
    PIP_NO_CACHE_DIR=1 \
    PYTHONPATH=/usr/local/lib/python3.11/site-packages/phagetermvirome
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/* \
 && pip install "phagetermvirome==4.3"
# the package imports its sibling modules without the package prefix, hence PYTHONPATH above; entry point is "phageterm"
WORKDIR /data
