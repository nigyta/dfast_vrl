# dfast_vrl image: prebuilt VADR base (no source compile), no conda, source COPYed in.
# See docs/superpowers/specs/2026-06-03-docker-slim-design.md
FROM --platform=linux/amd64 staphb/vadr:1.7-slim

LABEL software="DFAST_VRL"
LABEL description="Viral genome annotation -> DDBJ MSS (slim image)"
LABEL website="https://github.com/nigyta/dfast_vrl"
LABEL maintainer="Yasuhiro Tanizawa"

# dfast_vrl mounts VADR models from the host; override the base's VADRMODELDIR.
# All other VADR env vars / PATH are inherited from staphb/vadr:1.7-slim.
ENV VADRMODELDIR=/vadr_models \
    LC_ALL=C

# Runtime deps missing from the base image:
#   git                     -> pip-install dr_tools from GitHub
#   python3/pip/is-python3  -> scripts use `#!/usr/bin/env python` + Biopython
#   mafft                   -> MSA in dfv/detect_variants.py
#   openjdk-17-jre-headless -> runs snpEff (Java)
# (wget / unzip / curl already present in the base.)
RUN apt-get update && apt-get install -y --no-install-recommends \
      git \
      python3 \
      python3-pip \
      python-is-python3 \
      mafft \
      openjdk-17-jre-headless && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# snpEff 5.0 (pinned: committed snpEffectPredictor.bin DBs were built with 5.0).
# The core zip extracts to /opt/snpEff/snpEff.jar. Provide a `snpEff` command on
# PATH so dfv/detect_variants.py ("snpEff -c ... ") works unchanged.
RUN wget -q https://snpeff-public.s3.amazonaws.com/versions/snpEff_v5_0_core.zip -O /tmp/snpeff.zip && \
    unzip -q /tmp/snpeff.zip -d /opt && \
    rm /tmp/snpeff.zip && \
    test -f /opt/snpEff/snpEff.jar && \
    printf '#!/bin/sh\nexec java -jar /opt/snpEff/snpEff.jar "$@"\n' > /usr/local/bin/snpEff && \
    chmod +x /usr/local/bin/snpEff

# Python deps (no conda). pandas intentionally omitted (unused in code; pulled in
# transitively only if dr_tools requires it).
RUN python3 -m pip install --no-cache-dir --break-system-packages \
      biopython==1.84 \
      "git+https://github.com/ddbj/dr_tools.git"

# VADR model versions consumed at RUNTIME by the pipeline: dfv/vadr.py and
# dfv/vadr2mss_config.py shell-expand $VADR_*_MODELS_VERSION into the v-annotate.pl
# --mdir path, and dfv/reference_models.py reads VADR_SCOV2_MODELS_VERSION via
# os.getenv. Models themselves are NOT bundled (mounted at $VADRMODELDIR=/vadr_models).
# Placed before the COPY so the expensive apt/snpEff/pip layers above stay cached.
ENV VADR_SCOV2_MODELS_VERSION="1.3-2" \
    VADR_MPXV_MODELS_VERSION="1.4.2-1" \
    VADR_RSV_MODELS_VERSION="1.5-2" \
    VADR_COX1_MODELS_VERSION="1.2-2" \
    VADR_CORONA_MODELS_VERSION="1.3-3" \
    VADR_FLAVI_MODELS_VERSION="1.2-1" \
    VADR_CALCI_MODELS_VERSION="1.2-1"

# Copy the local source tree instead of cloning from GitHub: Docker invalidates
# this layer automatically whenever a copied file changes (no manual cache-bust
# knob), and the image reflects exactly this checkout. .dockerignore decides what
# lands here, so it must keep dfv/, refs/, and the three CLI scripts.
COPY . /dfast_vrl
RUN ln -s /dfast_vrl/dfast_vrl /usr/bin/dfast_vrl && \
    ln -s /dfast_vrl/vadr2mss.py /usr/bin/vadr2mss.py && \
    ln -s /dfast_vrl/cox1_to_ddbj.py /usr/bin/cox1_to_ddbj.py

WORKDIR /data

RUN dfast_vrl --version

CMD ["/bin/bash"]
