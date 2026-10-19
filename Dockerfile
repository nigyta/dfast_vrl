# dfast_vrl image: VADR built from source with no model libraries, no conda,
# source COPYed in. Models are mounted from the host at $VADRMODELDIR.
# See docs/superpowers/specs/2026-06-03-docker-slim-design.md

# ---------------------------------------------------------------------------
# Stage 1: build VADR. The compilers, autoconf and source trees stay here.
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ubuntu:noble AS vadr-build

ARG VADR_VERSION="1.7.3"
ENV VADRINSTALLDIR=/opt/vadr

# Build deps listed in vadr-install.sh (C/C++ compiler, make, perl, curl, unzip,
# git) plus autoconf (Bio-Easel's easel) and zlib (minimap2). python3 is
# deliberately NOT installed: vadr-install.sh then skips R2DT (used only by
# `v-annotate.pl --draw_r2dt`) with a warning, which is what we want.
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      autoconf \
      ca-certificates \
      curl \
      git \
      libinline-c-perl \
      perl \
      unzip \
      wget \
      zlib1g-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# `--models none` (added in VADR 1.7.3) installs the software only. dfast_vrl
# pins its own model versions and mounts them at runtime (download_vadr_models.sh).
# The cleanup keeps only what the runtime stage needs: the infernal binaries
# (not its source tree), and no tarballs, tests or docs.
WORKDIR ${VADRINSTALLDIR}
RUN curl -fsSL https://raw.githubusercontent.com/ncbi/vadr/vadr-${VADR_VERSION}/vadr-install.sh -o vadr-install.sh && \
    bash vadr-install.sh linux --models none
# Separate layer so the cleanup can be tuned without recompiling (layers in this
# stage do not reach the final image). Bio-Easel keeps blib/ (compiled module)
# and scripts/ (v-annotate.pl runs scripts/esl-ssplit.pl).
RUN rm -rf vadr-install.sh *.zip *.tar.gz \
      R2DT r2dt-venv traveler jiffy-infernal-hmmer-scripts r2dt-install.log \
      vadr/testfiles vadr/documentation minimap2/test && \
    find infernal -mindepth 1 -maxdepth 1 ! -name binaries -exec rm -rf {} + && \
    find fasta -mindepth 1 -maxdepth 1 ! -name bin -exec rm -rf {} + && \
    find minimap2 -mindepth 1 -maxdepth 1 ! -name minimap2 -exec rm -rf {} + && \
    find Bio-Easel -mindepth 1 -maxdepth 1 ! -name blib ! -name scripts -exec rm -rf {} + && \
    find /opt -type f -iname "*.fa" -delete

# ---------------------------------------------------------------------------
# Stage 2: runtime.
# ---------------------------------------------------------------------------
FROM --platform=linux/amd64 ubuntu:noble

LABEL software="DFAST_VRL"
LABEL description="Viral genome annotation -> DDBJ MSS (slim image)"
LABEL website="https://github.com/nigyta/dfast_vrl"
LABEL maintainer="Yasuhiro Tanizawa"

# VADR environment, as printed at the end of vadr-install.sh. VADRMODELDIR points
# at the host mount instead of $VADRINSTALLDIR (no models are installed there).
ENV LC_ALL=C \
    VADRINSTALLDIR=/opt/vadr \
    VADRMODELDIR=/vadr_models
ENV VADRSCRIPTSDIR=$VADRINSTALLDIR/vadr \
    VADRINFERNALDIR=$VADRINSTALLDIR/infernal/binaries \
    VADREASELDIR=$VADRINSTALLDIR/infernal/binaries \
    VADRHMMERDIR=$VADRINSTALLDIR/infernal/binaries \
    VADRBIOEASELDIR=$VADRINSTALLDIR/Bio-Easel \
    VADRSEQUIPDIR=$VADRINSTALLDIR/sequip \
    VADRBLASTDIR=$VADRINSTALLDIR/ncbi-blast/bin \
    VADRFASTADIR=$VADRINSTALLDIR/fasta/bin \
    VADRMINIMAP2DIR=$VADRINSTALLDIR/minimap2
ENV VADRCONFIGFILE=$VADRSCRIPTSDIR/vadr.config \
    VADRMINISCRIPTSDIR=$VADRSCRIPTSDIR/miniscripts \
    PERL5LIB=$VADRSCRIPTSDIR:$VADRSEQUIPDIR:$VADRBIOEASELDIR/blib/lib:$VADRBIOEASELDIR/blib/arch \
    PATH=$VADRSCRIPTSDIR:$VADRSCRIPTSDIR/miniscripts:$VADRINFERNALDIR:$VADRSEQUIPDIR:$VADRBLASTDIR:$VADRFASTADIR:$VADRMINIMAP2DIR:$PATH

# Runtime deps:
#   perl, libinline-perl        -> VADR scripts / Bio-Easel. Inline only loads the
#                                  .so compiled in stage 1, so Inline::C (which
#                                  drags in gcc) is not needed here
#   libgomp1                    -> OpenMP runtime for the BLAST+ binaries (no gcc here
#                                  to pull it in)
#   liblwp-protocol-https-perl,
#   libmozilla-ca-perl          -> VADR's HTTPS calls (Mozilla::CA needed since 1.6.4)
#   wget / curl / unzip         -> download_vadr_models.sh, snpEff install
#   git                         -> pip-install dr_tools from GitHub
#   python3/pip/is-python3      -> scripts use `#!/usr/bin/env python` + Biopython
#   mafft                       -> MSA in dfv/detect_variants.py
#   openjdk-17-jre-headless     -> runs snpEff (Java)
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates \
      curl \
      git \
      libgomp1 \
      libinline-perl \
      liblwp-protocol-https-perl \
      libmozilla-ca-perl \
      mafft \
      openjdk-17-jre-headless \
      perl \
      procps \
      python-is-python3 \
      python3 \
      python3-pip \
      unzip \
      wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

COPY --from=vadr-build /opt/vadr /opt/vadr

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
ENV VADR_SCOV2_MODELS_VERSION="1.6.3-1" \
    VADR_MPXV_MODELS_VERSION="1.7-2" \
    VADR_RSV_MODELS_VERSION="1.7-1" \
    VADR_COX1_MODELS_VERSION="1.2-2" \
    VADR_CORONA_MODELS_VERSION="1.3-3" \
    VADR_FLAVI_MODELS_VERSION="1.7-2" \
    VADR_CALCI_MODELS_VERSION="1.7-2" \
    VADR_FLU_MODELS_VERSION="1.7-2" \
    VADR_HMPV_MODELS_VERSION="1.01"

# Copy the local source tree instead of cloning from GitHub: Docker invalidates
# this layer automatically whenever a copied file changes (no manual cache-bust
# knob), and the image reflects exactly this checkout. .dockerignore decides what
# lands here, so it must keep dfv/, refs/, and the three CLI scripts.
COPY . /dfast_vrl
RUN ln -s /dfast_vrl/dfast_vrl /usr/bin/dfast_vrl && \
    ln -s /dfast_vrl/vadr2mss.py /usr/bin/vadr2mss.py && \
    ln -s /dfast_vrl/cox1_to_ddbj.py /usr/bin/cox1_to_ddbj.py

WORKDIR /data

# Fail the build if VADR (incl. Bio-Easel's compiled code) or dfast_vrl is broken.
# Fail the build on any VADR binary with an unresolved shared library (stage 2
# lacks the toolchain stage 1 linked against). v-annotate.pl refuses to start,
# even for -h, if $VADRMODELDIR does not exist, so create the (empty) mount point.
RUN ! find /opt/vadr -type f \( -perm -u+x -o -name "*.so" \) -exec ldd {} \; 2>/dev/null | grep "not found" && \
    mkdir -p ${VADRMODELDIR} && \
    v-annotate.pl -h > /dev/null && \
    dfast_vrl --version

CMD ["/bin/bash"]
