# Slim Docker Image Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `Dockerfile.slim` (+ `Singularity.slim.def`, `.dockerignore`) that builds the dfast_vrl container from the prebuilt `staphb/vadr:1.7-slim` base without compiling VADR and without conda/mamba, cutting build time and image size while keeping behavior unchanged.

**Architecture:** Base on `staphb/vadr:1.7-slim` (VADR 1.7 prebuilt, ubuntu:noble). Inherit the base VADR env/PATH, override only `VADRMODELDIR=/vadr_models` (host-mounted models). Replace the conda env with apt (`mafft`, `openjdk-17-jre-headless`, `python3`, `git`) + pip (`biopython`, `dr_tools`). Install snpEff 5.0 from its official S3 zip with a `snpEff` wrapper. Drop unused `pandas`. Keep the existing `Dockerfile`/`Singularity.def` untouched.

**Tech Stack:** Docker (BuildKit), Singularity/Apptainer, staphb/vadr:1.7-slim, VADR 1.7, snpEff 5.0 (Java 17), MAFFT, Biopython, dr_tools.

Spec: `docs/superpowers/specs/2026-06-03-docker-slim-design.md`. Branch: `docker-slim`.

---

## Context

`nigyta/dfast_vrl`'s current `Dockerfile` compiles VADR from source (`vadr-install.sh`) and installs a full Miniforge/conda env with no cache cleanup. That makes builds slow and the image large (build-essential + conda pkg cache + bundled JRE all retained, no `.dockerignore`). We want a slimmer, faster-building image without changing what the pipeline does.

Code facts established during design (see spec):
- `pandas` is **not referenced anywhere** → safe to drop.
- `mafft` and `snpEff` are invoked as shell commands in `dfv/detect_variants.py` (SARS-CoV-2 variant-detection path only). snpEff DBs (`nigvrl`/`nigvrl2`) are committed under `refs/snpeff/data/` and were built with snpEff **5.0**, so snpEff must be pinned to 5.0.
- `biopython` is used in 10 modules; `dr_tools` is imported in `dfv/common.py:mss2json` (guarded by try/except).
- The CLI scripts use `#!/usr/bin/env python`, so a `python` (→python3) on PATH is required.
- `staphb/vadr:1.7-slim` ships `wget/curl/unzip/perl/build-essential` and all VADR env vars, **already deletes all bundled models**, but does **not** include `git`, `python3`, `mafft`, or any JRE.

---

## File Structure

- Create: `Dockerfile.slim` — slim image definition (the core deliverable).
- Create: `Singularity.slim.def` — Singularity mirror of the same recipe.
- Create: `.dockerignore` — shrink the build context (body is `git clone`d, not COPYd).
- Leave unchanged: `Dockerfile`, `Singularity.def`, CI workflow, all `dfv/` code.

No application code changes are needed.

---

### Task 1: Add `.dockerignore`

**Files:**
- Create: `.dockerignore`

- [ ] **Step 1: Write `.dockerignore`**

The image clones the source from GitHub inside the build, so nothing from the working tree is needed in the context. Exclude everything heavy.

```
# dfast_vrl source is git-cloned inside the build; context only needs the Dockerfile.
.git
.github
.devcontainer
.pytest_cache
.ipynb_checkpoints
.DS_Store
**/.DS_Store

# local outputs / scratch / test fixtures
OUT
OUT_*
*_test
cox1_test
dev
docs
examples
tests
refs

# editor / OS cruft
*.pyc
__pycache__
```

- [ ] **Step 2: Verify the context shrinks**

Run:
```bash
docker build -f Dockerfile.slim -t dfast_vrl:slim-ctxcheck --target=nonexistent . 2>&1 | head -3 || true
```
This will fail (no such target) — that is expected; the point is the first line shows the transferred context size. Compare with `du -sh .` to confirm large dirs (`OUT_*`, `examples`, `.git`) are excluded. Expected: transferred context is a few KB–MB, not hundreds of MB.

> Note: `.dockerignore` applies to whichever Dockerfile is built from this context, including the existing `Dockerfile`. That is fine — the existing `Dockerfile` also clones the source and does not `COPY` the tree.

- [ ] **Step 3: Commit**

```bash
git add .dockerignore
git commit -m "build: add .dockerignore to shrink Docker build context"
```

---

### Task 2: Write `Dockerfile.slim`

**Files:**
- Create: `Dockerfile.slim`

- [ ] **Step 1: Write `Dockerfile.slim`**

```dockerfile
# Slim build of dfast_vrl: prebuilt VADR base, no source compile, no conda.
# Companion to the original ./Dockerfile (kept for reference). See
# docs/superpowers/specs/2026-06-03-docker-slim-design.md
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
#   git                     -> clone dfast_vrl
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

# VADR model versions consumed at RUNTIME (shell-expanded into the v-annotate.pl
# --mdir path by dfv/vadr.py & dfv/vadr2mss_config.py; read via os.getenv in
# dfv/reference_models.py). Models are mounted at $VADRMODELDIR, not bundled.
ENV VADR_SCOV2_MODELS_VERSION="1.3-2" \
    VADR_MPXV_MODELS_VERSION="1.4.2-1" \
    VADR_RSV_MODELS_VERSION="1.5-2" \
    VADR_COX1_MODELS_VERSION="1.2-1" \
    VADR_CORONA_MODELS_VERSION="1.3-3" \
    VADR_FLAVI_MODELS_VERSION="1.2-1" \
    VADR_CALCI_MODELS_VERSION="1.2-1"

# Cache-busting knob mirrored from the original Dockerfile: bump to force a fresh clone.
ARG INCREMENT_THIS_TO_DISABLE_CACHE_BELOW_THIS_LINE=1

RUN git clone https://github.com/nigyta/dfast_vrl.git && \
    ln -s /dfast_vrl/dfast_vrl /usr/bin/dfast_vrl && \
    ln -s /dfast_vrl/vadr2mss.py /usr/bin/vadr2mss.py

WORKDIR /data

RUN dfast_vrl --version

CMD ["/bin/bash"]
```

- [ ] **Step 2: Build the image**

Run:
```bash
DOCKER_BUILDKIT=1 docker build -f Dockerfile.slim -t dfast_vrl:slim .
```
Expected: build completes; the final `RUN dfast_vrl --version` prints the version without error.

If the `test -f /opt/snpEff/snpEff.jar` step fails, the zip's top-level dir differs — run `docker run --rm dfast_vrl:slim ls /opt` (after a build that stops before that line) or inspect the zip with `unzip -l /tmp/snpeff.zip`, then fix the jar path in both the `test` and the wrapper `printf` lines.

- [ ] **Step 3: Smoke-test the binaries are wired up**

Run:
```bash
docker run --rm dfast_vrl:slim sh -c \
  'dfast_vrl --version && mafft --version 2>&1 | head -1 && java -version 2>&1 | head -1 && snpEff -version 2>&1 | head -3 && python -c "import Bio; print(\"biopython\", Bio.__version__)"'
```
Expected: dfast_vrl version line; MAFFT version (v7.x); a JRE 17 line; snpEff prints `SnpEff 5.0...`; `biopython 1.84`.

- [ ] **Step 4: Commit**

```bash
git add Dockerfile.slim
git commit -m "build: add Dockerfile.slim (staphb/vadr:1.7-slim base, no conda)"
```

---

### Task 3: End-to-end validation of the SARS-CoV-2 path (exercises mafft + snpEff)

This is the strongest functional check: only `dfast_vrl` (SARS-CoV-2) runs variant detection, so it validates VADR 1.7 + mounted models + mafft + snpEff 5.0 on Java 17 together.

**Files:** none (validation only).

**Prerequisite:** A host directory with the SARS-CoV-2 VADR model to mount at `/vadr_models`. If not already present, fetch it (sarscov2 model, version pinned in the repo):
```bash
mkdir -p /tmp/vadr_models && cd /tmp/vadr_models && \
wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/sarscov2/1.3-2/vadr-models-sarscov2-1.3-2.tar.gz && \
tar -xf vadr-models-sarscov2-1.3-2.tar.gz && rm -f vadr-models-sarscov2-1.3-2.tar.gz && cd -
```

- [ ] **Step 1: Run the SARS-CoV-2 pipeline against the slim image**

Run:
```bash
docker run --rm \
  -v "$PWD":/data \
  -v /tmp/vadr_models:/vadr_models \
  dfast_vrl:slim \
  dfast_vrl -i /data/examples/LC570964-6.draft.contigs.fa -m /data/examples/metadata.txt -o /data/OUT_slim --force
```
Expected: pipeline runs to completion (VADR → vadr2ddbj → MSS → mss2json). No `snpEff` Java errors and no MAFFT errors in the logs.

- [ ] **Step 2: Inspect outputs**

Run:
```bash
ls -1 OUT_slim/ && head -5 OUT_slim/*.annt.tsv 2>/dev/null
```
Expected: MSS files (`*.annt.tsv`, `*.seq.fa`) and the DFAST Record JSON exist and are non-empty.

> If snpEff fails to read `refs/snpeff/data/nigvrl/snpEffectPredictor.bin` with a version error, the JRE/snpEff combo is the cause: change `openjdk-17-jre-headless` to `openjdk-11-jre-headless` in `Dockerfile.slim`, rebuild (Task 2 Step 2), and re-run. Record the working combination in the spec.

- [ ] **Step 3: (Optional) generic VADR path**

If a model for another virus is available under `/vadr_models`, also run a quick `vadr2mss.py` check:
```bash
docker run --rm -v "$PWD":/data -v /tmp/vadr_models:/vadr_models dfast_vrl:slim \
  vadr2mss.py -i /data/examples/mpox/LC756923.fasta -m /data/examples/mpox/metadata_mpox.txt -o /data/OUT_slim_mpox -M mpox --force
```
Expected: completes and writes MSS files. (Skip if the mpox model is not present.)

- [ ] **Step 4: No commit** (validation only — record results in the PR/spec notes).

---

### Task 4: Measure size & build-time improvement

**Files:** none (measurement only).

- [ ] **Step 1: Build the original image for comparison**

Run:
```bash
DOCKER_BUILDKIT=1 docker build -f Dockerfile -t dfast_vrl:orig . 2>&1 | tail -2
```
(Use the committed original `Dockerfile`. This may take a while — it compiles VADR.)

- [ ] **Step 2: Compare image sizes**

Run:
```bash
docker images --format '{{.Repository}}:{{.Tag}}\t{{.Size}}' | grep dfast_vrl
```
Expected: `dfast_vrl:slim` is meaningfully smaller than `dfast_vrl:orig`. Record both numbers.

- [ ] **Step 3: (Optional) compare clean build times**

Run each with a clean cache and time it:
```bash
docker builder prune -af >/dev/null 2>&1
time DOCKER_BUILDKIT=1 docker build -f Dockerfile.slim -t dfast_vrl:slim --no-cache . >/dev/null
```
Expected: slim build is dramatically faster (no VADR compile). Record the wall-clock times.

- [ ] **Step 4: No commit** (numbers go into the spec/PR description).

---

### Task 5: Write `Singularity.slim.def` (mirror)

**Files:**
- Create: `Singularity.slim.def`

- [ ] **Step 1: Write `Singularity.slim.def`**

```
Bootstrap: docker
From: staphb/vadr:1.7-slim

%labels
    software DFAST_VRL
    description Viral genome annotation -> DDBJ MSS (slim image)
    maintainer Yasuhiro Tanizawa

%environment
    # VADR env vars / PATH are inherited from the docker base image.
    export LC_ALL=C
    export VADRMODELDIR=/vadr_models
    # model versions are shell-expanded into v-annotate.pl --mdir at runtime
    export VADR_SCOV2_MODELS_VERSION="1.3-2"
    export VADR_MPXV_MODELS_VERSION="1.4.2-1"
    export VADR_RSV_MODELS_VERSION="1.5-2"
    export VADR_COX1_MODELS_VERSION="1.2-1"
    export VADR_CORONA_MODELS_VERSION="1.3-3"
    export VADR_FLAVI_MODELS_VERSION="1.2-1"
    export VADR_CALCI_MODELS_VERSION="1.2-1"

%post
    export DEBIAN_FRONTEND=noninteractive

    # runtime deps missing from the base (git, python3, mafft, JRE for snpEff)
    apt-get update && apt-get install -y --no-install-recommends \
        git python3 python3-pip python-is-python3 mafft openjdk-17-jre-headless
    apt-get clean && rm -rf /var/lib/apt/lists/*

    # snpEff 5.0 (matches committed snpEffectPredictor.bin databases)
    wget -q https://snpeff-public.s3.amazonaws.com/versions/snpEff_v5_0_core.zip -O /tmp/snpeff.zip
    unzip -q /tmp/snpeff.zip -d /opt
    rm /tmp/snpeff.zip
    test -f /opt/snpEff/snpEff.jar
    printf '#!/bin/sh\nexec java -jar /opt/snpEff/snpEff.jar "$@"\n' > /usr/local/bin/snpEff
    chmod +x /usr/local/bin/snpEff

    # python deps (no conda; pandas intentionally omitted)
    python3 -m pip install --no-cache-dir --break-system-packages \
        biopython==1.84 "git+https://github.com/ddbj/dr_tools.git"

    # dfast_vrl
    cd /
    git clone https://github.com/nigyta/dfast_vrl.git
    ln -s /dfast_vrl/dfast_vrl /usr/bin/dfast_vrl
    ln -s /dfast_vrl/vadr2mss.py /usr/bin/vadr2mss.py

    # verify
    dfast_vrl --version

%runscript
    exec /bin/bash "$@"

%startscript
    exec /bin/bash "$@"
```

- [ ] **Step 2: (Optional) build if Singularity/Apptainer is available**

Run (only if `singularity`/`apptainer` is installed):
```bash
apptainer build --fakeroot dfast_vrl_slim.sif Singularity.slim.def 2>&1 | tail -20
```
Expected: build succeeds; `apptainer exec dfast_vrl_slim.sif dfast_vrl --version` prints the version. If no Singularity is available locally, skip — the recipe mirrors `Dockerfile.slim` line-for-line and is validated by parity.

- [ ] **Step 3: Commit**

```bash
git add Singularity.slim.def
git commit -m "build: add Singularity.slim.def mirroring Dockerfile.slim"
```

---

### Task 6: Finish the branch

**Files:** none.

- [ ] **Step 1: Update the spec with measured results**

Edit `docs/superpowers/specs/2026-06-03-docker-slim-design.md`: fill in the recorded image sizes / build times and the confirmed JRE+snpEff combination (Java 17 vs 11).

- [ ] **Step 2: Commit and use the finishing-a-development-branch skill**

```bash
git add docs/superpowers/specs/2026-06-03-docker-slim-design.md
git commit -m "docs: record slim image size/build-time results"
```
Then invoke `superpowers:finishing-a-development-branch` to decide merge/PR/cleanup (CI Dockerfile switch is out of scope — discuss separately).

---

## Verification (end-to-end)

1. `DOCKER_BUILDKIT=1 docker build -f Dockerfile.slim -t dfast_vrl:slim .` succeeds and prints the dfast_vrl version at the end (Task 2).
2. Smoke test shows mafft, Java 17, snpEff 5.0, and Biopython 1.84 all resolve inside the container (Task 2 Step 3).
3. The SARS-CoV-2 `dfast_vrl` pipeline runs to completion against mounted models and produces non-empty MSS + JSON outputs — exercising mafft + snpEff together (Task 3).
4. `docker images` shows `dfast_vrl:slim` is smaller than `dfast_vrl:orig`, and a `--no-cache` slim build is far faster than the original (Task 4).
5. `Singularity.slim.def` mirrors the Docker recipe (and builds if Apptainer is available) (Task 5).

## Risks / fallbacks

- **VADR 1.6.4 → 1.7 model compatibility:** the repo's pinned model versions (SCOV2 1.3-2 etc.) must load under VADR 1.7. If `dfast_vrl`/`vadr2mss.py` errors on `.minfo`, bump the affected model version when downloading for the test and note it. (Task 3 surfaces this.)
- **snpEff 5.0 on Java 17:** if the committed `.bin` fails to load, switch to `openjdk-11-jre-headless` and rebuild (Task 3 Step 2 note).
- **snpEff zip layout:** the build's `test -f /opt/snpEff/snpEff.jar` fails fast if the extracted path differs; fix the jar path in the `test`, wrapper, and Singularity recipe together.
