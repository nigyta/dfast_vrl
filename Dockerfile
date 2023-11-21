# base image
FROM ubuntu:22.04

# Modified from docker://staphb/vadr

# metadata - optional, but highly recommended
LABEL base.image="ubuntu:22.04"
LABEL dockerfile.version="1"
LABEL software="VADR"
LABEL software.version=${VADR_VERSION}
LABEL description="This software does viral annotations"
LABEL website="https://github.com/ncbi/vadr"
LABEL license="https://github.com/ncbi/vadr/blob/master/LICENSE"
LABEL maintainer="Yasuhiro Tanizawa"

# install dependencies via apt-get. Clean up apt garbage 
RUN apt-get update && apt-get install -y wget perl curl unzip build-essential git autoconf && \
 apt-get install -y libinline-c-perl liblwp-protocol-https-perl zlib1g-dev


ENV VADR_VERSION="1.6"\
  LC_ALL=C \
  VADRINSTALLDIR=/opt/vadr

ENV VADRSCRIPTSDIR=$VADRINSTALLDIR/vadr \
 VADRMODELDIR=/vadr_models \
 VADRINFERNALDIR=$VADRINSTALLDIR/infernal/binaries \
 VADREASELDIR=$VADRINSTALLDIR/infernal/binaries \
 VADRHMMERDIR=$VADRINSTALLDIR/hmmer/binaries \
 VADRBIOEASELDIR=$VADRINSTALLDIR/Bio-Easel \
 VADRSEQUIPDIR=$VADRINSTALLDIR/sequip \
 VADRFASTADIR=$VADRINSTALLDIR/fasta/bin \
 VADRBLASTDIR=$VADRINSTALLDIR/ncbi-blast/bin \
 VADRMINIMAP2DIR=$VADRINSTALLDIR/minimap2


ENV PERL5LIB=$VADRSCRIPTSDIR:$VADRSEQUIPDIR:$VADRBIOEASELDIR/blib/lib:$VADRBIOEASELDIR/blib/arch:$PERL5LIB \
 PATH=$VADRSCRIPTSDIR:$VADRBLASTDIR:$PATH


# install VADR
RUN mkdir -p ${VADRINSTALLDIR} && \
 cd ${VADRINSTALLDIR} &&\
 wget https://raw.githubusercontent.com/ncbi/vadr/release-${VADR_VERSION}/vadr-install.sh && \
 bash vadr-install.sh linux && \
 rm -r $VADRINSTALLDIR/vadr-models-calici $VADRINSTALLDIR/vadr-models-flavi/




RUN  wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh && \
  sh Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p /miniconda3 && \
  eval "$(/miniconda3/bin/conda shell.bash hook)" && \
  conda install -y -c bioconda mafft=7.475 snpeff=5.0 biopython=1.78 pandas && \
  rm Miniconda3-py38_4.9.2-Linux-x86_64.sh
ENV PATH=/miniconda3/bin:$PATH


# install VADR virus models
ENV VADR_SCOV2_MODELS_VERSION="1.3-2" \
  VADR_MPXV_MODELS_VERSION="1.4.2-1" \
  VADR_RSV_MODELS_VERSION="1.5-2" \
  VADR_COX1_MODELS_VERSION="1.2-1" \
  VADR_CORONA_MODELS_VERSION="1.3-3" \
  VADR_FLAVI_MODELS_VERSION="1.2-1" \
  VADR_CALCI_MODELS_VERSION="1.2-1"

# For DFAST_VRL, the reference data will be downloaded under VADRINSTALLDIR.
# RUN cd ${VADRINSTALLDIR} && \
#  wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/sarscov2/${VADR_SCOV2_MODELS_VERSION}/vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz && \
#  tar -xf vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz && \
#  rm -f vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz
# As of 2023.10.28, reference data have been removed from the container

# For VADR2MSS. Currently, the reference data are not included in the container.
# RUN cd ${VADRINSTALLDIR} && \
#  wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/mpxv/${VADR_MPXV_MODELS_VERSION}/vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz && \
#  tar -xf vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz && \
#  rm -f vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz

# RUN cd ${VADRINSTALLDIR} && \
#  wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/rsv/${VADR_RSV_MODELS_VERSION}/vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz && \
#  tar -xf vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz && \
#  rm -f vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz

# RUN cd ${VADRINSTALLDIR} && \
#  wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/cox1/${VADR_COX1_MODELS_VERSION}/vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz && \
#  tar -xf vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz && \
#  rm -f vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz

# RUN cd ${VADRINSTALLDIR} && \
#  wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/${VADR_CORONA_MODELS_VERSION}/vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz && \
#  tar -xf vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz && \
#  rm -f vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz



ARG INCREMENT_THIS_TO_DISABLE_CACHE_BELOW_THIS_LINE=2

RUN cd / && \
  git clone https://github.com/nigyta/dfast_vrl.git && \
  cd /usr/bin && \
  ln -s /dfast_vrl/dfast_vrl && \
  ln -s /dfast_vrl/vadr2mss.py


# set working directory
WORKDIR /data

RUN dfast_vrl --version 


