#!/bin/sh

if [ -z ${VADRMODELDIR} ]; then
  echo "Environmental variable VADRMODELDIR is not set. Data will be downloaded into the current directory."

  # For local environment
  VADR_SCOV2_MODELS_VERSION="1.6.3-1"
  VADR_MPXV_MODELS_VERSION="1.7-2"
  VADR_RSV_MODELS_VERSION="1.7-1"
  VADR_COX1_MODELS_VERSION="1.2-2"
  VADR_CORONA_MODELS_VERSION="1.3-3"
  VADR_FLAVI_MODELS_VERSION="1.7-2"
  VADR_CALCI_MODELS_VERSION="1.7-2"
  VADR_FLU_MODELS_VERSION="1.7-2"
  VADR_HMPV_MODELS_VERSION="1.01"
else
  # For environment within the docker container
  echo "VADR reference data will be downloaded to ${VADRMODELDIR}"
  cd ${VADRMODELDIR} 
fi



 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/sarscov2/${VADR_SCOV2_MODELS_VERSION}/vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-sarscov2-${VADR_SCOV2_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/mpxv/${VADR_MPXV_MODELS_VERSION}/vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-mpxv-${VADR_MPXV_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/rsv/${VADR_RSV_MODELS_VERSION}/vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-rsv-${VADR_RSV_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/cox1/${VADR_COX1_MODELS_VERSION}/vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-cox1-${VADR_COX1_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/coronaviridae/${VADR_CORONA_MODELS_VERSION}/vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-corona-${VADR_CORONA_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flaviviridae/${VADR_FLAVI_MODELS_VERSION}/vadr-models-flavi-${VADR_FLAVI_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-flavi-${VADR_FLAVI_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-flavi-${VADR_FLAVI_MODELS_VERSION}.tar.gz


 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/caliciviridae/${VADR_CALCI_MODELS_VERSION}/vadr-models-calici-${VADR_CALCI_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-calici-${VADR_CALCI_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-calici-${VADR_CALCI_MODELS_VERSION}.tar.gz

 wget https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/flu/${VADR_FLU_MODELS_VERSION}/vadr-models-flu-${VADR_FLU_MODELS_VERSION}.tar.gz && \
 tar -xf vadr-models-flu-${VADR_FLU_MODELS_VERSION}.tar.gz && \
 rm -f vadr-models-flu-${VADR_FLU_MODELS_VERSION}.tar.gz

 # hMPV is distributed by the Greninger Lab on GitHub, not on the NCBI FTP site.
 # The tag archive holds the model files under hmpv/; flatten them into the
 # usual vadr-models-<name>-<version>/ layout.
 wget -O vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION}.tar.gz https://github.com/greninger-lab/vadr-models-hmpv/archive/refs/tags/v${VADR_HMPV_MODELS_VERSION}.tar.gz && \
 mkdir -p vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION} && \
 tar -xzf vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION}.tar.gz --strip-components=2 -C vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION} vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION}/hmpv && \
 rm -f vadr-models-hmpv-${VADR_HMPV_MODELS_VERSION}.tar.gz
