# DFAST_VRL
Viral genome annotation and data submission tool to DDBJ  
Mainly designed to process the SARS-Cov2 genome assembled using (META_VRL)[https://github.com/h-mori/meta_vrl]  

## Pipeline
1. preprocessing
2. annotation (using NCBI-VADR)
3. Format conversion (GenBank, DDBJ data submission format)

## Dependencies
[NCBI VADR](https://github.com/ncbi/vadr)

## Usage

## Docker container
  Includes VADR and Python modules  
  Modified from [staphb/vadr](https://hub.docker.com/r/staphb/vadr/)
