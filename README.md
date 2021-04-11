# DFAST_VRL
Viral genome annotation and data submission tool to DDBJ  
Mainly designed to process a SARS-Cov2 genome assembled using [META_VRL])https://github.com/h-mori/meta_vrl)

## Pipeline
1. Preprocessing
2. Annotation (using NCBI-VADR)
3. Output (GenBank format, DDBJ data submission format)

## Dependencies
Biopython
[NCBI VADR](https://github.com/ncbi/vadr)

## Usage
```
dfast_vrl -i input_genome.fasta [-o output_directory] [-m metadata_file.txt] [--enable_scaffolding]
```
  
Example:
```
./dfast_vrl -i examples/LC570964-6.draft.contigs.fa -m examples/metadata.txt --enable_scaffolding --force
```  

## Docker container
  Includes VADR and Python modules  
  Modified from [staphb/vadr](https://hub.docker.com/r/staphb/vadr/)
