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

## Docker/Singularity container
### Docker
  Available from [nigyta/dfast_vrl](https://hub.docker.com/r/nigyta/dfast_vrl/)  
  Includes VADR and Python modules  
  Modified from [staphb/vadr](https://hub.docker.com/r/staphb/vadr/)

Example:
Invoke `docker run` from the application root directory (where the dfast_vrl script exists).

```
docker run -v $PWD:/data -it --rm nigyta/dfast_vrl:latest dfast_vrl -i examples/LC570964-6.draft.contigs.fa -m examples/metadata.txt --enable_scaffolding --force -o OUTDOCKER
```

### Singularity

Sample data is available from github

- assembled contig FASTA
```
curl -O https://raw.githubusercontent.com/nigyta/dfast_vrl/main/examples/SRR10903401.meta_vrl.contig.fa
```

- metadata file
```
curl -O https://raw.githubusercontent.com/nigyta/dfast_vrl/main/examples/metadata.txt
```

Example:
```
singularity exec /path/to/dfast_vrl_0.0.2.sif dfast_vrl -i SRR10903401.meta_vrl.contig.fa -m metadata.txt -o dfv_result
```
