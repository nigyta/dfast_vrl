# DFAST_VRL
DFAST_VRL is a viral genome annotation and data submission tool to DDBJ.  
It internally uses NCBI-VADR and can generate data submission files to DDBJ.

## Pipeline
1. Preprocessing  
  Query genome is mapped to the reference to calculate the coverage against the reference. Runs of leading/trailing Ns are trimmed.   
  (Optional) Unmapped regions are removed, and mapped regions are scaffoled (ordered and concatenated with runs of Ns).  
2. Annotation  
  DFAST_VRL internally uses NCBI-VADR to annotate bilogical features (CDS, mature peptide, ribosomal slippage, stem loop).  
3. Variant detection  
  Based on the pairwise alignment with the reference genome using MAFFT
4. Output  
  GenBank format, VCF format, DDBJ data submission format  

## Dependencies
- [NCBI VADR](https://github.com/ncbi/vadr)
- NCBI-BLAST+ (included in the VADR package)
- snpEff
- MAFFT
- Biopython

Since DFAST_VRL depends on many external binaries and modules, it is recommended to use Docker/Singularity container as described below.

## Usage
```
dfast_vrl -i input_genome.fasta [-o output_directory] [-m metadata_file.txt]
```

Example:
```
./dfast_vrl -i examples/LC570964-6.draft.contigs.fa -m examples/metadata.txt
```

The input file must be a single-FASTA file consisting of one complete genome or a multi-FASTA file consisting of multiple sequences from a draft genome assebmly.

## Docker/Singularity container
### Docker  
  Available from [nigyta/dfast_vrl](https://hub.docker.com/r/nigyta/dfast_vrl/)  
  Includes NCBI VADR and other dependencies. Modified from [staphb/vadr](https://hub.docker.com/r/staphb/vadr/).

  
Example:  
In the example below, `docker run` is invoked from the application root directory (where the dfast_vrl script exists).

```
docker run -v $PWD:/data -it --rm nigyta/dfast_vrl:latest dfast_vrl -i /data/examples/LC570964-6.draft.contigs.fa -m /data/examples/metadata.txt -o /data/OUT
```

### Singularity

In the NIG-Super computer system, Singularity container is available at `/lustre6/public/vrl/dfast_vrl:latest.sif`.  
Sample data is available from github.

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
singularity exec /lustre6/public/vrl/dfast_vrl:latest.sif dfast_vrl -i SRR10903401.meta_vrl.contig.fa -m metadata.txt -o dfv_result
```
