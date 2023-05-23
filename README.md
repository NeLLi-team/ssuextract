# SSUEXTRACT
* Identify, extract, annotate SSU rRNA genes
* Inputs are a directory with fna files and a directory with cm files (e.g. from RFAM)
* Annotation will only be successful for SSU rRNA genes but identification and extraction of hits using any other cm should work fine
* Database was too large to upload here, copy it from 
```
cp /global/cfs/cdirs/nelli/frederik/github/nelli/ssuextract_db/* snakes/database/
```

## How to run it
* Conda env that has snakemake installed
```
snakemake -j 16 --use-conda --config modeldir="snakes/models" querydir="example"
snakemake -j 16 --use-conda --config modeldir="{dir with cms}" querydir="{dir with fna files}"
```
* Results can be found in a new subdir "cmsearch_out" in the query directory
