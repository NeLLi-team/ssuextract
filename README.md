# SSUEXTRACT
* Identify, extract, annotate SSU rRNA genes
* Annotation will only be successful for SSU rRNA genes but identification and extraction of hits using any other cm should work fine

## How to run it
* Run it from the repo base dir
* Inputs are a directory with fna files and a directory with cm files (e.g. from RFAM)
* Database was too large to upload here, copy it from 
```
cd ssuextract
wget -r -np -nH --cut-dirs=3 -P "snakes/database/" "https://portal.nersc.gov/cfs/nelli/ssuextract_db/"
```
* Conda env that has snakemake installed
```
snakemake -j 16 --use-conda --config modeldir="snakes/models" querydir="example" threads_per_job=4
snakemake -j 16 --use-conda --config modeldir="{dir with cms}" querydir="{dir with fna files}" threads_per_job=4
```
* Results can be found in a new subdir "cmsearch_out" in the query directory
