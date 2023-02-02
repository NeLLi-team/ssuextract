snakemake -j 16 --use-conda --config modeldir="snakes/models" querydir="example"
snakemake -j 16 --use-conda --config modeldir="{dir with cms}" querydir="{dir with fna files}"
