import subprocess
import glob
import os
import pathlib

modeldir = Path(config["modeldir"]) # dir with models ending with .cm
querydir = Path(config["querydir"]) # dir with assemblies ending with .fna
FNABASE = [x.stem for x in querydir.iterdir() if x.is_file() and x.suffix in [".fa", ".fna", ".fasta"]]
MODELSBASE = [x.stem for x in modeldir.iterdir() if x.is_file() and x.suffix in [".cm"]]
outdir = str(querydir) + "/cmsearch_out/"

rule all:
   input:
      expand(outdir + "out/{fnaf}_{cmodel}.out", fnaf=FNABASE, cmodel=MODELSBASE),
      expand(outdir + "extracted/{fnaf}_{cmodel}.fna", fnaf=FNABASE, cmodel=MODELSBASE),
      expand(outdir + "m8/{fnaf}_{cmodel}.m8", fnaf=FNABASE, cmodel=MODELSBASE),
      outdir + "m8/merged.m8",
      outdir + "cmsearch_summary.tab"


rule check_bins_fna:
  # check header format in fnaf and edit if necessary
  conda:
    "snakes/cmsearchclass.yml"
  input:
    str(querydir) + "/{fnaf}.fna"
  output:
    str(outdir) + "fna/{fnaf}.fna"
  shell:
    """
    python snakes/rename_fnaheaders.py {input} {output}
    """


rule run_cmsearch:
   # cmsearch for a set of covariance models
   conda:
      "snakes/cmsearchclass.yml"
   input:
      str(outdir) + "fna/{fnaf}.fna",
      str(modeldir) + "/{cmodel}.cm"
   output:
      outdir + "out/{fnaf}_{cmodel}.out"
   shell:
      """
      cmsearch --cut_ga --cpu 8 --tblout {output} {input[1]} {input[0]}
      """


rule extract_cmhits:
   # extract hits based on positions in cmsearchout, also rc if necessary
   conda:
      "snakes/cmsearchclass.yml"
   input:
      str(outdir) + "out/{fnaf}_{cmodel}.out",
      str(outdir) + "fna/{fnaf}.fna"
   output:
      str(outdir) + "extracted/{fnaf}_{cmodel}.fna"
   shell:
      """
      python snakes/cmsearchout_extract_by_position_size.py {input[0]} {input[1]} {output} 800
      """


rule annotate_cmhits:
   conda:
      "snakes/cmsearchclass.yml"
   input:
      str(outdir) + "extracted/{fnaf}_{cmodel}.fna"
   params:
      db = "snakes/database/silva-138-1_pr2-4-12"
   output:
      outdir + "m8/{fnaf}_{cmodel}.m8"
   shell:
      """
      blastn -outfmt 6 -db {params.db} -query {input} -max_target_seqs 5 -num_threads 8 -out {output}
      """

rule merge_cmhits:
  input:
    expand(outdir + "m8/{fnaf}_{cmodel}.m8", fnaf=FNABASE, cmodel=MODELSBASE)
  output:
    outdir + "m8/merged.m8"
  params:
    outdir + "m8/"
  shell:
    """
    cat {params}*.m8 > {output}
    """

rule process_cmhits:
   conda:
      "snakes/cmsearchclass.yml"
   input:
      outdir + "m8/merged.m8"
   output:
      outdir + "cmsearch_summary.tab"
   shell:
      """
      python snakes/cmprocessing.py fna {input} {output}
      """
