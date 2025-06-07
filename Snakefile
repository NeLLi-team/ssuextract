import subprocess
import glob
import os
import pathlib
from pathlib import Path

modeldir = Path(config["modeldir"]) # dir with models ending with .cm
querydir = Path(config["querydir"]) # dir with assemblies ending with .fna
cthreads = config["threads_per_job"]
FNABASE = [x.stem for x in querydir.iterdir() if x.is_file() and x.suffix in [".fa", ".fna", ".fasta"]]
MODELSBASE = [x.stem for x in modeldir.iterdir() if x.is_file() and x.suffix in [".cm"]]

# Create results directory name based on query directory
querydir_name = os.path.basename(str(querydir))
outdir = f"results/{querydir_name}/"

rule all:
   input:
      expand(outdir + "out/{fnaf}_{cmodel}.out", fnaf=FNABASE, cmodel=MODELSBASE),
      expand(outdir + "extracted/{fnaf}_{cmodel}.fna", fnaf=FNABASE, cmodel=MODELSBASE),
      expand(outdir + "m8/{fnaf}_{cmodel}.m8", fnaf=FNABASE, cmodel=MODELSBASE),
      outdir + "m8/merged.m8",
      outdir + "cmsearch_summary.tab",
      outdir + "cmsearch_summary.tsv"


rule check_bins_fna:
  # check header format in fnaf and edit if necessary
  conda:
    "config/environment.yml"
  input:
    str(querydir) + "/{fnaf}.fna"
  output:
    str(outdir) + "fna/{fnaf}.fna"
  shell:
    """
    python scripts/rename_fnaheaders.py {input} {output}
    """

rule run_cmsearch:
   # cmsearch for a set of covariance models
   conda:
      "config/environment.yml"
   input:
      str(outdir) + "fna/{fnaf}.fna",
      str(modeldir) + "/{cmodel}.cm"
   output:
      str(outdir) + "out/{fnaf}_{cmodel}.out"
   threads:
      cthreads
   shell:
      """
      cmsearch --anytrunc --cpu {threads} -o /dev/null --tblout {output} {input[1]} {input[0]}
      """

rule get_cmstats:
    # retrieve alignment stats to define coordinates to extract
   conda:
      "config/environment.yml"
   input:
      str(outdir) + "out/{fnaf}_{cmodel}.out"
   output:
      str(outdir) + "stats/{fnaf}_{cmodel}.seqmap",
      str(outdir) + "stats/{fnaf}_{cmodel}.seqsumt"
   shell:
      """
      python3 scripts/get_cmstats.py {input} {outdir}stats/{wildcards.fnaf}_{wildcards.cmodel} > {outdir}stats/{wildcards.fnaf}_{wildcards.cmodel}.log
      """


rule extract_cmhits:
   # extract hits based on positions in cmsearchout, also rc if necessary
   conda:
      "config/environment.yml"
   input:
      str(outdir) + "stats/{fnaf}_{cmodel}.seqmap",
      str(outdir) + "fna/{fnaf}.fna"
   output:
      str(outdir) + "extracted/{fnaf}_{cmodel}.fna"
   params:
      min_length = config.get("min_extract_length", 30)
   shell:
      """
      python3 scripts/get_cmsequences.py {input[1]} {input[0]} {output} {params.min_length}
      """


rule annotate_cmhits:
   conda:
      "config/environment.yml"
   input:
      str(outdir) + "extracted/{fnaf}_{cmodel}.fna"
   params:
      db = "resources/database/silva-138-1_pr2-4-12"
   output:
      outdir + "m8/{fnaf}_{cmodel}.m8"
   threads:
      cthreads
   shell:
      """
      blastn -outfmt 6 -db {params.db} -query {input} -max_target_seqs 5 -num_threads {threads} -out {output}
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
      "config/environment.yml"
   input:
      outdir + "m8/merged.m8"
   output:
      outdir + "cmsearch_summary.tab"
   shell:
      """
      python scripts/cmprocessing.py fna {input} {output}
      """

rule build_summary_tab:
   conda:
      "config/environment.yml"
   input:
      outdir + "cmsearch_summary.tab",
      expand(outdir + "extracted/{sample}_{model}.fna", sample=FNABASE, model=MODELSBASE),
      expand(outdir + "m8/{sample}_{model}.m8", sample=FNABASE, model=MODELSBASE),
      expand(outdir + "stats/{sample}_{model}.seqsumt", sample=FNABASE, model=MODELSBASE)
   output:
      outdir + "cmsearch_summary.tsv"
   run:
      import subprocess
      import shutil
      import os
      
      # Create a temporary symlink to mimic old structure
      temp_link = os.path.join(str(querydir), "cmsearch_out")
      if os.path.exists(temp_link):
          os.remove(temp_link)
      os.symlink(os.path.abspath(outdir), temp_link)
      
      try:
          # Run the script
          subprocess.run(["python", "scripts/get_table.py", str(querydir), str(modeldir)], check=True)
          
          # Move the output to the correct location
          shutil.move(os.path.join(outdir, "cmsearch_summary.tsv"), str(output[0]))
      finally:
          # Clean up symlink
          if os.path.exists(temp_link):
              os.remove(temp_link)
