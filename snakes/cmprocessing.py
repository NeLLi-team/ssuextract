import pandas as pd
import sys
from collections import defaultdict
from collections import Counter
import glob
import re

# file with all fna merged
bindir = sys.argv[1]
# m8 blastn output
blastnm8 = sys.argv[2]
output = sys.argv[3]


def get_taxa(blastnm8):
 taxa_inm8 = []
 with open(blastnm8, "r") as infile:
  for line in infile:
   if line.split("|")[0] not in taxa_inm8:
    taxa_inm8.append(line.split("|")[0])
 return taxa_inm8


category_mapping = {
    "Bacteria": "BacteriaSSU",
    "Mitochondria": "MitochondriaSSU",
    "Chloroplast": "PlastidSSU",
    "Holosporales": "HolosporalesSSU",
    "Patescibacteria": "PatescibacteriaSSU",
    "Legionellales": "LegionellalesSSU",
    "Rickettsiales": "RickettsialesSSU",
    "Dependentiae": "DependentiaeSSU",
    "Cyanobacteria": "CyanobacteriaSSU",
    "Archaea": "ArchaeaSSU",
    "Eukaryota": "EukaryotaSSU",
}

def get_categories(blastnm8, taxon):
    seen = []
    results = []
    with open(blastnm8, "r") as infile:
        for line in infile:
            hitsplit = re.split("-|_|;", line)
            seenstring = line.split("|")[0] + "_" + line.split("|")[1]
            if line.startswith(taxon + "|") and seenstring not in seen:
                for hit in hitsplit:
                    if hit in category_mapping:
                        results.append(category_mapping[hit])
                seen.append(seenstring)
    return results


def get_taxanames(bindir):
    taxa = []
    for binf in glob.glob(bindir + "/*.fna"):
        taxa.append(binf.split("/")[-1].split(".")[0])
    return taxa


resultsdict = defaultdict(list)
for taxon in get_taxa(blastnm8):
    resultsdict[taxon] = Counter(get_categories(blastnm8, taxon))

df = pd.DataFrame.from_dict(resultsdict).T
df.fillna(0, inplace=True)
df = df.astype(int)
df.to_csv(output, sep='\t')
