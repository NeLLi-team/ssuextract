import pandas as pd
import sys
import os
from collections import defaultdict
from collections import Counter
import glob
import re

# Constants
CATEGORY_MAPPING = {
    "Bacteria": "BacteriaSSU",
    "Mitochondria": "MitochondriaSSU", 
    "Chloroplast": "PlastidSSU",
    "Holosporales": "HolosporalesSSU",
    "Patescibacteria": "PatescibacteriaSSU",
    "Legionellales": "LegionellalesSSU",
}

def validate_inputs():
    """Validate command line arguments and input files."""
    if len(sys.argv) != 4:
        print("Usage: cmprocessing.py bindir blastnm8 output")
        sys.exit(1)
    
    bindir = sys.argv[1]
    blastnm8 = sys.argv[2] 
    output = sys.argv[3]
    
    # Validate input file exists
    if not os.path.exists(blastnm8):
        print(f"Error: Input file not found: {blastnm8}")
        sys.exit(1)
        
    return bindir, blastnm8, output

def get_taxa(blastnm8):
    """Extract unique taxa from blastn m8 output using set for efficiency."""
    taxa_inm8 = set()
    try:
        with open(blastnm8, "r") as infile:
            for line in infile:
                taxa_inm8.add(line.split("|")[0])
    except (IOError, IndexError) as e:
        print(f"Error processing taxa from {blastnm8}: {e}")
        sys.exit(1)
    return list(taxa_inm8)


# Update the constant with complete mapping
CATEGORY_MAPPING.update({
    "Rickettsiales": "RickettsialesSSU",
    "Dependentiae": "DependentiaeSSU", 
    "Cyanobacteria": "CyanobacteriaSSU",
    "Archaea": "ArchaeaSSU",
    "Eukaryota": "EukaryotaSSU",
})

def get_categories(blastnm8, taxon):
    """Extract categories for a given taxon using set for efficiency."""
    seen = set()
    results = []
    try:
        with open(blastnm8, "r") as infile:
            for line in infile:
                hitsplit = re.split("-|_|;", line)
                seenstring = line.split("|")[0] + "_" + line.split("|")[1]
                if line.startswith(taxon + "|") and seenstring not in seen:
                    for hit in hitsplit:
                        if hit in CATEGORY_MAPPING:
                            results.append(CATEGORY_MAPPING[hit])
                    seen.add(seenstring)
    except (IOError, IndexError) as e:
        print(f"Error processing categories for {taxon}: {e}")
        return []
    return results


def get_taxanames(bindir):
    taxa = []
    for binf in glob.glob(bindir + "/*.fna"):
        taxa.append(binf.split("/")[-1].split(".")[0])
    return taxa


def main():
    """Main function to process blastn output and generate summary table."""
    bindir, blastnm8, output = validate_inputs()
    
    # Build results dictionary first
    resultsdict = defaultdict(lambda: Counter())
    for taxon in get_taxa(blastnm8):
        resultsdict[taxon] = Counter(get_categories(blastnm8, taxon))
    
    # Convert to DataFrame
    df = pd.DataFrame.from_dict(resultsdict, orient='index')
    df.fillna(0, inplace=True)
    df = df.astype(int)
    df.to_csv(output, sep='\t')

if __name__ == "__main__":
    main()
