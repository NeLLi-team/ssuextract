#!/usr/bin/env python3

import sys
import pandas as pd
import pathlib
import numpy as np

try:
    querydir = sys.argv[1]
    modeldir = sys.argv[2]

except:
    print("""
Miguel F. Romero, 2024
github.com/miferg
get_table.py retrieve information from a cmsearch analysis.
Usage:
get_table.py QUERY_DIR MODEL_DIR
    """)
    sys.exit()

    
# DEFINE FUNCTIONS

def load_extracted_names_lengths(sample, model, inpath):
    """extract sequence headers and length from fasta file"""
    seqlens_dict = {'name':[], 'length':[]}
    with open(inpath +'extracted/'+ sample +'_'+ model +'.fna') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                seqlens_dict['name'].append(line.strip('>'))
            else:
                seqlens_dict['length'].append(len(line))
    return seqlens_dict

def load_m8_best(sample, model, inpath):
    """load blast best hits"""
    m8 = pd.read_csv(inpath +'m8/'+ sample +'_'+ model +'.m8', sep='\t', header=None)
    m8.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    bestidx = []
    for qseqid in list(m8['qseqid'].unique()):
        ctab = m8.loc[m8['qseqid']==qseqid]
        besthit = ctab.sort_values(by='bitscore').index.values[-1]
        bestidx.append(besthit)
    bestm8 = m8.loc[bestidx]
    return bestm8

def get_blast_best_stats(sample, model, inpath, df, load_m8_best=load_m8_best):
    """add blast best hits to a data frame"""
    besthit_d = {'blast_sseqid':[], 'blast_pident':[], 'blast_length':[], 'blast_bitscore':[]}
    try: # if m8 is not empty
        bestm8 = load_m8_best(sample, model, inpath)

        for i in df.index.values:
            seqname = df['name'][i]
            if seqname in list(bestm8['qseqid'].unique()):
                ctab = bestm8.loc[bestm8['qseqid']==seqname]
                besthit_d['blast_sseqid'].append(ctab['sseqid'].item())
                besthit_d['blast_pident'].append(round(ctab['pident'].item(), 2))
                besthit_d['blast_length'].append(ctab['length'].item())
                besthit_d['blast_bitscore'].append(ctab['bitscore'].item())
            else:
                besthit_d['blast_sseqid'].append(np.nan)
                besthit_d['blast_pident'].append(np.nan)
                besthit_d['blast_length'].append(np.nan)
                besthit_d['blast_bitscore'].append(np.nan)
            
    except pd.errors.EmptyDataError: # if m8 is empty
        for i in df.index.values:
            besthit_d['blast_sseqid'].append(np.nan)
            besthit_d['blast_pident'].append(np.nan)
            besthit_d['blast_length'].append(np.nan)
            besthit_d['blast_bitscore'].append(np.nan)

    for key in ['blast_sseqid', 'blast_pident', 'blast_length', 'blast_bitscore']:
        df[key] = besthit_d[key]
        
    return df

def check_assembled(sample, model, inpath, df):
    """declare if sequence was extracted from fragmented hits"""
    try:
        seqsumt = pd.read_csv(inpath +'stats/'+ sample +'_'+ model +'.seqsumt', sep='\t', header=None)
        seqsumt.columns = ['contig', 'length']

        isassembled_l = []
        for i in df.index.values:
            isassembled = False
            contig = df['contig_name'][i]
            length = df['length'][i]
            if len(seqsumt.loc[(seqsumt['contig']==contig)&(seqsumt['length']==length)].index.values) > 0:
                isassembled = True
            isassembled_l.append(isassembled)

    except pd.errors.EmptyDataError: # if seqsumt is empty
        isassembled_l = [False for i in range(len(df.index.values))]

    df['is_assembled'] = isassembled_l
    return df


# MAIN PROGRAM

# read input and prepare iterators
cmpath = querydir +'/cmsearch_out/'
querydir = pathlib.Path(querydir)
modeldir = pathlib.Path(modeldir)
FNABASE = [x.stem for x in querydir.iterdir() if x.is_file() and x.suffix in [".fa", ".fna", ".fasta"]]
MODELSBASE = [x.stem for x in modeldir.iterdir() if x.is_file() and x.suffix in [".cm"]]

# initialize output table
outfile = pathlib.Path(cmpath +'cmsearch_summary.tsv')

# iterate over all samples and models!
for sample in FNABASE:
    for model in MODELSBASE:
        # load sequence information
        seqlens_dict = load_extracted_names_lengths(sample, model, cmpath)

        # create data fame and add information found in sequence headers
        df = pd.DataFrame(seqlens_dict['name'])
        
        if df.shape[0] > 0: # if there were extracted sequences
            df.columns = ['name']
            df['sample'] = [sample for i in range(df.shape[0])]
            df['model'] = [model for i in range(df.shape[0])]
            df['length'] = seqlens_dict['length']
            df['coordinates'] = [df['name'][i].split('|')[-3] for i in df.index.values]
            df['strand'] = [df['name'][i].split('|')[-2].split('_')[-1] for i in df.index.values]
            df['sequence_type'] = [df['name'][i].split('|')[-1] for i in df.index.values]
            df['contig_name'] = ['|'.join(df['name'][i].split('|')[:-3]) for i in df.index.values] # in case original sequence name has pipes

            # add blast best hit stats
            df = get_blast_best_stats(sample, model, cmpath, df)

            # declare if sequence was assembled
            df = check_assembled(sample, model, cmpath, df)

            # append rows to table
            df.to_csv(outfile, sep='\t', index=False, mode='a', header=not outfile.exists())
