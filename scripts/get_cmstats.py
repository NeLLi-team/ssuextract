#!/usr/bin/env python3

import sys
import pandas as pd
import os
from pathlib import Path

def validate_inputs():
    """Validate command line arguments and input files."""
    if len(sys.argv) != 3:
        print("""
Retrieve information from a cmsearch output table.
Usage:
get_cmstats.py cmtable prefix
        """)
        sys.exit(1)
    
    cmtable = sys.argv[1]
    outpath = sys.argv[2]
    
    # Validate input file exists
    if not os.path.exists(cmtable):
        print(f"Error: Input file not found: {cmtable}")
        sys.exit(1)
    
    # Validate output directory exists
    output_dir = os.path.dirname(outpath)
    if output_dir and not os.path.exists(output_dir):
        print(f"Error: Output directory not found: {output_dir}")
        sys.exit(1)
    
    return cmtable, outpath

try:
    cmtable, outpath = validate_inputs()
except Exception as e:
    print(f"Error validating inputs: {e}")
    sys.exit(1)

    
## Function definitions

def get_clen(cmtable):
    """Retrieve the length of the covariance model"""
    cmfilename = cmtable.split('_')[-1].split('.')[0] +'.cm'
    with open('resources/models/'+ cmfilename) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('CLEN'):
                clen = int(line.split(' ')[-1])
                return clen

def get_singles(ctab, clen):
    """Find contigs with a a single alignment to the model or full alignments in a contig with multiple alignments"""
    singles = ctab['subject'].value_counts().loc[ctab['subject'].value_counts()==1].index.values # contigs with single alignments
    print('Single alignments')
    print(singles)
    print(len(list(ctab.loc[(ctab['subject'].isin(singles))&(ctab['include']=='!')].index.values)))
    singles_indexes = list(ctab.loc[(ctab['subject'].isin(singles))&(ctab['include']=='!')].index.values) # indexes with single alignments
    multi = ctab['subject'].value_counts().loc[ctab['subject'].value_counts()>1].index.values # multi alignments
    if len(multi) > 0: # if there are contigs with multiple hits
        print('Full alignments')
        print(list(ctab.loc[(ctab['subject'].isin(multi)) & (ctab['mdl-from']==1) & (ctab['mdl-to']==clen), 'subject'].unique()))
        print(len(list(ctab.loc[(ctab['subject'].isin(multi)) & (ctab['mdl-from']==1) & (ctab['mdl-to']==clen)].index.values)))
        singles_indexes += list(ctab.loc[(ctab['subject'].isin(multi)) & (ctab['mdl-from']==1) & (ctab['mdl-to']==clen)].index.values) # indexes with single + full alignments
    
    return singles_indexes

def get_sing_mdls_seqs(singles_indexes, outfmap, ctab):
    """Write the coordinates of single and full hits"""
    # outfile1 = open(outpath +'.singles.mdls','w') # length of each hit according to the model
    # outfile2 = open(outpath +'.singles.seqs','w') # length of each hit in the query sequence
    for i in singles_indexes:
        mdl = [ctab['mdl-from'][i], ctab['mdl-to'][i]]
        # outfile1.write(subject +'\t'+ str((max(mdl) - min(mdl))) +'\n')
        seq = [ctab['seq-from'][i], ctab['seq-to'][i]]
        # outfile2.write(subject +'\t'+ str((max(seq)) - min(seq)) +'\n')
        outfmap.write(ctab['subject'][i] +'\t'+ str(ctab['seq-from'][i]) +'\t'+ str(ctab['seq-to'][i]) +'\t'+ ctab['strand'][i] +'\tsimple\n')
    # outfile1.close()
    # outfile2.close()
    return outfmap
    
## GET SIZES FOR CONTIGS WITH MULTIPLE HITS

def get_repeated(ctab, singles_indexes):
    """Retrieve the alignment indexes for contigs with multiple hits excluding full alignments"""
    mctab = ctab.drop(singles_indexes, axis=0)
    repeated_raw = list(mctab.index.values)
    print('Multiple truncated alignments')
    print(list(mctab['subject'].unique()))
    print(len(repeated_raw))
    return repeated_raw

## GET COORDINATES 

def get_coords(repeated_raw, ctab):
    outfile = open(outpath +'.coords', 'w')
    
    for currentseq in list(ctab.loc[repeated_raw, 'subject'].unique()):
        outfile.write(currentseq +'\n')
        cmtab = ctab.loc[ctab['subject']==currentseq]

        index = 1
        aldictm = {}
        aldicts = {}

        for i in cmtab.index.values:
            strand = cmtab['strand'][i]
            mdl = [cmtab['mdl-from'][i], cmtab['mdl-to'][i]]
            if strand == '+':
                seq = [cmtab['seq-from'][i], cmtab['seq-to'][i]]
            else:
                seq = [-cmtab['seq-from'][i], -cmtab['seq-to'][i]]
            aldictm[mdl[0]] = '('*index + str(mdl[0])
            aldictm[mdl[1]] = str(mdl[1]) + ')'*index
            aldicts[seq[0]] = '('*index + str(seq[0])
            aldicts[seq[1]] = str(seq[1]) + ')'*index
            index += 1

        mdlal = list(aldictm.keys())
        mdlal.sort()

        mdlals = ''

        for i in mdlal:
            mdlals += '..'
            mdlals += aldictm[i]
            mdlals += '..'
    
        outfile.write('Model:    '+ mdlals +'\n')

        seqal = list(aldicts.keys())
        seqal.sort()

        seqals = ''

        for i in seqal:
            seqals += '..'
            seqals += aldicts[i]
            seqals += '..'
    
        outfile.write('Sequence: '+ seqals +'\n')
        outfile.write('\n')
    
    outfile.close()
    return 0

## GET PUTATIVE DUPLICATIONS/PARALOGS

def get_paralg_seqsumt_ins(repeated_raw, outfmap, ctab):
    outfile1 = open(outpath +'.paralg', 'w')    # sequences with putative duplications
    outfile2 = open(outpath +'.seqsumt', 'w')   # total lengths of composie hits including insertions
    outfile3 = open(outpath +'.insert', 'w')    # lenghts of each insertion between hit pairs in each contig

    for currentseq in list(ctab.loc[repeated_raw, 'subject'].unique()): 
        cmtab = ctab.loc[ctab['subject']==currentseq]

        index = 1
        index2 = 1
        aldictm = {}
        aldicts = {}
        pocket = ''
        first_strand = 0

        for i in cmtab.index.values:
            strand = cmtab['strand'][i]
            if first_strand == 0:
                first_strand = strand # 
            mdl = [cmtab['mdl-from'][i], cmtab['mdl-to'][i]]
            if strand == '+':
                seq = [cmtab['seq-from'][i], cmtab['seq-to'][i]]
            else:
                seq = [-cmtab['seq-from'][i], -cmtab['seq-to'][i]]
            aldictm[mdl[0]] = index
            aldictm[mdl[1]] = index+1
            aldicts[seq[0]] = index
            aldicts[seq[1]] = index+1
            index += 2
            if index2 < 3:
                pocket += currentseq +'_d_'+ str(index2) +'\t'+ str(cmtab['seq-from'][i]) +'\t'+ str(cmtab['seq-to'][i]) +'\t'+ cmtab['strand'][i] +'\tsimple\n' # pocket = all hits are independent (paralogs)
                index2 += 1

        mdlal = list(aldictm.keys())
        mdlal.sort()
        mdlalc = []
        for i in mdlal:
            mdlalc.append(aldictm[i]) # order of start and end of hits acording to the model
    
        seqal = list(aldicts.keys())
        seqal.sort()
        seqalc = []
        for i in seqal:
            seqalc.append(aldicts[i]) # order of start and end of hits acording in the subject sequence
    
        print(currentseq) # add info about the current contig to stdout
        print(mdlalc)
        print(seqalc)
        
        if mdlalc != seqalc: # if the order of start and end of hits is different, there might be duplications
            print('Putative duplication')
            outfile1.write(currentseq +'\n'+'mdl: '+ '.'.join(str(x) for x in mdlalc) +'\n'+'seq: '+ '.'.join(str(x) for x in seqalc) +'\n'+'\n')
            if len(mdlalc) != len(seqalc):
                print('Paralogs with exact overlap')
                outfmap.write(pocket) # paralogs 
            else:
                print('Paralogs with unexact overlap')
                outfmap.write(pocket) # paralogs
        
        elif len(seqalc) > 2: # if the order of start and end of hits is the same, it is a fragmented sequence. Single alignments are too short, thus discarded.
            print('Most likely fragmented and will be assembled')
            index = 1
            outfile2.write(currentseq +'\t'+ str(seqal[-1]-seqal[0]) +'\n') # seqsumt
            outfmap.write(currentseq +'\t'+ str(abs(seqal[0])) +'\t'+ str(abs(seqal[-1])) +'\t'+ first_strand +'\tassembled\n')
            nseqal = seqal
            while len(nseqal) > 2:
                outfile3.write(currentseq +'_'+ str(index) +'\t'+ str(abs(nseqal[2]-nseqal[1])) +'\n')
                nseqal = nseqal[2:]
                index += 1

    outfile1.close()
    outfile2.close()
    outfile3.close()
    return outfmap


## MAIN PROGRAM

print(cmtable)

## START OUTPUT FILE FOR SEQUENCES TO EXTRACT

outfmap = open(outpath +'.seqmap', 'w')

## LOAD CM LENGTH

clen = get_clen(cmtable)
            
## LOAD CMSEARCH OUTPUT

try: # there could be no hits
    ctab = pd.read_csv(cmtable, sep=r'\s+', comment='#', header=None)
    ctab.columns = ['subject','subject-acc','model','model-acc','model-type','mdl-from','mdl-to','seq-from','seq-to', 'strand','trunc','pass','gc','bias','score','evalue','include','description']
except pd.errors.EmptyDataError: 
    outfmap.close()
    Path(outpath + '.seqsumt').touch()
    Path(outpath + '.paralg').touch()
    Path(outpath + '.insert').touch()
    sys.exit()

## GET SIZES FOR CONTIGS WITH ONE HIT

singles_raw =  get_singles(ctab, clen)

print('Number of simple alignments hit: '+ str(len(singles_raw)))

outfmap = get_sing_mdls_seqs(singles_raw, outfmap, ctab)
    
## GET SIZES FOR CONTIGS WITH MULTIPLE HITS

repeated_raw = get_repeated(ctab, singles_raw)

print('Number of contigs with more than one hit: '+ str(len(repeated_raw)))

## GET COORDINATES 

get_coords(repeated_raw, ctab)

## GET PUTATIVE DUPLICATIONS/PARALOGS

outfmap = get_paralg_seqsumt_ins(repeated_raw, outfmap, ctab)

outfmap.close()