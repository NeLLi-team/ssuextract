#!/usr/bin/env python3

import sys
import pandas as pd

try:
    cmtable = sys.argv[1]
    outpath = sys.argv[2]

except:
    print("""
Retrieve information from a cmsearch output table.
Usage:
get_cmstats.py cmtable prefix
    """)
    sys.exit()

    
## Function definitions

def get_singles(ctab):
    singles = ctab['subject'].value_counts().loc[ctab['subject'].value_counts()==1].index.values
    passed_singles = list(ctab.loc[(ctab['subject'].isin(singles))&(ctab['include']=='!'), 'subject'])
    return(passed_singles)

def get_sing_mdls_seqs(singles_raw, outfmap, ctab):
    # outfile1 = open(outpath +'.singles.mdls','w') # length of each hit according to the model
    # outfile2 = open(outpath +'.singles.seqs','w') # length of each hit in the query sequence
    for subject in singles_raw:
        cmtab = ctab.loc[ctab['subject']==subject]
        mdl = [cmtab['mdl-from'].item(), cmtab['mdl-to'].item()]
        # outfile1.write(subject +'\t'+ str((max(mdl) - min(mdl))) +'\n')
        seq = [cmtab['seq-from'].item(), cmtab['seq-to'].item()]
        # outfile2.write(subject +'\t'+ str((max(seq)) - min(seq)) +'\n')
        outfmap.write(subject +'\t'+ str(cmtab['seq-from'].item()) +'\t'+ str(cmtab['seq-to'].item()) +'\t'+ cmtab['strand'].item() +'\n')    
    # outfile1.close()
    # outfile2.close()
    return outfmap
    
## GET SIZES FOR CONTIGS WITH MULTIPLE HITS

def get_repeated(ctab):
    repeated_raw = list(ctab['subject'].value_counts().loc[ctab['subject'].value_counts()>1].index.values)
    return repeated_raw

def get_rep_mdls_seqs_sums(repeated_raw, ctab):
    mdlsum = 0
    seqsum = 0
    # outfile1 = open(outpath +'.mdls','w')    # length of each hit according to the model
    # outfile2 = open(outpath +'.seqs','w')    # length of each hit in the query sequence
    # outfile3 = open(outpath +'.mdlsums','w') # summed lengths of hits in each contig according to the model
    # outfile4 = open(outpath +'.seqsums','w') # summed lengths of hits in each contig

    for subject in repeated_raw:
        cmtab = ctab.loc[ctab['subject']==subject]
        index=1
        for i in cmtab.index.values:
            mdl = [cmtab['mdl-from'][i], cmtab['mdl-to'][i]]
            mdlsum += (max(mdl) - min(mdl))
            # outfile1.write(subject +'_'+ str(index) +'\t'+ str(max(mdl) - min(mdl))+'\n')
            seq = [cmtab['seq-from'][i], cmtab['seq-to'][i]]
            seqsum += ((max(seq)) - min(seq))
            # outfile2.write(subject +'_'+ str(index) +'\t'+ str(max(seq) - min(seq))+'\n')
            index += 1
    
        # outfile3.write(subject +'\t'+ str(mdlsum) +'\n')
        mdlsum = 0
        # outfile4.write(subject +'\t'+ str(seqsum) +'\n')
        seqsum = 0

    # outfile1.close()
    # outfile2.close()
    # outfile3.close()
    # outfile4.close()
    return 0

## GET COORDINATES 

def get_coords(repeated_raw, ctab):
    # outfile = open(outpath +'.coords', 'w')
    
    for currentseq in repeated_raw:
        # outfile.write(currentseq +'\n')
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
    
        # outfile.write('Model:    '+ mdlals +'\n')

        seqal = list(aldicts.keys())
        seqal.sort()

        seqals = ''

        for i in seqal:
            seqals += '..'
            seqals += aldicts[i]
            seqals += '..'
    
        # outfile.write('Sequence: '+ seqals +'\n')
        # outfile.write('\n')
    
    # outfile.close()
    return 0

## GET PUTATIVE DUPLICATIONS/PARALOGS

def get_paralg_seqsumt_ins(repeated_raw, outfmap, ctab):
    # outfile1 = open(outpath +'.paralg', 'w')    # sequences with putative duplications
    # outfile2 = open(outpath +'.seqsumt', 'w')   # total lengths of composie hits including insertions
    # outfile3 = open(outpath +'.insert', 'w')    # lenghts of each insertion between hit pairs in each contig

    for currentseq in repeated_raw: 
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
                first_strand = strand
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
                pocket += currentseq +'_d_'+ str(index2) +'\t'+ str(cmtab['seq-from'][i]) +'\t'+ str(cmtab['seq-to'][i]) +'\t'+ cmtab['strand'][i] +'\n'
                index2 += 1

        mdlal = list(aldictm.keys())
        mdlal.sort()
        mdlalc = ''
        for i in mdlal:
            mdlalc += str(aldictm[i])
    
        seqal = list(aldicts.keys())
        seqal.sort()
        seqalc = ''
        for i in seqal:
            seqalc += str(aldicts[i])
    
        if mdlalc != seqalc:
            # outfile1.write(currentseq +'\n'+'mdl: '+ mdlalc +'\n'+'seq: '+ seqalc +'\n'+'\n')
            if mdlalc == '1342' or mdlalc == '342' or mdlalc == '34' or mdlalc == '134':
                outfmap.write(pocket)
            elif mdlalc == '3142' or mdlalc == '1324':
                # outfile2.write(currentseq +'\t'+ str(seqal[-1]-seqal[0]) +'\n')
                outfmap.write(currentseq +'\t'+ str(abs(seqal[0])) +'\t'+ str(abs(seqal[-1])) +'\t'+ first_strand +'\n')
        else:
            index = 1
            # outfile2.write(currentseq +'\t'+ str(seqal[-1]-seqal[0]) +'\n')
            outfmap.write(currentseq +'\t'+ str(abs(seqal[0])) +'\t'+ str(abs(seqal[-1])) +'\t'+ first_strand +'\n')
            nseqal = seqal
            while len(nseqal) > 2:
                # outfile3.write(currentseq +'_'+ str(index) +'\t'+ str(abs(nseqal[2]-nseqal[1])) +'\n')
                nseqal = nseqal[2:]
                index += 1

    # outfile1.close()
    # outfile2.close()
    # outfile3.close()
    return outfmap


## MAIN PROGRAM

## START OUTPUT FILE FOR SEQUENCES TO EXTRACT

outfmap = open(outpath +'.seqmap', 'w')

## LOAD CMSEARCH OUTPUT

try: # there could be no hits
    ctab = pd.read_csv(cmtable, sep='\s+', comment='#', header=None)
    ctab.columns = ['subject','subject-acc','model','model-acc','model-type','mdl-from','mdl-to','seq-from','seq-to', 'strand','trunc','pass','gc','bias','score','evalue','include','description']
except pd.errors.EmptyDataError: 
    outfmap.close()
    sys.exit()

## GET SIZES FOR CONTIGS WITH ONE HIT

singles_raw =  get_singles(ctab)

print('Number of contigs with one single hit: '+ str(len(singles_raw)))

outfmap = get_sing_mdls_seqs(singles_raw, outfmap, ctab)
    
## GET SIZES FOR CONTIGS WITH MULTIPLE HITS

repeated_raw = get_repeated(ctab)

print('Number of contigs with more than one hit: '+ str(len(repeated_raw)))

get_rep_mdls_seqs_sums(repeated_raw, ctab)

## GET COORDINATES 

get_coords(repeated_raw, ctab)

## GET PUTATIVE DUPLICATIONS/PARALOGS

outfmap = get_paralg_seqsumt_ins(repeated_raw, outfmap, ctab)

outfmap.close()