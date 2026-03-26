#!/usr/bin/env python3

import sys

try:
    fasta_name = sys.argv[1]
    seqmap_name = sys.argv[2]
    outfilename = sys.argv[3]
    minlength = sys.argv[4]
    
except:
    print("""
Retrieve sequences from a "cmtable.seqmap" file.
Usage:
get_cmsequences.py SEQUENCES.fna CMTABLE.seqmap OUTPUT_NAME.fna min_length
    """)
    sys.exit()
    

# WATSON-CRICK DICTIONARY
    
wcdict = {
    'A':'T',
    'T':'A',
    'U':'A',
    'C':'G',
    'G':'C',
    'a':'t',
    't':'a',
    'u':'a',
    'c':'g',
    'g':'c',
    'N':'N',
    'n':'n',
    'Y':'R',
    'R':'Y',
    'W':'W',
    'S':'S',
    'K':'M',
    'M':'K',
    'B':'V',
    'V':'B',
    'D':'H',
    'H':'D'
    }
    
# FUNCTION DEFINITION

def get_seqmap_dict(seqmap):
    seqmap_dict = {}
    for line in seqmap:
        linesep = line.split('\t')
        if '_d_' in linesep[0]:
            linesepz = linesep[0].split('_d_')
            key = linesepz[0]
            linesep[0] = "_".join(linesepz)
        else:
            key = linesep[0]
        if key not in seqmap_dict.keys():
            seqmap_dict[key] = []
        seqmap_dict[key].append(linesep)
    return seqmap_dict

def get_sequences(seqmap_dict, fasta_name, outfile, minlength):
    # First pass: load all sequences that appear in the seqmap into memory.
    # This correctly handles multi-line FASTA files by concatenating all
    # sequence lines belonging to each header.
    sequences = {}
    currkey = ''
    collecting = False
    with open(fasta_name) as fileobject:
        for line in fileobject:
            line = line.rstrip('\n')
            if line.startswith('>'):
                header = line[1:]
                if header in seqmap_dict:
                    currkey = header
                    sequences[currkey] = []
                    collecting = True
                else:
                    collecting = False
            elif collecting:
                sequences[currkey].append(line)
    # Join the collected lines into full sequences
    for key in sequences:
        sequences[key] = ''.join(sequences[key])

    # Second pass: extract subsequences using the seqmap coordinates
    for currkey, c_seq in sequences.items():
        for seqlist in seqmap_dict[currkey]:
            if seqlist[3] == '+':
                d_seq = c_seq[int(seqlist[1]) : int(seqlist[2])]
            else:
                sub = c_seq[int(seqlist[2]) : int(seqlist[1])]
                sub = sub[::-1]
                d_seq = ''
                for nuc in sub:
                    d_seq += wcdict[nuc]
            if len(d_seq) >= int(minlength):
                seqfrom = str(min([int(i) for i in seqlist[1:3]]))
                seqto = str(max([int(i) for i in seqlist[1:3]]))
                outfile.write('>'+ currkey +'|'+ seqfrom +'-'+ seqto +'|strand_'+ seqlist[3] +'|'+ seqlist[4] +'\n'+ d_seq +'\n')
    return outfile
            

# MAIN PROGRAM

# Prepare output
outfile = open(outfilename, 'w')

# Load coordinates
handle = open(seqmap_name, 'r')
seqmap = handle.read().split('\n')[:-1]
handle.close()

seqmap_dict = get_seqmap_dict(seqmap)

# Extract sequences
outfile = get_sequences(seqmap_dict, fasta_name, outfile, minlength)
outfile.close()
        
