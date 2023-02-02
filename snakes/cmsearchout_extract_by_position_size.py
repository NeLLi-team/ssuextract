import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

cmsearchoutput = sys.argv[1]
contigsfile = sys.argv[2]
seqout = sys.argv[3]
minlength = int(sys.argv[4])

def get_hits(cmsearchhits, minlength):
    filteredhits_dict = {}
    for line in cmsearchhits:
        if not line.startswith("#"):
            contigid = line.split()[0]
            start = line.split()[7]
            stop = line.split()[8]
            strand  = line.split()[9]
            start, stop = int(start), int(stop)
            if (start > stop and start - stop > minlength) or (stop > start and stop - start > minlength):
                filteredhits_dict[contigid] = [min(start, stop), max(start, stop), strand]
    return filteredhits_dict

def extract_hits(contigsfile, filteredhits_dict, outfile):
    for seq_record in SeqIO.parse(contigsfile, "fasta"):
        if seq_record.id in filteredhits_dict:
            start, stop, strand = filteredhits_dict[seq_record.id]
            header = ">{}|{}-{}|strand {}".format(seq_record.id, start, stop, strand)
            sequence = seq_record.seq[start:stop]
            if strand == "-":
                sequence = sequence.reverse_complement()
            outfile.write("{}\n{}\n".format(header, sequence))

with open(cmsearchoutput) as cmsearchhits, open(seqout, "w") as outfile:
    extract_hits(contigsfile, get_hits(cmsearchhits, minlength), outfile)
