import os
import sys
import glob
from Bio import SeqIO


# directory with fna files
fnain = sys.argv[1]
fnaout = sys.argv[2]


def rename_headers_merge(fnain, fnaout):
    """
    Add filename to header in first field (separated by "|").
    Write all faa into a single file.
    """
    with open(fnaout, "w") as outfile:
        for seq_record in SeqIO.parse(fnain, "fasta"):
            filename = fnain.split("/")[-1].split(".")[0]

            if "|" in seq_record.id and seq_record.id.split("|")[0] != filename:
                header = "{}|{}".format(filename, seq_record.id.split()[0])
            else:
                header = seq_record.id.split()[0]

            outfile.write(">{}\n{}\n".format(header, seq_record.seq))


rename_headers_merge(fnain, fnaout)
