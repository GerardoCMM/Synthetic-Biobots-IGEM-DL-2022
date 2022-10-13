from Bio import SeqIO

import pandas as pd

genes = ["Pn8.2617","Pn2.84","Pn1.1317","Pn7.1626","Pn16.1198","Pn4.3222","Pn2.2377","Pn6.2477"]

seqs = SeqIO.parse("../Genome/Piper_nigrum.cds",format="fasta")

for seq in seqs:
    if seq.id in genes:
        print(">"+seq.id)
        print(seq.seq)