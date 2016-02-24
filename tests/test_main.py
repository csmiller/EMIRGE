#!/usr/bin/env python

import sys
from os import path
sys.path.insert(1, path.dirname(path.dirname(path.abspath(__file__))))
from emirge_amplicon import main

base = path.dirname(path.abspath(__file__))

outdir = path.join(base, "ten_seq_out")

base = path.join(base, "test_data")
in1 = path.join(base, "ten_seq_community_000_50K_L150_I350.1.fastq")
in2 = path.join(base, "ten_seq_community_000_50K_L150_I350.2.fastq")
ref = path.join(base, "twenty_seq_database.fasta")
bt = path.join(base, "twenty_seq_database")

main(argv=[
    outdir,
    "-1", in1, "-2", in2, "-f", ref, "-b", bt,
    "--phred33",
    "-i", "350",
    "-s", "500",
    "-l", "150",
    "-j", "0.97",
    "-n", "5",
    "-a", "8",
    "-o", "indels"
])


