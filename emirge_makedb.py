#!/usr/bin/env python

import hashlib
import os
import random
import re
import subprocess
import sys
from optparse import OptionParser

from Emirge.download import DownloadException, SilvaDownloader

USAGE = """usage: %prog [OPTIONS]

%prog creates a reference database and the necessay indices for use by
EMIRGE from an rRNA reference database. Without extra parameters, %prog
will 1) download the most recent SILVA SSU database, 2) filter it by sequence
length, 3) cluster at 97% sequence identity, 4) replace ambiguous bases
with random characters and 5) create a bowtie index.
"""


def cluster_fasta(vsearch_bin, filein, minlen, maxlen, clusterid):
    """Runs vsearch cluster_fast on `filein`, considering only sequences
    with `minlen` <= sequence length <= `maxlen`
    """
    filein_split = filein.split(".")
    fileout = ".".join(
        filein_split[0:-2] +
        ["ge{0}bp".format(minlen), "le{0}bp".format(maxlen), str(clusterid)] +
        filein_split[-2:-1]
    )

    if os.path.isfile(fileout):
        print ("Found existing file \"{0}\". Skipping clustering stage."
               .format(fileout))
        return fileout

    cmd = [vsearch_bin,
           "--minseqlength", str(minlen),
           "--maxseqlength", str(maxlen),
           "--fasta_width", "0",
           "--notrunclabels",
           "--centroids", fileout,
           "--cluster_fast", filein,
           "--id", str(clusterid)]

    print (" ".join(["Running: "] + cmd))
    subprocess.call(cmd)
    print "Done"

    return fileout


def pairs(lst):
    """Creates pairwise iterator on `lst`.
    I.e. if lst returns 1..5, pairs(lst) will return (1,2), (3,4), (5,"")
    """
    it = iter(lst)
    for item in it:
        try:
            yield item, it.next()
        except StopIteration:
            yield item, ""


iupac_map = {
    "": " ",
    "R": "AG",    # Purine (A or G)
    "Y": "CT",    # Pyrimidine (C, T, or U)
    "M": "CA",    # C or A
    "K": "TG",    # T, U, or G
    "W": "TA",    # T, U, or A
    "S": "CG",    # C or G
    "B": "CTG",   # C, T, U, or G (not A)
    "D": "ATG",   # A, T, U, or G (not C)
    "H": "ATC",   # A, T, U, or C (not G)
    "V": "ACG",   # A, C, or G (not T, not U)
    "N": "ACTG"   # Any base (A, C, G, T, or U)
}


def randomize_ambiguous(seq):
    """Replaces IUPAC ambiguous characters in seq with random characters,
    respecting options, i.e. "R" replaced with A or G.
    """
    # slower (but clearer) implementation:
    # return "".join([choice(iupac_map.get(x, "ACTG")) for x in seq.upper()])
    seq = seq.upper().replace("U", "T")
    list = [ok + random.choice(iupac_map.get(fix, "ACGT"))
            for ok, fix in pairs(re.split("([^ACGT])", seq))]
    return "".join(list).rstrip(" ")


def randomize_ambiguous_fasta(filein, folder=None):
    """Replaces IUPAC ambiguous character in `filein` with random choice
    from allowed characters.
    Returns output filename.
    """
    filein_split = filein.split(".")
    fileout = ".".join(filein_split[0:-1] + ["fixed"] + filein_split[-1:])
    if folder is not None:
        fileout = os.path.join(folder, os.path.basename(fileout))

    if os.path.exists(fileout):
        print ("Found existing file \"{0}\". Skipping randomization stage."
               .format(fileout))
        return fileout

    processed = 0
    total = os.path.getsize(filein)
    dots = 0
    linewidth = 77
    print "Randomizing ambiguous bases"
    print "|" + "-" * (linewidth-2) + "|"
    with open(filein, "rb") as inf, open(fileout, "wb") as outf:
        for line in inf:
            if line[0] == '>':
                outf.write(line)
            else:
                outf.write(randomize_ambiguous(line.rstrip("\n")))
                outf.write("\n")
            processed += len(line)
            dotstoprint = int(linewidth * processed / total) - dots
            sys.stdout.write("."*dotstoprint)
            sys.stdout.flush()
            dots += dotstoprint

    return fileout


def build_bowtie_index(bowtie_bin, filein):
    """Calls bowtie-build on `filein` to compute bowtie index"""
    fileout = ".".join(filein.split(".")[:-1])
    cmd = [bowtie_bin, "-o", "0", filein, fileout]
    print "Running: " + " ".join(cmd)
    subprocess.call(cmd)


def main(argv=sys.argv[1:]):
    parser = OptionParser(USAGE)
    parser.add_option(
        "-g", "--gene", dest="gene", metavar="[SSU|LSU]", type="string",
        default="SSU",
        help="build database from this gene"
    )
    parser.add_option(
        "-t", "--tmpdir", dest="tmpdir", metavar="DIR", type="string",
        default="/tmp",
        help="work directory for temporary files"
    )
    parser.add_option(
        "-r", "--release", dest="release", metavar="N", type="string",
        default="current",
        help="SILVA release number"
    )
    parser.add_option(
        "-m", "--min-len", dest="min_len", metavar="LEN", type="int",
        default=1200,
        help="minimum reference sequence length"
    )
    parser.add_option(
        "-M", "--max-len", dest="max_len", metavar="LEN", type="int",
        default=2000,
        help="maximum reference sequence length"
    )
    parser.add_option(
        "-i", "--id", dest="clusterid", metavar="FLOAT", type="float",
        default=0.97,
        help="Cluster at this identity level"
    )
    parser.add_option(
        "-k", "--keep", dest="keep", action="store_true",
        help="keep intermediary files"
    )
    parser.add_option(
        "-V", "--vsearch", dest="vsearch", metavar="FILE", type="string",
        default="vsearch",
        help="path to vsearch binary"
    )
    parser.add_option(
        "-B", "--bowtie-build", dest="bowtie", metavar="FILE", type="string",
        default="bowtie-build",
        help="path to bowtie-build binary"
    )

    (options, args) = parser.parse_args(argv)

    try:
        downloader = SilvaDownloader()
        silva_fasta = downloader.run(options.gene, options.release,
                                     options.tmpdir)
        clustered_fasta = cluster_fasta(options.vsearch, silva_fasta,
                                        options.min_len, options.max_len,
                                        options.clusterid)
        randomized_fasta = randomize_ambiguous_fasta(clustered_fasta,
                                                     folder="")
        build_bowtie_index(options.bowtie, randomized_fasta)
        if (not options.keep):
            os.unlink(silva_fasta)
            os.unlink(clustered_fasta)

    except DownloadException as e:
        print(e.args[0])
        exit(1)


if __name__ == '__main__':
    main()
