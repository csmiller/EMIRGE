#!/usr/bin/env python

import urllib2
import urllib
import re
import sys, os
import hashlib
from optparse import OptionParser, OptionGroup
import random
import subprocess

USAGE = \
"""usage: %prog [OPTIONS]

%prog creates a reference database and the necessay indices for use by EMIRGE
from an rRNA reference database. Without extra parameters, %prog will 1) download
the most recent SILVA SSU database, 2) filter it by sequence length, 3) cluster at
97% sequence identity, 4) replace ambiguous bases with random characters and 5) 
create a bowtie index.
"""

SILVA_BASEURL = "https://ftp.arb-silva.de"
SILVA_SSUFILE = "SILVA_{0}_SSURef_Nr99_tax_silva_trunc.fasta.gz"
SILVA_LSUFILE = "SILVA_{0}_LSURef_tax_silva_trunc.fasta.gz"

iupac_map = {
    "": " ",
    "R": "AG", #	Purine (A or G)
    "Y": "CT", #	Pyrimidine (C, T, or U)
    "M": "CA", #	C or A
    "K": "TG", #	T, U, or G
    "W": "TA", #	T, U, or A
    "S": "CG", #	C or G
    "B": "CTG",#	C, T, U, or G (not A)
    "D": "ATG",#	A, T, U, or G (not C)
    "H": "ATC",#	A, T, U, or C (not G)
    "V": "ACG",#	A, C, or G (not T, not U)
    "N": "ACTG" }#	Any base (A, C, G, T, or U)

class DownloadException(Exception):
    pass

def fetch_url(url):
    try:
        return urllib2.urlopen(url).read()
    except urllib2.URLError as e:
        raise DownloadException(
            "Unable to fetch \"{0}\":\n \"{1}\"".format(license_url, e.reason)
            )

def print_progress(block, blocksize, total):
    blocks = int((total-1)/blocksize) + 1
    linewidth = min(77, blocks)
    if block == 0:
        print "Downloading file of size {0}:".format(total)
        if (linewidth > 1):
            print "|" + "-" * (linewidth-2) + "|"
    else:
        if (block % (blocks/linewidth) == 0):
            sys.stdout.write(".")
            sys.stdout.flush()
    
def download_url(url, folder):
    print "Downloading \"{0}\" to \"{1}\"".format(url, folder)
        
    filename = os.path.join(folder, url.split("/")[-1])
    if os.path.isfile(filename):
        existing_file_size = os.path.getsize(filename)
        remote_file_size = urllib.urlopen(url).info().getheaders("Content-Length")[0]
        if int(existing_file_size) == int(remote_file_size):
            print "Found existing file matching remote size. Skipping download."
            return filename
        else:
            msg = "Local file with size {0} does not match remote file size {1}."
            print msg.format(existing_file_size, remote_file_size)
            
    (filename, header) = urllib.urlretrieve(url, filename, print_progress)
    
    return filename

def compute_file_md5(filename):
    hash = hashlib.md5()
    with open(filename, "rb") as file:
        for chunk in iter(lambda: file.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()
    
def silva_find_most_recent_version():
    url = "/".join([SILVA_BASEURL, "current", "Exports"])
    files = re.findall("<a href=\"SILVA_([[0-9_]*)_",fetch_url(url))
    return files[0]

def silva_confirm_license(release):
    license_url = "/".join([SILVA_BASEURL, "release_" + release.replace(".","_"),
                            "Exports", "LICENSE.txt"])
    license_txt = fetch_url(license_url)
                           
    print("""
The SILVA database is published under a custom license. To proceed, you need to read this
license and agree with its terms:
    """)
    print ("> "+license_txt.replace("\n", "\n> ").rstrip("\n> "))
    print ""

    answer = raw_input("Do you agree to these terms? [yes|NO]")
    if (answer != "yes"):
        raise DownloadException("Unable to continue")
        

def silva_download_fasta(release, gene, tmpdir):
    filename = SILVA_SSUFILE if gene.upper() == "SSU" else SILVA_LSUFILE
    filename = filename.format(release)
    database_url = "/".join([SILVA_BASEURL, "release_" + release.replace(".","_"),
                             "Exports", filename])
    silva_fasta = download_url(database_url, tmpdir)
    silva_md5   = download_url(database_url + ".md5", tmpdir)
    print "Verifying file...",
    with open(silva_md5, "rb") as file:
          silva_md5 = file.read().split(" ")[0]
    local_md5 = compute_file_md5(silva_fasta)
    if (local_md5 == silva_md5):
        print "OK"
        return silva_fasta
    else:
        print "FAILED"
        raise DownloadException("Corrupted File?!")

def cluster_fasta(binary, filein, minlen, maxlen, clusterid):
    filein_split = filein.split(".")
    fileout = ".".join(filein_split[0:-2] +
                       ["ge{0}bp".format(minlen), "le{0}bp".format(maxlen), str(clusterid)] +
                       filein_split[-2:-1]
                      )

    if os.path.isfile(fileout):
        print "Found existing file \"{0}\". Skipping clustering stage.".format(fileout)
        return fileout
    
    cmd = [binary,
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
    it   = iter(lst)
    for item in it:
        try:
            yield item, it.next()
        except StopIteration:
            yield item, ""

def randomize_ambiguous(seq):
    #return "".join([choice(iupac_map.get(x, "ACTG")) for x in seq.upper()])
    seq = seq.upper().replace("U","T")
    list = [ ok + random.choice(iupac_map.get(fix, "ACGT"))
             for ok, fix in pairs(re.split("([^ACGT])", seq)) ]
    return "".join(list).rstrip(" ")

def randomize_ambiguous_fasta(filein, folder = None):
    filein_split = filein.split(".")
    fileout = ".".join(filein_split[0:-1] + ["fixed"] + filein_split[-1:])
    if folder != None:
        fileout = os.path.join(folder, os.path.basename(fileout))

    if os.path.exists(fileout):
        print "Found existing file \"{0}\". Skipping randomization stage.".format(fileout)
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
            dots+=dotstoprint
            
    return fileout

def build_bowtie_index(binary, filein):
    fileout = ".".join(filein.split(".")[:-1])
    cmd = [binary, "-o", "0", filein, fileout]
    print "Running: " + " ".join(cmd)
    subprocess.call(cmd)

def main(argv = sys.argv[1:]):
    parser = OptionParser(USAGE)
    parser.add_option("-g", "--gene", dest="gene", metavar="[SSU|LSU]", type="string",
                      default="SSU",
                      help="build database from this gene")
    parser.add_option("-t", "--tmpdir", dest="tmpdir", metavar="DIR", type="string",
                      default="/tmp",
                      help="work directory for temporary files")
    parser.add_option("-r", "--release", dest="release", metavar="N", type="string",
                      default="current",
                      help="SILVA release number")
    parser.add_option("-m", "--min-len", dest="min_len", metavar="LEN", type="int",
                      default=1200,
                      help="minimum reference sequence length")
    parser.add_option("-M", "--max-len", dest="max_len", metavar="LEN", type="int",
                      default=2000,
                      help="maximum reference sequence length")
    parser.add_option("-i", "--id", dest="clusterid", metavar="FLOAT", type="float",
                      default=0.97,
                      help="Cluster at this identity level")
    parser.add_option("-k", "--keep", dest="keep", action="store_true",
                      help="keep intermediary files")
    parser.add_option("-V", "--vsearch", dest="vsearch", metavar="FILE", type="string",
                      default="vsearch",
                      help="path to vsearch binary")
    parser.add_option("-B", "--bowtie-build", dest="bowtie", metavar="FILE", type="string",
                      default="bowtie-build",
                      help="path to bowtie-build binary")

    (options, args) = parser.parse_args(argv)

    if (options.release == "current"):
        options.release = silva_find_most_recent_version()

    try:
        silva_confirm_license(options.release)
        silva_fasta      = silva_download_fasta(options.release, options.gene, options.tmpdir)
        clustered_fasta  = cluster_fasta(options.vsearch, silva_fasta,
                                         options.min_len, options.max_len, options.clusterid)
        randomized_fasta = randomize_ambiguous_fasta(clustered_fasta, folder="")
        build_bowtie_index(options.bowtie, randomized_fasta)
        if (not options.keep):
            os.unlink(silva_fasta)
            os.unlink(silva_fasta+".md5")
            os.unlink(clustered_fasta)

    except DownloadException as e:
        print(e.args[0])
        exit(1)



if __name__ == '__main__':
    main()
    
