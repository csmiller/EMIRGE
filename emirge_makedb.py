#!/usr/bin/env python

import hashlib
import os
import random
import re
import subprocess
import sys
import urllib.request
import urllib.parse
import urllib.error
from optparse import OptionParser

USAGE = """usage: %prog [OPTIONS]

%prog creates a reference database and the necessay indices for use by
EMIRGE from an rRNA reference database. Without extra parameters, %prog
will 1) download the most recent SILVA SSU database, 2) filter it by sequence
length, 3) cluster at 97% sequence identity, 4) replace ambiguous bases
with random characters and 5) create a bowtie index.

Requires vsearch executable can be found in path for clustering.
https://github.com/torognes/vsearch

Requires bowtie-build (from bowtie version 1) can be found in path

"""
def INFO(s):
    sys.stdout.write(s+'\n')
    return

# these are back-ported from EMIRGE2 for the time being
class DownloadException(Exception):
    """Used to pass the URL downwards for urllib(2) exceptions"""
    pass


class BaseDownloader(object):
    """Base class for downloaders"""

    target_folder = "."

    @staticmethod
    def _print_progress(block, block_size, total):
        """Print progress of download to stderr.

        This is a callback function for urllib.

        Args:
            block:     number of block downloaded
            block_size: size of a block in bytes
            total:     total file size in bytes
        """
        blocks = int((total - 1) / block_size) + 1
        line_width = min(77, blocks)
        if block == 0:
            print(("Downloading file of size {0}:".format(total)))
            if line_width > 1:
                print(("|" + "-" * (line_width - 2) + "|"))
        else:
            if block % (blocks / line_width) == 0:
                sys.stderr.write(".")
                sys.stderr.flush()

    @staticmethod
    def compute_file_md5(filename):
        """Computes the MD5 checksum of a file"""
        INFO("Computing MD5 checksum of {}...".format(filename))

        hash = hashlib.md5()
        with open(filename, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash.update(chunk)
        return hash.hexdigest()

    @staticmethod
    def fetch_url(url):
        """Fetch document at url

        Args:
            url: URL to fetch.

        Returns:
            contents of file at url

        Raises:
            DownloadException: if download failed
        """
        INFO("Fetching {url}\n".format(url=url))

        try:
            return urllib.request.urlopen(url).read()
        except urllib.error.URLError as e:
            raise DownloadException(
                "Unable to fetch \"{0}\":\n \"{1}\"".format(url, e.reason)
            )

    def download_file(self, url, folder=None, filename=None):
        """Downloads URL as file into folder. Returns filename.

        Skips download if the file already exists and has the same size.

        Args:
            url:      URL to download
            folder:   target folder  (default ".")
            filename: specify target filename (default: guess from url)

        Returns:
            name of file

        """
        INFO("Downloading {url}".format(url=url))
        if folder is None:
            folder = self.target_folder
        if filename is None:
            filename = url.split("/")[-1].split("?")[0]
        local_file = os.path.join(folder, filename)

        INFO("Local filename is \"{1}\"".format(url, local_file))

        # check if file exists
        if os.path.isfile(filename):
            existing_file_size = int(os.path.getsize(filename))
            url_info = urllib.request.urlopen(url).info()
            try:
                remote_file_size = int(url_info.getheaders("Content-Length")[0])
            except AttributeError:
                remote_file_size = int([
                    h[1] for h in
                    url_info._headers
                    if h[0] == 'Content-Length'][0])
            if existing_file_size == 0:
                INFO("Found existing file of size 0. Re-downloading...")
            elif existing_file_size == remote_file_size:
                INFO("Found existing file matching remote size. "
                     "Skipping download.")
                self.check_file(url, filename)
                return filename
            else:
                INFO("Local file with size {0} does not match "
                     "remote file size {1}."
                     .format(existing_file_size, remote_file_size))

        (filename, _) = urllib.request.urlretrieve(url, filename,
                                           self._print_progress)
        sys.stderr.write("\n")  # clear last endline from progress bar
        self.check_file(url, filename)
        return filename

    def check_file(self, url, filename):
        url = url.split("?", 1)
        url[0] += ".md5"
        url = "?".join(url)
        try:
            remote_md5 = self.fetch_url(url).decode('utf8').split(" ")[0].strip()
        except DownloadException:
            INFO("No MD5 file found on remote")
            return
        local_md5 = self.compute_file_md5(filename)
        if local_md5 != remote_md5:
            raise DownloadException(
                "MD5 sum mismatch: remote {} - {} != local {} - {}".format(
                    url, remote_md5,
                    filename, local_md5
            ))
        else:
            INFO("Verified MD5 sum for {}".format(filename))


class SilvaDownloader(BaseDownloader):
    # ftp://ftp.arb-silva.de/release_123_1/Exports/
    # SILVA_123.1_SSURef_Nr99_tax_silva_trunc.fasta.gz
    # SILVA_123.1_LSURef_tax_silva_trunc.fasta.gz
    BASEURL = "https://ftp.arb-silva.de/{reldir}/Exports/"
    FILENAMES = {
        "SSU": "SILVA_{rel}_SSURef_Nr99_tax_silva_trunc.fasta.gz",
        "LSU": "SILVA_{rel}_LSURef_tax_silva_trunc.fasta.gz",
        "LICENSE": "LICENSE.txt",
        "LISTING": ""
    }

    def __init__(self):
        self.license_confirmed = False

    def get_url(self, name, release=None):
        try:
            return (self.BASEURL + self.FILENAMES[name.upper()]).format(
                reldir="release_" + release.replace(".", "_") if release \
                    else "current",
                rel=release
            )
        except KeyError:
            raise DownloadException("Gene \"{}\" not available at SILVA")

    def download_file(self, gene, release):
        return super(SilvaDownloader, self).download_file(
            self.get_url(gene, release)
        )

    def get_current_version(self):
        listing_url = self.get_url("LISTING")
        listing = self.fetch_url(listing_url).decode('utf8')
        pattern = self.FILENAMES["SSU"].format(rel='([0-9.]+)')
        try:
            version = re.search(pattern, listing).group(1)
        except AttributeError:
            raise DownloadException(
                "Could not find entry matching regex {} on {}".
                    format(pattern, listing_url))
        return version

    def confirm_license(self, release=None):
        if self.license_confirmed == True:
            return

        license_url = self.get_url("LICENSE", release)
        license = self.fetch_url(license_url)

        print(("""
The SILVA database is published under a custom license. To proceed,
you need to read this license and agree with its terms:

Contents of \"{url}\":
            """.format(url=license_url)))
        print(("> " + license.decode('utf8').replace("\n", "\n> ").rstrip("\n> ")))
        print("")

        answer = input("Do you agree to these terms? [yes|NO]")
        if (answer.lower() != "yes"):
            raise DownloadException(
                "Unable to continue -- license not accepted")
        self.license_confirmed = True

    def run(self, gene="SSU", release="current", tmpdir=None):
        if release == "current":
            release = self.get_current_version()
        self.confirm_license(release)
        return self.download_file(gene, release)

def cluster_fasta(vsearch_bin, filein, minlen, maxlen, clusterid, threads):
    """Runs vsearch cluster_fast on `filein`, considering only sequences
    with `minlen` <= sequence length <= `maxlen`
    """
    filein_split = filein.split(".")
    fileout = ".".join(
        filein_split[0:-2] +
        ["ge{0}bp".format(minlen), "le{0}bp".format(maxlen), str(clusterid)] +
        filein_split[-2:-1]
    )

    if os.path.isfile(fileout) and os.path.getsize(fileout) >0:
        print(("Found existing file \"{0}\". Skipping clustering stage."
               .format(fileout)))
        return fileout

    cmd = [vsearch_bin,
           "--threads", str(threads),
           "--minseqlength", str(minlen),
           "--maxseqlength", str(maxlen),
           "--fasta_width", "0",
           "--notrunclabels",
           "--centroids", fileout,
           "--cluster_fast", filein,
           "--id", str(clusterid)]

    print((" ".join(["Running: "] + cmd)))
    subprocess.call(cmd)
    print("Done")

    return fileout


def pairs(lst):
    """Creates pairwise iterator on `lst`.
    I.e. if lst returns 1..5, pairs(lst) will return (1,2), (3,4), (5,"")
    """
    it = iter(lst)
    for item in it:
        try:
            yield item, next(it)
        except StopIteration:
            yield item, ""



def randomize_ambiguous(seq):
    """Replaces IUPAC ambiguous characters in seq with random characters,
    respecting options, i.e. "R" replaced with A or G.
    """
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

    if os.path.exists(fileout) and os.path.getsize(fileout) > 0:
        print(("Found existing file \"{0}\". Skipping randomization stage."
               .format(fileout)))
        return fileout

    processed = 0
    total = os.path.getsize(filein)
    dots = 0
    linewidth = 77
    print("Randomizing ambiguous bases")
    print("|" + "-" * (linewidth-2) + "|")
    with open(filein, "rb") as inf, open(fileout, "wb") as outf:
        for line in inf:
            if line[0] == '>':
                outf.write(line)
            else:
                outf.write(randomize_ambiguous(line.rstrip("\n")))
                outf.write("\n")
            processed += len(line)
            dotstoprint = int(linewidth * processed / total) - dots
            sys.stderr.write("."*dotstoprint)
            sys.stderr.flush()
            dots += dotstoprint
    sys.stderr.write("\n")
    return fileout


def build_bowtie_index(bowtie_bin, filein):
    """Calls bowtie-build on `filein` to compute bowtie index"""
    fileout = ".".join(filein.split(".")[:-1])
    cmd = [bowtie_bin, "-o", "0", filein, fileout]
    print("Running: " + " ".join(cmd))
    subprocess.call(cmd)


def main(argv=sys.argv[1:]):
    parser = OptionParser(USAGE)
    parser.add_option(
        "-g", "--gene", dest="gene", metavar="[SSU|LSU]", type="string",
        default="SSU",
        help="build database from this gene (SSU=Small Subunit rRNA; LSU=Large Subunit rRNA)\ndefault = %default"
    )
    parser.add_option(
        "-p", "--threads", dest="threads", type="int",
        default=0,
        help="number of threads to use for vsearch clustering of database (default = use all available)"
    )
    parser.add_option(
        "-t", "--tmpdir", dest="tmpdir", metavar="DIR", type="string",
        default="/tmp",
        help="working directory for temporary files (default = %default)"
    )
    parser.add_option(
        "-r", "--release", dest="release", metavar="N", type="string",
        default="current",
        help="SILVA release number (default: current SILVA release)"
    )
    parser.add_option(
        "-m", "--min-len", dest="min_len", metavar="LEN", type="int",
        default=1200,
        help="minimum reference sequence length (default = %default)"
    )
    parser.add_option(
        "-M", "--max-len", dest="max_len", metavar="LEN", type="int",
        default=2000,
        help="maximum reference sequence length (default = %default)"
    )
    parser.add_option(
        "-i", "--id", dest="clusterid", metavar="FLOAT", type="float",
        default=0.97,
        help="Cluster at this fractional identity level (default = %default)"
    )
    parser.add_option(
        "-k", "--keep", dest="keep", action="store_true",
        help="keep intermediary files (default: do not keep)"
    )
    parser.add_option(
        "-V", "--vsearch", dest="vsearch", metavar="FILE", type="string",
        default="vsearch",
        help="path to vsearch binary (default: look in $PATH)"
    )
    parser.add_option(
        "-B", "--bowtie-build", dest="bowtie", metavar="FILE", type="string",
        default="bowtie-build",
        help="path to bowtie-build binary (default: look in $PATH)"
    )
    parser.add_option(
        "--silva-license-accepted", dest="license",
        action="store_true", default=False,
        help="I have read and accepted the SILVA license."
    )

    (options, args) = parser.parse_args(argv)

    try:
        downloader = SilvaDownloader()
        downloader.license_confirmed = options.license
        silva_fasta = downloader.run(options.gene, options.release,
                                     options.tmpdir)
        clustered_fasta = cluster_fasta(options.vsearch, silva_fasta,
                                        options.min_len, options.max_len,
                                        options.clusterid, options.threads)
        randomized_fasta = randomize_ambiguous_fasta(clustered_fasta,
                                                     folder="")
        build_bowtie_index(options.bowtie, randomized_fasta)
        if (not options.keep):
            os.unlink(silva_fasta)
            os.unlink(clustered_fasta)

    except DownloadException as e:
        print((e.args[0]))
        exit(1)


if __name__ == '__main__':
    main()
