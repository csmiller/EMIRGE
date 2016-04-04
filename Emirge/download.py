"""
helpers for downloading things
"""

import os
import re
import sys
import urllib
import urllib2
from string import lower

import time

from Emirge.log import INFO


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
        blocks = int((total-1) / block_size) + 1
        line_width = min(77, blocks)
        if block == 0:
            print("Downloading file of size {0}:".format(total))
            if line_width > 1:
                print("|" + "-" * (line_width-2) + "|")
        else:
            if block % (blocks/line_width) == 0:
                sys.stderr.write(".")
                sys.stderr.flush()

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
        try:
            return urllib2.urlopen(url).read()
        except urllib2.URLError as e:
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
        if folder is None:
            folder = self.target_folder
        if filename is None:
            filename = url.split("/")[-1].split("?")[0]
        local_file = os.path.join(folder, filename)

        INFO('Downloading "{0}" -> "{1}"'.format(url, local_file))

        # check if file exists
        if os.path.isfile(filename):
            existing_file_size = int(os.path.getsize(filename))
            url_info = urllib.urlopen(url).info()
            remote_file_size = int(url_info.getheaders("Content-Length")[0])
            if existing_file_size == 0:
                INFO("Found existing file of size 0. Re-downloading...")
            elif existing_file_size == remote_file_size:
                INFO("Found existing file matching remote size. "
                     "Skipping download.")
                return filename
            else:
                INFO("Local file with size {0} does not match "
                     "remote file size {1}."
                     .format(existing_file_size, remote_file_size))

        (filename, _) = urllib.urlretrieve(url, filename,
                                                self._print_progress)

        return filename

class SourceForgeDownloader(BaseDownloader):
    FILESURL = "https://sourceforge.net/projects/{0.project}/files/{0.path}/"
    FILES_RE = 'latest version.*?title=./{0.path}([^/:]*)'
    DOWNLOADURL = "http://downloads.sourceforge.net/project/" \
                  "{0.project}/{0.path}/{0.filename}" \
                  "?use_mirror=autoselect&ts={0.ts}"

    project = None
    tool = None
    version = None

    @property
    def path(self):
        path = ""
        if self.tool:
            path += self.tool + "/"
        print(path)
        if self.version:
            path += self.version
        print(path)
        return path

    @property
    def filename(self):
        return "{0.tool}-{0.version}-{0.osname}-x86_64.zip" \
               .format(self)

    @property
    def osname(self):
        osname=os.uname()[0]
        if osname == "Darwin":
            return "macos"
        elif osname == "Linux":
            return "linux"
        else:
            return "unknown"

    @property
    def ts(self):
        return str(int(time.time()))

    def get_current_version(self):
        doc = self.fetch_url(self.FILESURL.format(self))
        match = re.search(self.FILES_RE.format(self),
                          doc.replace("\n",""))
        if match is None:
            raise DownloadException("couldn't find latest version of {}"
                                    .format(self.project + self.path))
        return match.groups()[0]

    def download_file(self):
        return super(SourceForgeDownloader, self).download_file(
            self.DOWNLOADURL.format(self),

        )

    def run(self):
        INFO("Retrieving version informatiom")
        self.version = self.get_current_version()
        INFO("Most recent version is: {}".format(self.version))

        return self.download_file()


class Bowtie2Downloader(SourceForgeDownloader):
    project = "bowtie-bio"
    tool = "bowtie2"

class BBMapDownloader(SourceForgeDownloader):
    DOWNLOADURL = "http://downloads.sourceforge.net/project/" \
                  "{0.project}/{0.path}" \
                  "?use_mirror=autoselect&ts={0.ts}"
    project = "bbmap"

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Valid arguments: bowtie2 bbmap")
        sys.exit()
    cmd = lower(sys.argv[1])
    if cmd == 'bowtie2':
        downloader = Bowtie2Downloader()
    elif cmd == 'bbmap':
        downloader = BBMapDownloader()
    else:
        print("can't download that")
        exit()

    downloader.run()
