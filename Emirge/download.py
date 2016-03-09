"""
helpers for downloading things
"""

import os
import sys
import urllib
import urllib2
from string import lower

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

        (filename, header) = urllib.urlretrieve(url, filename,
                                                self._print_progress)

        return filename

class SourceForgeDownloader(BaseDownloader):
    FILESURL = "https://sourceforge.net/projects/{project}/files/{path}"
    FILES_RE = 'latest version.*?title=./{path}([^/:]*)'
    DOWNLOADURL = "http://downloads.sourceforge.net/project/{project}/{" \
                  "path}?use_mirror=autoselect&ts={ts}"

    project = None
    tool = None

    def get_current_version(self):
        if path is not "" and path[-1] != "/":
            path += "/"
        doc = self.fetch_url(self.FILESURL.format(**self))
        match = re.search(self.FILES_RE.format(path),
                          doc.replace("\n",""))
        if match is None:
            raise DownloadException("couldn't find latest version of {}"
                                    .format(project + path))
        return match.groups()[0]

    def download_file(self):
        return super(SourceForgeDownloader, self).download_file(
            self.DOWNLOADURL.format(**self),

        )

    def run(self):
        self.version = self.get_current_version()
        return self.download_file()


class Bowtie2Downloader(SourceForgeDownloader):
    project = "bowtie-bio"
    tool = "bowtie2"

class BBMapDownloader(SourceForgeDownloader)
    project = "bbmap"

if __name__ == '__main__':
    downloader = None
    if lower(sys.argv[1]) == 'bowtie2':
        downloader = Bowtie2Downloader()

    if downloader is not None:
        downloader.run()
