import os
import sys
import urllib
import urllib2

from Emirge.log import INFO


class DownloadException(Exception):
    """Used to pass the URL downwards for urllib(2) exceptions"""
    pass


def fetch_url(url):
    """Fetches URL contents into variable

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


def print_progress(block, blocksize, total):
    """Callback function to print progress during download

    Args:
        block:     number of block downloaded
        blocksize: size of a block in bytes
        total:     total file size in bytes

    Returns:
    """

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
    """Downloads URL as file into folder. Returns filename.

    Skips download if the file already exists and has the same size.

    Args:
        url:    URL to download
        folder: target folder

    Returns:
        name of file

    """
    print "Downloading \"{0}\" to \"{1}\"".format(url, folder)

    filename = os.path.join(folder, url.split("/")[-1])
    if os.path.isfile(filename):
        existing_file_size = os.path.getsize(filename)
        urlinfo = urllib.urlopen(url).info()
        remote_file_size = urlinfo.getheaders("Content-Length")[0]
        if int(existing_file_size) == int(remote_file_size):
            INFO(
                "Found existing file matching remote size. Skipping download."
            )
            return filename
        else:
            INFO(
                "Local file with size {0} does not match remote file size {1}."
                .format(existing_file_size, remote_file_size)
            )

    (filename, header) = urllib.urlretrieve(url, filename, print_progress)

    return filename
