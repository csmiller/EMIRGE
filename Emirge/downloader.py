"""
helpers for downloading things
"""

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


def print_progress(block, block_size, total):
    """Callback function to print progress during download

    Args:
        block:     number of block downloaded
        block_size: size of a block in bytes
        total:     total file size in bytes

    Returns:
    """

    blocks = int((total-1) / block_size) + 1
    line_width = min(77, blocks)
    if block == 0:
        print "Downloading file of size {0}:".format(total)
        if line_width > 1:
            print "|" + "-" * (line_width-2) + "|"
    else:
        if block % (blocks/line_width) == 0:
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
        url_info = urllib.urlopen(url).info()
        remote_file_size = url_info.getheaders("Content-Length")[0]
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
