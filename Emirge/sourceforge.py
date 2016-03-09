"""
download from sourceforge
"""
import json
import re

import time

import os

from Emirge.download import fetch_url, DownloadException, download_url
from Emirge.log import INFO



d

def list_files(project, path=""):
    js = fetch_url(BASEURL.format(project=project, path=path)
                     + "/list")
    js = json.loads(js)
    for line in js:
        INFO(line)

def download_file(project, path, folder=".", filename=None):
    file = download_url(
        "{base}{project}/{path}?use_mirror=autoselect&ts={ts}"
        .format(base=DOWNLOADURL,
                project=project,
                path=path,
                ts=str(time.time())),
        folder, filename=filename
    )
    # TODO: verify md5
    return file

def get_bowtie2(osname = None):
    if osname is None:
        osname = os.uname()[0]
    if osname == "Darwin":
        osname = "macos"
    elif osname == "Linux":
        osname = "linux"

    project = "bowtie-bio"
    tool = "bowtie2"
    version = get_version(project, tool)
    filename = "{tool}-{version}-{osname}-x86_64.zip"
    path = "{tool}/{version}/"
    return download_file(
        project,
        (path + filename).format(tool=tool, version=version, osname=osname),
        filename=filename.format(tool=tool, version=version, osname=osname)
    )

def get_bbmap():
    project="bbmap"
    version = get_version(project)
    return download_file(project, version, filename=version)

def get_bwa():
    pass
# http://downloads.sourceforge.net/project/
# bowtie-bio/bowtie2/2.2.7/bowtie2-2.2.7-linux-x86_64.zip

def test():
    INFO(get_bowtie2())

if __name__ == '__main__':
    test()
