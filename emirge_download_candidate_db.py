#!/usr/bin/env python2
import sys
import os
import urllib2
import email.utils
import datetime
import hashlib

urlbase     = 'https://googledrive.com/host/0B7hz7JVEE15dbUtkRmxKVlhtd1U'
db_filename = 'SSURef_111_candidate_db.fasta.gz'
md5sum      = '4f9a12fe76dd33cf8101ad66d079c63e'

url = '%s/%s'%(urlbase, db_filename)

try:
    local_file_utc_timestamp = os.path.getmtime(os.path.join(os.getcwd(), db_filename))
except OSError:  # no such file
    local_file_utc_timestamp = 1  # seconds from epoch.  Make this really old file, so url file always newer.

try:
    r = urllib2.urlopen(url)
    url_last_modified_string = r.info().dict['last-modified']
    url_last_modified_utc_timestamp = email.utils.mktime_tz(email.utils.parsedate_tz(url_last_modified_string))

    if local_file_utc_timestamp < url_last_modified_utc_timestamp:  # url file newer
        print "downloading default candidate db from %s... "%(url)
        with open(db_filename, "wb") as local_file:
            local_file.write(r.read())
    else:
        print "NOTE: Skipping download of default candidate db because up-to-date file %s already found in local directory."%db_filename

    h = hashlib.md5(open(db_filename, "rb").read()).hexdigest()
    if h != md5sum:
        print >> sys.stderr, "*** WARNING: Local file %s appears to be corrupted!\n*** md5sum   = %s\n*** expected = %s\n*** Try to manually download from\n*** %s"%(db_filename, h, md5sum, url)
    
except urllib2.HTTPError, e:
    print >> sys.stderr,  "*** ERROR, candidate db file not downloaded! Try to manually download from\n%s\n\nHTTP Error:"%url, e.code, url
except urllib2.URLError, e:
    print >> sys.stderr,  "*** ERROR, candidate db file not downloaded! Try to manually download from\n%s\n\nURL Error:"%url, e.reason, url
except:
    print >> sys.stderr,  "*** ERROR, candidate db file not downloaded! Try to manually download from\n%s"%url
    
