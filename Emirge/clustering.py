import re
from subprocess import check_output, PIPE, Popen

from Emirge.io import command_avail, check_call
from Emirge.log import CRITICAL


class vsearch(object):
    binary = "vsearch"
    @classmethod
    def available(cls):
        """
        Check if vsearch is present in correct version
        """
        if not command_avail(cls.binary):
            return False

        out = Popen([cls.binary, "--version"], stdout = PIPE, stderr = PIPE)
        match = re.search('vsearch.* v([0-9]*)\.([0-9]*)\.([0-9]*)', out.communicate()[1])
        version = tuple(int(x) for x in match.groups())

        if version < (1,1,0):
            CRITICAL(
                "FATAL: vsearch version found was {}.{}.{}. "
                "EMIRGE requires vsearch version 1.1.0 or above."
                "".format(*version)
            )
            return False

        return True

    def search(self, query, db, maxaccepts=0, maxrejects=0, id=1):
        out = fasta + ".us.txt"
        cmd = [
            self.binary,
            "--usearch_global", query,
            "--db", db,
            "--query_cov", "0.5",
            "--strand", "plus",
            "--userfields", "query+target+id+caln+qlo+qhi+tlo+thi",
            "--quiet",
            "--threads", str(self.n_cpus),
            "--maxaccepts", str(maxaccepts),
            "--maxrejects", str(maxrejects),
            "--id", str(id),
            "--userout", out,
        ]
        check_call(cmd)
        return out

    def cluster(self, fasta, id):
        cmd = [
            self.binary,
            "--cluster_smallmem", fasta,
            "--usersort",
            "--notrunclabels",
            "--id", str(id),
            "--centroids", centroids,
            "--uc", uc
        ]
