import optparse, sys, pysam
from Emirge.rep_finder import RepFinder
from collections import Counter

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-k", metavar="K", type="int", default=16,
                      help="%default")
    parser.add_option("-l", metavar="length",type="int", default=50,
                      help="%default")
    parser.add_option("-v", action="store_true", default=False)
    (options, args) = parser.parse_args(sys.argv)

    f = RepFinder(k=options.k)

    matches_by_distance = Counter()
    matches_by_length = Counter()
    smatches_by_count = Counter()
    smatches_by_tax = Counter()
    taxsize = Counter()

    for filename in args:
        with pysam.FastxFile(filename) as fh:
            for i,entry in enumerate(fh):
                if not entry.comment: continue
                res = f.check(entry.sequence, minlen=options.l)
                tax = ";".join(entry.comment.split(";")[:-1])
                taxsize[tax] += 1
                if len(res):
                    print("{:>12}: {}".format(
                        entry.name,
                        " ".join([
                            "{},{}-{},{}-{}".format(l, s1, s1+l, s2, s2+l)
                            for l, s1, s2 in res
                            ])
                        ))
                    for l, s1, s2 in res:
                        if options.v:
                            print(entry.sequence[s1:s1+l])
                        matches_by_distance[s2 - l - s1] += 1
                        matches_by_length[l] += 1
                    smatches_by_count[len(res)] += 1
                    smatches_by_tax[tax] += 1

    print("\n"
          "Matches by Distance\n"
          "===================")
    vals = matches_by_distance.items()
    vals.sort()
    for elem, val in vals:
            if val>0:
                print("{:4}: {:4}".format(elem, val))

    print("\n"
          "Matches by Length\n"
          "=================")
    vals = matches_by_length.items()
    vals.sort()
    for elem, val in vals:
            if val>0:
                print("{:4}: {:4}".format(elem, val))

    print("\n"
          "Matches per Sequence\n"
          "====================")
    vals = smatches_by_count.items()
    vals.sort()
    for elem, val in vals:
            if val>0:
                print("{:4}: {:4}".format(elem, val))

    print("\n"
          "Matches by taxonomy\n"
          "===================")
    vals = smatches_by_tax.items()
    vals.sort()
    for elem, val in vals:
            if val>1:
                print("{:4}: {:4.1f}%".format(elem, 100.0*val/taxsize[elem]))


