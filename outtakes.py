# bowtie2 version:
def do_initial_mapping_bt2(em, working_dir, options):
    """
#    IN:  takes the working directory and an OptionParser options object
#
#    does the initial 1-reference-per-read bowtie mapping to initialize the algorithm
#    OUT:  path to the bam file from this initial mapping
#    """
    initial_mapping_dir = os.path.join(working_dir, "initial_mapping")
    if not os.path.exists(initial_mapping_dir):
        os.mkdir(initial_mapping_dir)

    minins = max((options.insert_mean - 3 * options.insert_stddev), options.max_read_length)
    maxins = options.insert_mean + 3 * options.insert_stddev
    bampath_prefix = os.path.join(initial_mapping_dir, "initial_bowtie_mapping.PE")

    nicestring = ""
    if options.nice_mapping is not None:
        nicestring = "nice -n %d"%(options.nice_mapping)  # TODO: fix this so
        # it isn't such a hack and will work in bash.  Need to rewrite all
        #  subprocess code, really (shell=False)
    reads_ascii_offset = {False: 64, True: 33}[options.phred33]
    if options.fastq_reads_1.endswith(".gz"):
        option_strings = ["gzip -dc "]
    else:
        option_strings = ["cat "]
    # shared regardless of whether paired mapping or not
    option_strings.extend([options.fastq_reads_1, nicestring,
                           reads_ascii_offset, options.processors])

    # PAIRED END MAPPING
    if options.fastq_reads_2 is not None:
        option_strings.extend([minins, maxins, options.bowtie2_db,
                               options.fastq_reads_2, bampath_prefix])
        cmd = "%s %s | %s bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --no-unal" \
              " --phred%d -t -p %s -I %d -X %d --no-mixed --no-discordant " \
              "-x %s -1 - -2 %s | samtools view -b -S -u - > %s.u.bam "\
              %tuple(option_strings)
    # SINGLE END MAPPING
    else:
        option_strings.extend([options.bowtie2_db, bampath_prefix])
        cmd = "%s %s | %s bowtie2 -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --no-unal" \
              " --phred%d -t -p %s -x %s -U - | samtools view -b -S -u - > " \
              "%s.u.bam "%tuple(option_strings)


    log.warning("Performing initial mapping with command:\n%s" % cmd)
    p = Popen(cmd, shell=True, stdout=sys.stdout, stderr=PIPE, close_fds=True)
    p.wait()
    stderr_string = p.stderr.read()
    em.fragments_mapped = em.get_frags_mapped_bt2(stderr_string)
    # re-print this to stdout, since we stole it.
    sys.stdout.write(stderr_string)
    sys.stdout.flush()

    return bampath_prefix + ".u.bam"


    # original bowtie1 version:
def do_initial_mapping(em, working_dir, options):

    """
    IN:  takes the em object, working directory and an OptionParser options object

    does the initial 1-reference-per-read bowtie mapping to initialize the algorithm
    OUT:  path to the bam file from this initial mapping

    TODO:  move this to em method.  A bit clumsy right now.
    """
    initial_mapping_dir = os.path.join(working_dir, "initial_mapping")
    if not os.path.exists(initial_mapping_dir):
        os.mkdir(initial_mapping_dir)

    # minins = max((options.insert_mean - 3*options.insert_stddev), options.max_read_length) # see comments above.  This assumes there might be adapter sequence in dovetailed reads
    minins = max((options.insert_mean - 3*options.insert_stddev), 0) # Puts burden on user to make sure there are no adapter sequences in input reads.
    maxins = options.insert_mean + 3*options.insert_stddev
    bampath_prefix = os.path.join(initial_mapping_dir, "initial_bowtie_mapping.PE")

    nicestring = ""
    if options.nice_mapping is not None:
        nicestring = "nice -n %d"%(options.nice_mapping)  # TODO: fix this so it isn't such a hack and will work in non-bash shells.  Need to rewrite all subprocess code, really (shell=False)
    reads_ascii_offset = {False: 64, True: 33}[options.phred33]
    if options.fastq_reads_1.endswith(".gz"):
        option_strings = ["gzip -dc "]
    else:
        option_strings = ["cat "]
    # shared regardless of whether paired mapping or not
    option_strings.extend([options.fastq_reads_1, nicestring, reads_ascii_offset, options.processors, BOWTIE_l, BOWTIE_e])
    samtools_cmd    = "samtools view -S -h -u -b -F 0x0004 -"  # -F instead of piping to awk?    |  awk '{if ($3!="*") print }'

    # PAIRED END MAPPING
    if options.fastq_reads_2 is not None:
        option_strings.extend([minins, maxins, options.bowtie_db, options.fastq_reads_2, samtools_cmd, bampath_prefix])
        cmd = """%s %s | %s bowtie --phred%d-quals -t -p %s -n 3 -l %s -e %s --best --sam --chunkmbs 512 --minins %s --maxins %s %s -1 - -2 %s | %s > %s.u.bam """%tuple(option_strings)
    # SINGLE END MAPPING
    else:
        option_strings.extend([options.bowtie_db, samtools_cmd, bampath_prefix])
        cmd = """%s %s | %s bowtie --phred%d-quals -t -p %s -n 3 -l %s -e %s --best --sam --chunkmbs 512  %s - | %s > %s.u.bam """%tuple(option_strings)


    log.warning("Performing initial mapping with command:\n%s" % cmd)
    p = Popen(cmd, shell=True, stdout=sys.stdout, stderr=PIPE, close_fds=True)
    p.wait()
    stderr_string = p.stderr.read()
    em.fragments_mapped = em.get_frags_mapped_bt1(stderr_string)
    # re-print this to stdout, since we stole it.
    sys.stdout.write(stderr_string)
    sys.stdout.flush()

    return bampath_prefix + ".u.bam"



def do_premapping(pre_mapping_dir, options):
    """
    Do only if metagenomic reads
    IN: the metagenomic reads
    OUT: a new set of reads to be used in the iterations
    """

    if not os.path.exists(pre_mapping_dir):
        os.mkdir(pre_mapping_dir)

    nicestring = ""
    if options.nice_mapping is not None:
        nicestring = "nice -n %d"%(options.nice_mapping)  # TODO: fix this so it isn't such a hack and will work in non-bash shells.  Need to rewrite all subprocess code, really (s$

    reads_ascii_offset = {False: 64, True: 33}[options.phred33]

    premapSam = os.path.join(pre_mapping_dir, "bowtie_pre_mapping.PE.sam")
    premap_reads_1 = os.path.join(pre_mapping_dir, "pre_mapped_reads.1.fastq")
    premap_reads_2 = os.path.join(pre_mapping_dir, "pre_mapped_reads.2.fastq")

    if options.fastq_reads_1.endswith(".gz"):
        option_strings = ["gzip -dc "]
    else:
        option_strings = ["cat "]


    # premapping done as single reads regardless of whether paired mapping or not  CAN'T DEAL W/GZIPPED READ 2 HERE
    if options.fastq_reads_2 is not None:
        option_strings.extend([options.fastq_reads_1,options.fastq_reads_2,
                               nicestring, reads_ascii_offset,
                               options.processors,options.bowtie2_db,
                               premapSam])
        cmd = """%s %s %s | %s bowtie2 --very-sensitive-local --phred%d -t -p %s -k1 -x %s --no-unal -S %s -U - """ %tuple(option_strings)

    log.warning("Performing pre mapping with command:\n%s" % cmd)
    p = Popen(cmd, shell=True, stdout=sys.stdout, stderr=PIPE, close_fds=True)
    p.wait()
    stderr_string = p.stderr.read()
    sys.stdout.write(stderr_string)
    sys.stdout.flush()


    # get IDs of mapped reads
    sam = pysam.AlignmentFile(premapSam, 'r')
    keepers = set()
    for read in sam.fetch():
        keepers.add(read.query_name)

    ## use IDs to create new fastq files of reads where at least one of the pair mapped in the premapping step
    r1_in = (options.fastq_reads_1)
    r2_in = (options.fastq_reads_2)
    p1_out = open(premap_reads_1, 'w')
    p2_out = open(premap_reads_2, 'w')

    # this code from fastq_pull_sequence:
    for infile, outfile in [(r1_in, p1_out), (r2_in, p2_out)]:
        ks = Kseq(infile)
        keptseqs = 0
        total_seqs = 0
        while 1:
            t = ks.read_sequence_and_quals()  # (name, sequence, qualstring)
            if t is None:
                break

            # get name up to first whitespace, check if in file
            if t[0].split()[0] in keepers:
                outfile.write("@%s\n%s\n+\n%s\n" % (t[0], t[1], t[2]))
                keptseqs += 1
            total_seqs += 1

        log.warning("Total records read: %s" % (total_seqs))
        log.warning("Total premapped records written: %s" % (keptseqs))

@log.timed("Mapping reads for iteration {self.iteration_i}")
    def do_mapping(self, full_fasta_path, nice=None):
        """
        IN:  path of fasta file to map reads to
        run external mapping program to produce bam file
        right now this is bowtie

        should also set self.n_alignments and self.current_bam_filename
        """

        self.do_mapping_bowtie2(full_fasta_path, nice=nice)
        # self.do_mapping_bowtie(full_fasta_path, nice = nice)

    # bowtie2 version:
    def do_mapping_bowtie2(self, full_fasta_path, nice=None):
        """
        run bowtie2 to produce bam file for next iteration

        """
        bowtie_logfile = os.path.join(self.iterdir, "bowtie.iter.%02d.log"%(self.iteration_i))
        bowtie_index   = os.path.join(self.iterdir, "bowtie.index.iter.%02d"%(self.iteration_i))
        # 1. build index
        cmd = "bowtie2-build -o 3 %s %s > %s 2>&1"%(full_fasta_path , bowtie_index, bowtie_logfile) # -o 3 for speed? magnitude of speedup untested!
        log.info("\tbowtie2-build command:")
        log.info("\t%s" % cmd)
        check_call(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)

        # 2. run bowtie
        nicestring = ""
        if nice is not None:
            nicestring = "nice -n %d" % (nice)

        if self.reads1_filepath.endswith(".gz"):
            cat_cmd = "gzip -dc "
        else:
            cat_cmd = "cat "


        # minins = max((self.insert_mean - 3*self.insert_sd), self.max_read_length) # this was to keep "dovetailing" (overlapping reads that extend past each other) reads from mapping, but causes problems with modern overlapped reads common on MiSeq, etc.  Goal was to keep adapter sequence bases out of EMIRGE's view.
        minins = max((self.insert_mean - 3*self.insert_sd), 0) # Puts burden on user to make sure there are no adapter sequences in input reads.
        maxins = self.insert_mean + 3*self.insert_sd
        output_prefix = os.path.join(self.iterdir, "bowtie.iter.%02d"%(self.iteration_i))
        output_filename = "%s.PE.u.bam"%output_prefix


        # these are used for single reads too.  Bowtie2 ignores the paired-end options when single reads are used

        shared_bowtie_params = "-D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -k 20 " \
                               "--no-unal --phred%d -t -p %s" \
                               % (self.reads_ascii_offset, self.n_cpus)

        #Build Bowtie2 command depending if reads2 was given in Emirge command line parameters

        if self.reads2_filepath is not None:
            bowtie_command = "%s %s | %s bowtie2 %s -I %d -X %d --no-mixed " \
                             "--no-discordant -x %s -1 - -2 %s | " \
                             "samtools view -b -S - > %s 2> %s " \
                             % (cat_cmd,
                                self.reads1_filepath,
                                nicestring,
                                shared_bowtie_params,
                                minins, maxins,
                                bowtie_index,
                                self.reads2_filepath,
                                output_filename,
                                bowtie_logfile)
        else: # single reads
            bowtie_command = "%s %s | %s bowtie2 %s -x %s -U - | " \
                             "samtools view -b -S - > %s 2> %s " \
                             %(cat_cmd,
                               self.reads1_filepath,
                               nicestring,
                               shared_bowtie_params,
                               bowtie_index,
                               output_filename,
                               bowtie_logfile)

        log.info("\tbowtie command:")
        log.info("\t%s" % bowtie_command)

        p = Popen(bowtie_command, shell=True, stdout=sys.stdout, stderr=PIPE, close_fds=True)
        p.wait()
        stderr_string = p.stderr.read()
        self.fragments_mapped = self.get_frags_mapped_bt2(stderr_string)
        sys.stdout.write(stderr_string)
        sys.stdout.flush()

        log.info("\tFinished Bowtie for iteration %02d" % (self.iteration_i))

        # 3. clean up
        # check_call("samtools index %s.sort.PE.bam"%(output_prefix), shell=True, stdout = sys.stdout, stderr = sys.stderr)
        check_call("gzip -f %s" % (bowtie_logfile), shell=True)

        assert self.iterdir != '/'
        for filename in os.listdir(self.iterdir):
            assert (len(os.path.basename(bowtie_index)) >= 20)  # weak check that I'm not doing anything dumb.
            if os.path.basename(bowtie_index) in filename:
                os.remove(os.path.join(self.iterdir, filename))

        self.current_bam_filename = output_filename  # do this last.


   # original bowtie1 version:
    def do_mapping_bowtie(self, full_fasta_path, nice=None):
        """
        run bowtie to produce bam file for next iteration

        sets self.n_alignments
        sets self.current_bam_filename
        """
        bowtie_index   = os.path.join(self.iterdir, "bowtie.index.iter.%02d"%(self.iteration_i))
        bowtie_logfile = os.path.join(self.iterdir, "bowtie.iter.%02d.log"%(self.iteration_i))
        # 1. build index
        cmd = "bowtie-build -o 3 %s %s > %s"%(full_fasta_path , bowtie_index, bowtie_logfile) # -o 3 for speed? magnitude of speedup untested!
        # note: just send stdout to log file, as stderr captured in emirge stderr
        log.info("\tbowtie-build command:")
        log.info("\t%s" % cmd)
        check_call(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
        sys.stdout.flush()
        sys.stderr.flush()

        # 2. run bowtie
        nicestring = ""
        if nice is not None:
            nicestring = "nice -n %d" % (nice)

        if self.reads1_filepath.endswith(".gz"):
            cat_cmd = "gzip -dc "
        else:
            cat_cmd = "cat "

        # these are used for single reads too.
        shared_bowtie_params = "--phred%d-quals -t -p %s  -n 3 -l %s -e %s  --best --strata --all --sam --chunkmbs 512"%(self.reads_ascii_offset, self.n_cpus, BOWTIE_l, BOWTIE_e)

        # minins = max((self.insert_mean - 3*self.insert_sd), self.max_read_length) # this was to keep "dovetailing" (overlapping reads that extend past each other) reads from mapping, but causes problems with modern overlapped reads common on MiSeq, etc.  Goal was to keep adapter sequence bases out of EMIRGE's view.
        minins = max((self.insert_mean - 3*self.insert_sd), 0) # Puts burden on user to make sure there are no adapter sequences in input reads.
        maxins = self.insert_mean + 3*self.insert_sd
        output_prefix = os.path.join(self.iterdir, "bowtie.iter.%02d"%(self.iteration_i))
        output_filename = "%s.PE.u.bam"%output_prefix
        samtools_cmd    = "samtools view -S -h -u -b -F 0x0004 -"  # -F instead of piping to awk?    |  awk '{if ($3!="*") print }'

        if self.reads2_filepath is not None:
            bowtie_command = "%s %s | " \
                             "%s bowtie %s --minins %d --maxins %d %s " \
                             "-1 - -2 %s | %s > %s" \
                             %(cat_cmd,
                               self.reads1_filepath,
                               nicestring,
                               shared_bowtie_params,
                               minins, maxins,
                               bowtie_index,
                               self.reads2_filepath,
                               samtools_cmd,
                               output_filename)
        else:  # single reads
            bowtie_command = "%s %s | %s bowtie %s %s - | %s > %s" \
                             % (cat_cmd,
                                self.reads1_filepath,
                                nicestring,
                                shared_bowtie_params,
                                bowtie_index,
                                samtools_cmd,
                                output_filename)

        log.info("\tbowtie command:")
        log.info("\t%s" % bowtie_command)

        p = Popen(bowtie_command, shell=True, stdout=sys.stdout, stderr=PIPE, close_fds=True)
        p.wait()
        stderr_string = p.stderr.read()
        self.fragments_mapped = self.get_frags_mapped_bt1(stderr_string)


        # re-print this to stdout, since we stole it from bowtie
        sys.stdout.write(stderr_string)
        sys.stdout.flush()
        # and now put in separate bowtie logfile
        of = open(bowtie_logfile, 'w')
        of.write("\nBOWTIE STDERR:\n")
        of.write(stderr_string)
        of.write("\n")
        of.close()

        log.info("\tFinished Bowtie for iteration %02d" % (self.iteration_i))

        # 3. clean up
        # check_call("samtools index %s.sort.PE.bam"%(output_prefix), shell=True, stdout = sys.stdout, stderr = sys.stderr)
        if os.path.exists(bowtie_logfile):
            check_call("gzip -f %s" % (bowtie_logfile), shell=True)

        assert self.iterdir != '/'
        for filename in os.listdir(self.iterdir):
            assert (len(os.path.basename(bowtie_index)) >= 20)  # weak check that I'm not doing anything dumb.
            if os.path.basename(bowtie_index) in filename:
                os.remove(os.path.join(self.iterdir, filename))

        self.current_bam_filename = output_filename  # do this last.

 def get_frags_mapped_bt1(self, stderr_string):
        r = ""
        try:
            r = re.findall(r'Reported ([0-9]+) (paired-end )?alignments', stderr_string)
            if r[0][1] != '':  # "paired-end" string matched -- two lines in samfile per paired-end aln
                return int(r[0][0]) * 2
            else:  # single-end -- one line in samfile per alignment
                return int(r[0][0])
        except IndexError:
            log.error("OOPS, we didn't get number of reads from bowtie:")
            log.error(stderr_string)
            log.error(r)
            raise

    def get_frags_mapped_bt2(self, stderr_string):
        """
        49227 reads; of these:
            49227 (100.00%) were paired; of these:
                0 (0.00%) aligned concordantly 0 times
                36694 (74.54%) aligned concordantly exactly 1 time
                12533 (25.46%) aligned concordantly >1 times
        100.00% overall alignment rate
"""
        r = ""
        try:
            r = re.findall(r'([0-9]+) reads;', stderr_string)
            rt = re.findall(r'([0-9]+\.[0-9]+)% overall', stderr_string)
            frags = (float(rt[0]) / 100) * int(r[0])
            return int(frags)
        except IndexError:
            log.error("OOPS, we didn't get number of reads from bowtie:")
            log.error(stderr_string)
            log.error(r)
            rais