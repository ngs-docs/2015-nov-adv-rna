Measuring transcriptome quality with transrate
==============================================

How do you measure the quality of your transcriptome? In some of the
beginner workshops, we suggested mapping your RNAseq reads back to
the transcriptome and counting the fraction that mapped.  Transrate
takes this kind of idea quite a bit further and measures several
read-based metrics.

Transrate Web site + docs: http://hibberdlab.com/transrate/

Transrate preprint: http://biorxiv.org/content/early/2015/06/27/021626

Starting up a machine
---------------------

Start up an m3.xlarge running blank Ubuntu 14.04.  (This gives you 15 GB of
RAM, plus lots of working disk space on /mnt.)

Log in with MobaXterm or ssh.  (See `using Amazon docs
<http://angus.readthedocs.org/en/2015/amazon/>`__ for help.)

Install the necessary software::

   sudo apt-get update && \
   sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
        python-matplotlib python-pip python-virtualenv sysstat fastqc \
        trimmomatic bowtie samtools blast2 cmake libboost-all-dev liblzma-dev \
        r-bioc-edgeR hmmer ncbi-blast+-legacy emboss

Install `khmer <http://khmer.readthedocs.org/en/v2.0/>`__::

   cd ~/
   python2.7 -m virtualenv work
   source work/bin/activate
   pip install -U setuptools
   git clone --branch v2.0 https://github.com/dib-lab/khmer.git
   cd khmer
   make install

Now, grab some test data::

   sudo chmod a+rwxt /mnt
   mkdir /mnt/data
   cd /mnt/data/
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/nema-subset.tar.gz
   tar xzf nema-subset.tar.gz

Install transrate
-----------------

Grab and install transrate::

   cd
   curl -O -L https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
   tar xzf transrate-1.0.1-linux-x86_64.tar.gz

   export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64
   echo 'export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64' >> ~/.bashrc
   export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64/bin
   echo 'export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64/bin' >> ~/.bashrc

   transrate --install-deps ref

Copy over the data
------------------

Create a working directory::

   mkdir /mnt/transrate
   cd /mnt/transrate

Copy in the data, fixing any long headers that might be leftover from
Trinity or Cufflinks or whatnot::

   sed -e '/^>.\{81\}/ s/^\(.\{80\}\).*$/\1/' /mnt/data/nema.fa > nema.fa

Running an initial evaluation: contig metrics
---------------------------------------------

Run an initial evaluation of your assembly, without using any reads or
reference transcriptome::

   transrate -a nema.fa

You should see the following output::

   [ INFO] 2015-11-01 15:13:17 : Contig metrics:
   [ INFO] 2015-11-01 15:13:17 : -----------------------------------
   [ INFO] 2015-11-01 15:13:17 : n seqs                       198151
   [ INFO] 2015-11-01 15:13:17 : smallest                        201
   [ INFO] 2015-11-01 15:13:17 : largest                       17655
   [ INFO] 2015-11-01 15:13:17 : n bases                   137744672
   [ INFO] 2015-11-01 15:13:17 : mean len                     695.15
   [ INFO] 2015-11-01 15:13:17 : n under 200                       0
   [ INFO] 2015-11-01 15:13:17 : n over 1k                     37271
   [ INFO] 2015-11-01 15:13:17 : n over 10k                       64
   [ INFO] 2015-11-01 15:13:17 : n with orf                    46134
   [ INFO] 2015-11-01 15:13:17 : mean orf percent              63.77
   [ INFO] 2015-11-01 15:13:17 : n90                             252
   [ INFO] 2015-11-01 15:13:17 : n70                             573
   [ INFO] 2015-11-01 15:13:17 : n50                            1315
   [ INFO] 2015-11-01 15:13:17 : n30                            2271
   [ INFO] 2015-11-01 15:13:17 : n10                            4111
   [ INFO] 2015-11-01 15:13:17 : gc                             0.44
   [ INFO] 2015-11-01 15:13:17 : gc skew                        0.01
   [ INFO] 2015-11-01 15:13:17 : at skew                         0.0
   [ INFO] 2015-11-01 15:13:17 : cpg ratio                      1.73
   [ INFO] 2015-11-01 15:13:17 : bases n                           0
   [ INFO] 2015-11-01 15:13:17 : proportion n                    0.0
   [ INFO] 2015-11-01 15:13:17 : linguistic complexity          0.13
   [ INFO] 2015-11-01 15:13:17 : Contig metrics done in 35 seconds

You can read more about the `contig metrics, here <http://hibberdlab.com/transrate/metrics.html#contig-metrics>`__.

Running a reference analysis: comparative metrics
-------------------------------------------------

Let's download the existing reference transcriptome and see how our
own assembled transcriptome compares::

   curl -O http://cnidarians.bu.edu/stellabase/assembly/NvT1.fasta
   transrate -a nema.fa --reference NvT1.fasta

(This will take about 20 minutes, note.)

Results::

   [ INFO] 2015-11-01 15:42:28 : Comparative metrics:
   [ INFO] 2015-11-01 15:42:28 : -----------------------------------
   [ INFO] 2015-11-01 15:42:28 : CRBB hits                    106203
   [ INFO] 2015-11-01 15:42:28 : n contigs with CRBB          106203
   [ INFO] 2015-11-01 15:42:28 : p contigs with CRBB            0.54
   [ INFO] 2015-11-01 15:42:28 : rbh per reference              0.92
   [ INFO] 2015-11-01 15:42:28 : n refs with CRBB              44743
   [ INFO] 2015-11-01 15:42:28 : p refs with CRBB               0.39
   [ INFO] 2015-11-01 15:42:28 : cov25                         19091
   [ INFO] 2015-11-01 15:42:28 : p cov25                        0.17
   [ INFO] 2015-11-01 15:42:28 : cov50                         13278
   [ INFO] 2015-11-01 15:42:28 : p cov50                        0.11
   [ INFO] 2015-11-01 15:42:28 : cov75                          8519
   [ INFO] 2015-11-01 15:42:28 : p cov75                        0.07
   [ INFO] 2015-11-01 15:42:28 : cov85                          6695
   [ INFO] 2015-11-01 15:42:28 : p cov85                        0.06
   [ INFO] 2015-11-01 15:42:28 : cov95                          4786
   [ INFO] 2015-11-01 15:42:28 : p cov95                        0.04
   [ INFO] 2015-11-01 15:42:28 : reference coverage             0.16
   [ INFO] 2015-11-01 15:42:28 : Comparative metrics done in 1377 seconds
   [ INFO] 2015-11-01 15:42:28 : -----------------------------------

You can read more about the `comparative metrics, here <http://hibberdlab.com/transrate/metrics.html#comparative-metrics>`__.

A really important note: this analysis can be done not only with a
DNA/RNA file of transcripts from your organism, but also with a
**peptide file from a nearby reference organism.**

Running a read-based analysis: read mapping metrics
---------------------------------------------------

The most powerful metrics that transrate offers are the `read mapping
metrics
<http://hibberdlab.com/transrate/metrics.html#read-mapping-metrics>`__.
These look at how the reads actually map to your transcriptome, and how
well the transcripts in your transcriptome are supported by the reads.

Next, let's evaluate against reads, prepared as in `salmon.rst <salmon.rst>`__::

   ln -fs ../data/*.?.fq .

   LIST1=$(ls -1 *.1.fq | head -1 | sort -n | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/')
   LIST2=$(ls -1 *.2.fq | head -1 | sort -n | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/')

   transrate -a nema.fa --left=$LIST1 --right=$LIST2

Results::

   [ INFO] 2015-11-03 16:24:39 : Read mapping metrics:
   [ INFO] 2015-11-03 16:24:39 : -----------------------------------
   [ INFO] 2015-11-03 16:24:39 : fragments                     50000
   [ INFO] 2015-11-03 16:24:39 : fragments mapped              46378
   [ INFO] 2015-11-03 16:24:39 : p fragments mapped             0.93
   [ INFO] 2015-11-03 16:24:39 : good mappings                 42113
   [ INFO] 2015-11-03 16:24:39 : p good mapping                 0.84
   [ INFO] 2015-11-03 16:24:39 : bad mappings                   4265
   [ INFO] 2015-11-03 16:24:39 : potential bridges               146
   [ INFO] 2015-11-03 16:24:39 : bases uncovered           133517286
   [ INFO] 2015-11-03 16:24:39 : p bases uncovered              0.97
   [ INFO] 2015-11-03 16:24:39 : contigs uncovbase            198150
   [ INFO] 2015-11-03 16:24:39 : p contigs uncovbase             1.0
   [ INFO] 2015-11-03 16:24:39 : contigs uncovered            197905
   [ INFO] 2015-11-03 16:24:39 : p contigs uncovered             1.0
   [ INFO] 2015-11-03 16:24:39 : contigs lowcovered           198133
   [ INFO] 2015-11-03 16:24:39 : p contigs lowcovered            1.0
   [ INFO] 2015-11-03 16:24:39 : contigs segmented                12
   [ INFO] 2015-11-03 16:24:39 : p contigs segmented             0.0
   [ INFO] 2015-11-03 16:24:39 : Read metrics done in 196 seconds
   [ INFO] 2015-11-03 16:24:39 : No reference provided, skipping comparative diagnostics
   [ INFO] 2015-11-03 16:24:39 : TRANSRATE ASSEMBLY SCORE      0.009
   [ INFO] 2015-11-03 16:24:39 : -----------------------------------
   [ INFO] 2015-11-03 16:24:39 : TRANSRATE OPTIMAL SCORE       0.184
   [ INFO] 2015-11-03 16:24:39 : TRANSRATE OPTIMAL CUTOFF     0.1432
   [ INFO] 2015-11-03 16:24:40 : good contigs                   3332
   [ INFO] 2015-11-03 16:24:40 : p good contigs                 0.02
   [ INFO] 2015-11-03 16:24:40 : Writing contig metrics for each contig to /mnt/transrate/transrate_results/nema/contigs.csv
   [ INFO] 2015-11-03 16:25:24 : Writing analysis results to assemblies.csv
   
Of particular note, this analysis may be the analysis you want to try before
deciding if you should generate a new transcriptome.
   
Challenge exercise
------------------
   
Repeat the above analyses with the transcriptome published in `Tulin et
al., 2013 <http://www.evodevojournal.com/content/4/1/16>`__::
   
   curl -L https://darchive.mblwhoilibrary.org/bitstream/handle/1912/5613/Trinity.fasta > tulin-2013-long.fa
      
You'll need to run the sed command, above, to convert
tulin-2013-long.fa into tulin-2013.fa.

Is the Tulin transcriptome better or worse than the more recently assembled
one (nema.fa, above)?

----

`Return to agenda <AGENDA.md>`__
