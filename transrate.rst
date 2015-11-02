Measuring transcriptome quality with transrate
==============================================

How do you measure the quality of your transcriptome? In some of the
beginner workshops, we suggested mapping your RNAseq reads back to
the transcriptome and counting the fraction that mapped.  Transrate
takes this kind of idea quite a bit further and measures several
read-based metrics.

Transrate Web site + docs: http://hibberdlab.com/transrate/

Transrate preprint: http://biorxiv.org/content/early/2015/06/27/021626

We'll start from the m4.xlarge Amazon machine booted & configured in
`salmon.rst <salmon.rst>`__.  If you are just running this, you'll need
to run the apt-get commands, install khmer, and mount the data snapshot
before continuing.

Install transrate
-----------------

Next, grab transrate::

   cd
   curl -O -L https://bintray.com/artifact/download/blahah/generic/transrate-1.0.1-linux-x86_64.tar.gz
   tar xzf transrate-1.0.1-linux-x86_64.tar.gz

   export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64
   echo 'export PATH=$PATH:$HOME/transrate-1.0.1-linux-x86_64' >> ~/.bashrc

Get the data
------------

Create a working directory::

   mkdir /mnt/transrate
   cd /mnt/transrate

Copy in the data, fixing any long headers that might be leftover from
Trinity or Cufflinks or whatnot::

   sed -e '/^>.{81}/ s/^(.{80}).*$/\1/' /mnt/data/nema.fa > nema.fa

Run an initial evaluation of your assembly, without using any reads::

   transrate -a nema.fa

See::

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

####

Then::

   curl -O http://cnidarians.bu.edu/stellabase/assembly/NvT1.fasta
   transrate -a nema.fa --reference NvT1.fasta

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

Grab the transcriptome from `Tulin et al., 2013
<http://www.evodevojournal.com/content/4/1/16>`__::

   curl -L https://darchive.mblwhoilibrary.org/bitstream/handle/1912/5613/Trinity.fasta > tulin-2013.fa

Next, let's evaluate against reads, prepared as in salmon.rst::

   ln -fs ../quant/*.?.fq

   LIST1=$(ls -1 *.1.fq | sort -n | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/')
   LIST2=$(ls -1 *.2.fq | sort -n | awk -vORS=, '{ print $1 }' | sed 's/,$/\n/')


Run! :::

   transrate -a nema.fa --left=$LIST1 --right=$LIST2

Challenge:

* evaluate reference transcriptome, published transcriptome
