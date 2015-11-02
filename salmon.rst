Quantifying RNAseq data with salmon
===================================

Salmon is one of a breed of new, very fast RNAseq counting packages.
Like Kallisto and Sailfish, Salmon counts fragments without doing
up-front read mapping.  Salmon can be used with edgeR and others
to do differential expression analysis.

Salmon preprint: http://biorxiv.org/content/early/2015/06/27/021592

Salmon documentation: http://salmon.readthedocs.org/en/latest/

Starting up a machine
---------------------

We're going to use Amazon Web Services for this.

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

Clone salmon from github (use CTB's version that's fixed for Boost 1.5.4)::

   cd
   git clone https://github.com/ctb/salmon.git -b boost1.54
   cd salmon
   cmake .
   make
   make install

Put Salmon in your path::

   export PATH=$PATH:$HOME/salmon/bin
   echo 'export PATH=$PATH:$HOME/salmon/bin' >> ~/.bashrc

   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/salmon/lib
   echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/salmon/lib' >> ~/.bashrc

Getting the data
----------------

.. Now, create an EBS volume from snapshot snap-a84c2ee7, attach it to
   your machine, and mount it as /mnt/data.  Also make sure /mnt/ is
   writeable::

      sudo mkdir /mnt/data
      sudo mount /dev/xvdf /mnt/data
      sudo chmod a+rwxt /mnt

Do::

   sudo chmod a+rwxt /mnt
   mkdir /mnt/data
   cd /mnt/data/
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/nema-subset.tar.gz
   tar xzf nema-subset.tar.gz

(This is data from `Tulin et al., 2013
<http://www.evodevojournal.com/content/4/1/16>`__ that was processed
and assembled with `the khmer protocols steps 1-3
<http://khmer-protocols.readthedocs.org/en/ctb/mrnaseq/index.html>`__.)

Make a directory to work in::

   mkdir /mnt/quant

Copy in the transcriptome from the snapshot::

   cd /mnt/quant
   cp /mnt/data/nema.fa .

Index it with salmon::

   salmon index --index nema_index --transcripts nema.fa --type quasi   

Grab the reads and split them up::

   for i in /mnt/data/*.pe.qc.fq.gz;
   do
      BASE=$(basename $i .pe.qc.fq.gz);
      zcat $i |
           head -400000 |
           split-paired-reads.py -1 $BASE.1.fq -2 $BASE.2.fq;
   done

(Note, here we're taking only the first 100,000 reads; remove the
``head -400000`` line if you want all of them, but that will take much
longer.)

Now, quantify the reads against the reference using Salmon::

   for i in /mnt/data/*.pe.qc.fq.gz;
   do
      BASE=$(basename $i .pe.qc.fq.gz)
      salmon quant -i nema_index -1 $BASE.1.fq -2 $BASE.2.fq \
             -o $BASE.quant --libType IU;
   done

This will create a bunch of directories named something like
``0Hour_ATCACG_L002001.quant``, containing a bunch of files.  Take a look
at what files there are::

   find 0Hour_ATCACG_L002001.quant -type f

The two most interesting files are ``salmon_quant.log`` and ``quant.sf``.
The latter contains the counts; the former contains the log information
from running things.  Take a look at the log -- ::

   more 0Hour_ATCACG_L002001.quant/libFormatCounts.txt

and see what you think it means... (Use 'q' to quit out of more.)

One last thing before we move on to quantification -- the ``quant.sf`` files
contain the mapping rates per library, which might be of interest...::

   find . -name \*.sf -exec grep -H "mapping rate" {} \;

Working with the counts
-----------------------

Now, the ``quant.sf`` files actually contain the relevant information about
expression -- take a look::

   head -20 0Hour_ATCACG_L002001.quant/quant.sf

The first column contains the transcript names, and the
fourth column is what edgeR etc will want - the "raw counts".
However, they're not in a convenient location / format for edgeR to use;
let's fix that.

Download the ``gather-counts.py`` script::

   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/files/gather-counts.py

and run it::

   python ./gather-counts.py

This will give you a bunch of .counts files, processed from the quant.sf files
and named for the directory they are in.

Now, run an edgeR script (`nema.salmon.R
<https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/files/nema.salmon.R>`__)
that loads all this in and calculates a few plots -- ::

   curl -O -L https://raw.githubusercontent.com/ngs-docs/2015-nov-adv-rna/master/files/nema.salmon.R
   Rscript nema.salmon.R

These will produce two plots, nema-edgeR-MDS.pdf and nema-edgeR-MA-plot.pdf.
Try downloading them to your computer using either MobaXTerm or CyberDuck.

----

You can see the plot outputs for the whole data set (all the reads) here:

* `nema-edgeR-MDS.pdf <https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/files/nema-edgeR-MDS.pdf>`__
* `nema-edgeR-MA-plot.pdf <https://github.com/ngs-docs/2015-nov-adv-rna/blob/master/files/nema-edgeR-MA-plot.pdf>`__ (0 vs 6 hour)

A challenge exercise
--------------------

Download the entire counts data set::

  mkdir /mnt/fullquant
  cd /mnt/fullquant
  curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/files/nema-counts.tar.gz
  tar xzf nema-counts.tar.gz

and run edgeR differential expression etc on it, as above.

Then, create an MA plot comparing 6 Hour vs 12 Hour.

----

`Return to agenda <AGENDA.md>`__
