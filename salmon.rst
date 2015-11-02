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

Start up an m4.xlarge running blank Ubuntu 14.04.  (This gives you 15 GB of
RAM.)

Log in with MobaXterm or ssh.  (See `using Amazon docs
<http://angus.readthedocs.org/en/2015/amazon/>`__ for help.)

Install the necessary software::

   sudo apt-get update && \
   sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
         default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
         python-matplotlib python-pip python-virtualenv sysstat fastqc \
         trimmomatic bowtie samtools blast2 cmake libboost-all-dev liblzma-dev

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

Now, create an EBS volume from snapshot snap-a84c2ee7, attach it to
your machine, and mount it as /mnt/data.  Also make sure /mnt/ is
writeable::

   sudo mount /dev/xvdf /mnt/data
   sudo chmod a+rwxt /mnt

(This snapshot contains data from `Tulin et al., 2013
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
``head -400000`` line if you want all of them.)

Quantify them::

   for i in /mnt/data/*.pe.qc.fq.gz;
   do
      BASE=$(basename $i .pe.qc.fq.gz)
      salmon quant -i nema_index -1 $BASE.1.fq -2 $BASE.2.fq \
             -o $BASE.quant --libType IU;
   done

::

   find . -name \*.sf -exec grep "mapping rate"