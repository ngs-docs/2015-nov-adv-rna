Using Salmon Quantification to Reduce Transcriptome Complexity
==============================================================

In many assemblies - especially de novo assemblies, done without a
reference genome - there are many transcripts that don't have significant
expression levels.  This can interfere with quantification and annotation.
We can use Salmon to eliminate many of these transcripts, by estimating
their transcription levels and removing low-expressing transcripts.

----

We'll start from the m3.xlarge Amazon machine booted & configured in
`salmon.rst <salmon.rst>`__.  If you are just running this, you'll need
to run the apt-get commands, install khmer, and mount the data snapshot
before continuing.

First, create a working directory and download a full set of Salmon
quant results (done on all the reads)::

   mkdir /mnt/reduced
   cd /mnt/reduced

   curl -O http://athyra.idyll.org/~t/nema-salmon-quant.tar.gz
   tar xzf nema-salmon-quant.tar.gz

Reducing transcriptome numbers
------------------------------

Reduce transcriptome complexity by removing low coverage contigs. Use awk to sum the counts in the 
numReads column of each quant.sf file and print contignames that have counts above threshold 
(10 total counts)::

   awk '{x[$1] += $4} END {for(y in x) if (x[y] >10) print y}' *.quant/quant.sf | sort > contigs_over_threshold.txt

Then, use this list of contigs to subset the fasta reference, by getting (and running) a python script to subset the fasta file, keeping all entries that match the contig names above::

   curl -O https://raw.githubusercontent.com/bluegenes/fasta_tools/master/extract_fasta.py
   python extract_fasta.py /mnt/data/nema.fa contigs_over_threshold.txt

You'll see the results here::

   ls -la /mnt/data/nema*

Evaluate the full and reduced transcriptomes using BUSCO.
---------------------------------------------------------

Install BUSCO::

   cd $HOME
   curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
   tar -zxf BUSCO_v1.1b1.tar.gz
   cd BUSCO_v1.1b1/
   chmod +x BUSCO_v1.1b1.py
   export PATH=$PATH:$(pwd)

Run BUSCO on the filtered assembly::

   mkdir /mnt/busco
   cd /mnt/busco

   #Download metazoa busco database
   curl -LO http://busco.ezlab.org/files/metazoa_buscos.tar.gz
   tar -zxf metazoa_buscos.tar.gz

   python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans -in /mnt/data/nema_inList.fa \
      --cpu 16 -l metazoa -f -o nema_inList_metazoaBusco

   python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -m trans -in /mnt/data/nema.fa \
      --cpu 16 -l metazoa -f -o nema_all_metazoaBusco

   less run*/short*


Visualize transcriptome coverage in R::

   curl -O -L https://raw.githubusercontent.com/ngs-docs/2015-nov-adv-rna/master/plotTPM.R
   Rscript plotTPM.R

