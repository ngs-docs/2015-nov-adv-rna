Using Salmon Quantification to Reduce Transcriptome Complexity
===================================

Reduce transcriptome complexity by removing low coverage contigs. Use awk to sum the counts in the 
numReads column of each quant.sf file and print contignames that have counts above threshold 
(10 total counts)::

   awk '{x[$1] += $4} END {for(y in x) if (x[y] >10) print y}' /mnt/data/*.quant/quant.sf | sort > /mnt/contigs_over_threshold.txt


Use this list of contigs to subset the fasta reference.

Get (and run) a python script to subset the fasta file, keeping all entries that match the contig names above::

   curl -O https://raw.githubusercontent.com/bluegenes/fasta_tools/master/extract_fasta.py
   python extract_fasta.py /mnt/nema.fa /mnt/contigs_over_threshold.txt


Evaluate the full and reduced transcriptomes using BUSCO.

Install BUSCO::

   cd $HOME
   curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
   tar -zxf BUSCO_v1.1b1.tar.gz
   cd BUSCO_v1.1b1/
   PATH=$PATH:$(pwd)

Run BUSCO on the filtered assembly::

   mkdir /mnt/busco
   cd /mnt/busco

   #Download metazoa busco database
   tmux new -s busco

   curl -LO http://busco.ezlab.org/files/metazoa_buscos.tar.gz
   tar -zxf metazoa_buscos.tar.gz

   python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
   -m trans -in /mnt/nema_inList.fa \
   --cpu 16 -l metazoa -o nema_inList_metazoaBusco

   less run*/short*

   Control-b d #to exit tmux



Visualize transcriptome coverage in R::

   Rscript plotTPM.R

