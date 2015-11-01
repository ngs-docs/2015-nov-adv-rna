Using Salmon Quantification to Reduce Transcriptome Complexity
===================================

Reduce transcriptome complexity by removing low coverage contigs. Use awk to sum the counts in the 
numReads column of each quant.sf file and print contignames, summed counts that are above threshold 
(10 total counts)::

   awk '{x[$1] += $4} END {for(y in x) if (x[y] >10) print y,x[y]}' *.sf > counts_over_threshold.txt

Once we're satified this has the right functionality, we can just print the contignames::

   awk '{x[$1] += $4} END {for(y in x) if (x[y] >10) print y}' *.sf | sort > contigs_over_threshold.txt


Use this list of contigs to subset the fasta reference::

   extract_fasta.py nema.fa contigs_over_threshold.txt nema_over_10.fa nema_under_10.fa


Evaluate the full and reduced transcriptomes using BUSCO.

Install BUSCO::

   cd $HOME
   curl -O http://busco.ezlab.org/files/BUSCO_v1.1b1.tar.gz
   tar -zxf BUSCO_v1.1b1.tar.gz
   cd BUSCO_v1.1b1/
   PATH=$PATH:$(pwd)

Run BUSCO for assemblies: There are Eukaryote, Metazoa, Arthropod, Vertebrate, Plant 
references for use with other genomes::

   mkdir /mnt/busco
   cd /mnt/busco

   #Download busco database
   tmux new -s busco

   curl -LO http://busco.ezlab.org/files/vertebrata_buscos.tar.gz
   tar -zxf vertebrata_buscos.tar.gz

   python3 /home/ubuntu/BUSCO_v1.1b1/BUSCO_v1.1b1.py \
   -m trans -in /mnt/assembly/trinity_out_dir/Trinity.fasta \
   --cpu 16 -l vertebrata -o trin.assemblty

   less run*/short*

   Control-b d #to exit tmux





