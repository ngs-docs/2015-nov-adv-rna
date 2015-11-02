Blasting with DIAMOND 
===================================

Make a new directory for diamond::
   
   cd
   mkdir diamond

Download diamond::
   
   wget http://github.com/bbuchfink/diamond/releases/download/v0.7.9/diamond-linux64.tar.gz
   tar xzf diamond-linux64.tar.gz

Put diamond in your path:
   
   export PATH=$PATH:$HOME/diamond
   echo 'export PATH=$PATH:$HOME/diamond' >> ~/.bashrc

Download the swissprot database::
   
   wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
   gunzip uniprot_sprot.fasta.gz

Index the database::

   diamond makedb --in  uniprot_sprot.fasta -d  uniprot_sprot

Link in the fasta::

   ln -s /salmon/nema_inList.fa
   mkdir ./tempDir

Run Diamond::

   diamond  blastx -d uniprot_sprot -q nema.fa -a matches -t tempDir -p 16

Convert to blast tab-separated output format::

   diamond view -a matches.daa -o matches.m8
