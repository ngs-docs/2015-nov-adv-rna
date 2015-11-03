Visualize KEGG Pathways
===================================

We'll start from the m3.xlarge Amazon machine booted & configured in
`salmon.rst <salmon.rst>`__.  If you are just running this, you'll need
to run the apt-get commands, install khmer, and mount the data snapshot
before continuing. 

----

First, create a working directory and download the swiss-prot to kegg orthology
database linkage file:: 

   cd 
   mkdir kegg
   cd kegg
   curl -O https://raw.githubusercontent.com/ngs-docs/2015-nov-adv-rna/master/swiss_orthology.list


Note: this file was downloaded on 11/02/2015, but this database is updated regularly. In order to
download a recent version, go to http://www.genome.jp/linkdb/ and click on "SWISS-PROT", then "Orthology," 
then click the "Download" button just under the database linkage map.

Getting Kegg Orthology Numbers
------------------------------

We need the blastx results from `diamond_blast.rst <diamond_blast.rst>`__, which you can download with:: 

   curl -O https://s3-us-west-1.amazonaws.com/dib-training.ucdavis.edu/mRNAseq-non-2015-05-04/nema.x.swissprot.diamond.m8.gz
   gunzip nema.x.swissprot.diamond.m8.gz 

Download and run a python script to relate your blastx results to their KEGG Orthology numbers
using the linkage database information::

   curl -O https://raw.githubusercontent.com/ngs-docs/2015-nov-adv-rna/master/getKOfromBlastx.py

   python getKOfromBlastx.py --blastx nema.x.swissprot.diamond.m8 --spToKO swiss_orthology.list 


See the results of the file with::

   less nema.x.swissprot.diamond_keggOrthology.txt

Type 'q' to exit out of the file viewer.


Visualize Kegg Pathways Present in the Data 
-------------------------------------------

Install some extra python modules to visualize KEGG maps::

   sudo pip install reportlab
   git clone https://github.com/bluegenes/KGML.git


Run the python script to visualize the KEGG annotations. ko01100 is the KEGG Metabolic Pathways overview map:: 
 
   python ./KGML/simpleDrawKeggMap.py --path ko01100 --blastKO nema.x.swissprot.diamond_keggOrthology.txt 


Copy the pdf file to your computer so you can open it with a pdf viewer. The pathways that were found in the
nema data are displayed as thickened dark gray lines::

   # This code uses the same information we used to login:
   # KEY = your private key
   # YOUR-AMAZON-MACHINE = the "Public DNS" of your amazon machine
   scp -i KEY.pem ubuntu@YOUR-AMAZON-MACHINE:~/kegg/ko01100.pdf ~/Desktop

This script can be used to visualize any KEGG pathway. For example, TCA cycle::

   #run script
   python ./KGML/simpleDrawKeggMap.py --path ko00020 --blastKO nema.x.swissprot.diamond_keggOrthology.txt
   #download pdf file(s)
   scp -i your-key.pem ubuntu@your-amazon-machine.amazonaws.com:~/kegg/*.pdf ~/Desktop


In this map, orthologs are represented as boxes, rather than lines (as in the metabolic map, above).
Again, orthologs that are colored dark gray are present in your data.











