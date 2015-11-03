Get KEGG Annotations
===================================

We'll start from the m3.xlarge Amazon machine booted & configured in
`salmon.rst <salmon.rst>`__.  If you are just running this, you'll need
to run the apt-get commands, install khmer, and mount the data snapshot
before continuing. You should also have completed the blastx run from
`diamond_blast.rst <diamond_blast.rst>`__.

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

Link the blastx results file into your current working directory::

   ln -s $HOME/diamond/matches.m8 ./


Download and run a python script to relate your blastx results to their KEGG Orthology numbers
using the linkage database information::

   curl -O https://raw.githubusercontent.com/ngs-docs/2015-nov-adv-rna/master/getKOfromBlastx.py

   python getKOfromBlastx.py --blastx matches.m8 --spToKO swiss_orthology.list 


See the results of the file with::

   less matches_keggOrthology.txt

Type 'q' to exit out of the less file viewer.


Visualize Kegg Pathways Present in the Data 
-------------------------------------------

Download and run a python script to visualize the KEGG annotations found in the metabolic pathways map::

   #code coming once it's finished
   #KEGG REST api is down - can't test the visualization 


See https://github.com/widdowquinn/notebooks/blob/master/Biopython_KGML_intro.ipynb for a full
introduction to the Biopython KGML rendering module.








