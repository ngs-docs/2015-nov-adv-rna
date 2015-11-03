# Python script to get KEGG Orthology for transcriptome
import sys, re, os

from optparse import OptionParser

desc = """ get Kegg Orthology from Blast file"""
parser = OptionParser(description = desc)

### input ###
parser.add_option("--blastx", help = "name of the blastx output table (matches.m8)" , action="store", type="string", dest="blastx")
parser.add_option("--spToKO", help = "name of the swiss_prot to kegg orthology conversion file (from http://www.genome.jp/linkdb/)" , action="store", type="string", dest="spKO")

#note: spToKO file must be downloaded via http://www.genome.jp/linkdb/ web interface. On the map, first click "swiss-prot," then click "Orthology."
#Then click "Download" (just below the map) to download the file containing swiss-prot to kegg orthology conversion data.

(opts, args) = parser.parse_args()
#get sp2ko conversion
sp2ko = {trans.split(":")[1]: sp.split(":")[1] for (trans,sp,extra) in (line.strip().split("\t") for line in (open(opts.spKO, 'r')))}


header= ['Transcript_id', 'BLASTX','Kegg_Orthology']
out = open(opts.blastx.split('.m8')[0] + '_keggOrthology.txt', 'w')
out.write("\t".join(header) + "\n")


with open(opts.blastx) as f:
    prev = ''
    for line in f:
	line = line.strip().split('\t')
	contig = line[0]
	if contig == prev:
	    continue  #use the top e-val result, if multiple (=first hit)
	else:
	    prev = contig
	    spHit = line[1]
	    spShort = spHit.split('|')[1]
        hitValue = sp2ko.get(spShort, None)
	    #if spShort in sp2ko.keys():
        if hitValue is not None:
            keggOrthology = sp2ko[spShort]
            out.write(contig + '\t' + spHit + '\t' + keggOrthology + '\n')
    f.close()


