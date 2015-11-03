# 7.1.2015

# Builtins
import os, sys
# from armchair biology blog (now implemented in biopython)
from KGML_parser import read
from KGML_vis import KGMLCanvas
from KGML_scrape import retrieve_kgml_to_file, retrieve_KEGG_pathway
import argparse
############################
psr = argparse.ArgumentParser(description="Color KEGG maps with the Kegg Orthology entries that exist in your annotated transcriptome.")
# kegg pathway name
psr.add_argument('--path',help='Kegg Pathway Name', dest="path")
# folder for map output
psr.add_argument('--output',help='name of output folder', default = './', dest="outDir")
# KO (for just presence/absence)
psr.add_argument('--blastKO', help='Kegg Orthology from blast', dest='KO')
##############################
args = psr.parse_args()
##############################

if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)

pathway = retrieve_KEGG_pathway(args.path) #pathway of interest
#pathway = read(args.path)
#import pdb;pdb.set_trace()

#get KO numbers from previously created file
ko = ['ko:'+ ko for (trans,sp,ko) in (line.strip().split("\t") for line in (open(args.KO, 'r')))]

#knownKOSet = set([e for e in pathway.entries.values() if len(set(e.name.split()).intersection(ko))])
knownKOSet = set([e for e in pathway.orthologs if len(set(e.name.split()).intersection(ko))])

def colorMapItems(geneSet, color, width):
    for e in geneSet:
	for g in e.graphics:
	    if g.type == 'line':
		g.fgcolor = color
		g.width = width
	    g.bgcolor = color

#main
enhanceSet = knownKOSet 
#notDE = set([e for e in pathway.entries.values() if not len(set(e.name.split()).intersection(enhanceSet)) and e.type != 'map'])
notDE = set([e for e in pathway.orthologs if not len(set(e.name.split()).intersection(enhanceSet))])

kgml_map = KGMLCanvas(pathway, show_maps=True)
kgml_map.import_imagemap = True  
#kgml_map.show_maps = False
kgml_map.show_orthologs = True
kgml_map.draw_relations = False
kgml_map.show_compounds = False
kgml_map.show_genes = False

colorMapItems(notDE,'#D3D3D3', 1)
os.chdir(args.outDir)
koInMap = open(args.path + '_KO.txt', 'w')

colorMapItems(knownKOSet,'#666666', 10)
for k in knownKOSet:
    koInMap.write(k.name + '\t' + 'present' + '\n')
koInMap.close()

# And rendering elements as an overlay
#kgml_map.show_compounds = True
#kgml_map.show_genes = True
#kgml_map.show_orthologs = True
#kgml_map.draw_relations = True

kgml_map.draw(args.path + '.pdf')



