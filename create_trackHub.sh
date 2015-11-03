UCSCgenome="$1"
hub_name="$2"
shortlabel="$3"
longlabel="$4"
email="$5"

## create the hub file
echo "hub $hub_name" >> hub_$shortlabel.txt
echo "shortLabel $shortlabel" >> hub_$shortlabel.txt
echo "longLabel $longlabel" >> hub_$shortlabel.txt
echo "genomesFile genomes_$shortlabel.txt" >> hub_$shortlabel.txt
echo "email $email" >> hub_$shortlabel.txt


## create the genomes file
echo "genome $UCSCgenome" >> genomes_$shortlabel.txt
echo "trackDb $UCSCgenome/trackDb_$shortlabel.txt" >> genomes_$shortlabel.txt




