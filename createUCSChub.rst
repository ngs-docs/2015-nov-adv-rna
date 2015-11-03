Make and share your reference based transcriptome
=================================================

:author: Tamer Mansour
:date: Nov 3, 2015

This tutorial shows how to build a new transcriptome with Cufflinks
and share it as a custom track on the UCSC Genome Browser, using github
to store the files.

Background documentation:

   https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html

Launch an EC2 instance
----------------------
Go ahead and open the EC2 service to lunch a suitable instance. 
Change your location on the top right corner of the page to be US East (N. Virginia).
Push the **Launch Instance** button then follow these instructions:

1. On "Step 1: Choose an Amazon Machine Image (AMI)": The protocol was prepared to run on "Ubuntu Server 14.04 LTS (HVM)"
2. On "Step 2: Choose an Instance Type": You can run this protocol on **m3.xlarge instance** (4 vCPU and 15 GiB memory).
3. On "Step 3: Configure Instance Details": accept the default settings
4. On "Step 4: Add storage": accept the default settings
5. On "Step 5: Tag Instance": Give your instance a name.
6. On "Step 6: Configure Security Group": Create a new security group and add security rules to enable ports 22, 80, and 443 (SSH, HTTP, and HTTPS).
7. On "Step 7: Review Instance Launch": Review the information of your instance. You will see an alarm about the security of your instance. This is ok, once you click launch, you will be able to make a security key-pair (or use yours if you already have one)
8. After you launch your instance, the confirmation page will show your instance ID. Click the instance ID to watch your instance status. In the lower half of the page, you will description info of your instance. Copy the public IP address to use in the next step.  

Logging into your new cloud instance (Windows version)
------------------------------------------------------
You need an SSH client to connect from you own computer to the Amazon instance you just started in the cloud:

1. Download PuTTY and PuTTYgen from: http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html
2. Run PuTTYgen, Find and load your '.pem' file, and Now, "save private key"
3. Run PuTTY, paste the public DNS in the host name. In the left category panel, expand the SSH entry, select Auth, find your private key and click open
4. click yes if prompted then Log in as "ubuntu"
5. create a folder to contain all the files of this experiment and define this path

Install software
----------------

Define your working directory::

    mkdir evalTrans
    workingPath=$"/home/ubuntu/evalTrans"
    export workingPath

Install prerequisites
::
   
   cd ~
   sudo apt-get update
   sudo apt-get install gcc g++ pkg-config wget make
   
install cufflinks/2.1.1  
::

   cd /usr/src/
   sudo wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.1.1.Linux_x86_64.tar.gz
   sudo tar -zxvf cufflinks-2.1.1.Linux_x86_64.tar.gz
   cd /usr/bin
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/cufflinks .
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/cuffmerge .
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/cuffcompare .
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/cuffdiff .
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/gffread .
   sudo ln -s /usr/src/cufflinks-2.1.1.Linux_x86_64/gtf_to_sam .
   
Get the data
------------
::

   cd $workingPath
   mkdir data && cd data
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/control1.tar.gz
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/control2.tar.gz
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/control3.tar.gz
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/test1.tar.gz
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/test2.tar.gz
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/cufflinkData/test3.tar.gz
   for f in *.tar.gz; do tar zxvf $f;done


Run Cufflinks
-------------
::

   for f in $(pwd)/*.bam; do
     label=${f%.bam}
     mkdir ${label}_cufflinks && cd ${label}_cufflinks
     cufflinks --label ${label} --num-threads 4 $f
   done

Run Cuffmerge
-------------
::
  
   cd ../
   for dir in $(pwd)/*_cufflinks; do
     echo $dir/transcripts.gtf; done > assemblies.txt

   cuffmerge -o isofrac0.05 --num-threads 4 --min-isoform-fraction 0.05 assemblies.txt 2&> cuffisofrac0.05.log
   cuffmerge -o isofrac0.2 --num-threads 4 --min-isoform-fraction 0.2 assemblies.txt 2&> cuffisofrac0.2.log
   cuffmerge -o isofrac0.5 --num-threads 4 --min-isoform-fraction 0.5 assemblies.txt 2&> cuffisofrac0.5.log
   cuffmerge -o isofrac0.7 --num-threads 4 --min-isoform-fraction 0.7 assemblies.txt 2&> cuffisofrac0.7.log
   cuffmerge -o isofrac0.9 --num-threads 4 --min-isoform-fraction 0.9 assemblies.txt 2&> cuffisofrac0.9.log


Download UCSC tools & couple custom scripts
-------------------------------------------
:: 

   mkdir $workingPath/UCSC_kent_commands
   cd $workingPath/UCSC_kent_commands
   wget -r --no-directories ftp://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
   wget -r --no-directories ftp://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed
   wget -r --no-directories ftp://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
   wget -r --no-directories ftp://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes
   chmod 755 *

   cd $workingPath
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/create_trackHub.sh
   curl -L -O https://github.com/ngs-docs/2015-nov-adv-rna/raw/master/edit_trackDb.sh

Initiate the basic structure for horse track hubs
------------------------------------------------
::

   UCSCgenome=$"equCab2"
   hub_name=$"testhub1"
   shortlabel=$"CompIsoformFrac"
   longlabel=$"UCSC track hub to compare selection isoform fraction"
   email=$"youremail@somthing.com"
   mkdir -p $workingPath/track_hub/$UCSCgenome/BigBed   
   cd track_hub
   bash $workingPath/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"


Convert GTF files to BigBed files
---------------------------------
::

   cd $workingPath/data
   $workingPath/UCSC_kent_commands/fetchChromSizes $UCSCgenome > chromSizes.txt
   > $workingPath/data/UCSC_assemblies.txt
   for assembly in $(pwd)/isofrac*; do
     echo $assembly
     cd $assembly
     $workingPath/UCSC_kent_commands/gtfToGenePred merged.gtf merged.gpred
     $workingPath/UCSC_kent_commands/genePredToBed merged.gpred merged.bed
     sort -k1,1 -k2,2n merged.bed > merged_sorted.bed
     $workingPath/UCSC_kent_commands/bedToBigBed merged_sorted.bed $workingPath/data/chromSizes.txt merged.BigBed
     identifier=$(basename $assembly)
     cp merged.BigBed $workingPath/track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
     echo $identifier >> $workingPath/data/UCSC_assemblies.txt
   done


Populate the track DB  
---------------------
::
  
   trackDb=$workingPath/track_hub/$UCSCgenome/trackDb_$shortlabel.txt 
   > $trackDb
   bash $workingPath/edit_trackDb.sh "$trackDb" "$workingPath/data/UCSC_assemblies.txt"


Upload your data to github
--------------------------
Reference: https://help.github.com/articles/adding-an-existing-project-to-github-using-the-command-line/
:
 
1. Make a github account at https://github.com/
::
    
2. Create a new repository (https://help.github.com/articles/creating-a-new-repository/)
::
    
3. Install Git to your Amazon instance
::
    sudo apt-get install git

4. Initialize the local directory as a Git repository
::
    cd $workingPath
    git init 
 
5. Add the files in your new local repository
::
    git add track_hub
 
6. Commit the files that you've staged in your local repository.
::
    git commit -m "upload track hub"
 
7. At the top of your GitHub repository's Quick Setup page, copy the remote repository URL(the HTTPS one).
::
    

8. In Terminal, add the URL for the remote repository 
::
    git remote add origin <remote repository URL>  
 
9. Push the hub directory in your local repository to GitHub.
::
    git push origin master
 
Visualize your tracks in UCSC
-----------------------------
get the URL of the raw hub_CompIsoformFrac.txt and add to your tracks on UCSC



