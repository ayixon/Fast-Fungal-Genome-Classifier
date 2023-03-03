#!/bin/bash 
	
#"FAST FUNGAL GENOME CLASSIFIER":  A TOOL FOR QUICKLY AND ACCURATELY CLASSIFYING FUNGAL GENOMES
####################################################################################################################

####################################################################################################################
#                                                                                                                  #
#COPYRIGHT="Copyright (C) 2022 Ayixon Sánchez Reyes"                                                               #
#                                                                                                                  #
# This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU         #
# General Public License as published by the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.This program is distributed WITHOUT ANY WARRANTY.                                              #
#                                                                                                                  #
####################################################################################################################
####################################################################################################################
# DEPENDENCIES: JolyTree, apcalc, ncbi-entrez-direct, orthoani, Blast, Biopython, bPTP,mPTP                        #  
#                                                                                                                  #
##########Before you begin, install the following:                                                                 #
#
##########	sudo apt install apcalc
##########	sudo apt install ncbi-entrez-direct
##########  sudo apt install ncbi-blast+ 
########## 	conda install -c bioconda fastani
########## 	conda install -c bioconda jolytree
##########  https://github.com/Pas-Kapli/mptp
##########  conda install -c bfurneaux bptp     
  
####################################################################################################################
#                                                                                                                  #
####################################################################################################################
# Usage: ./fclassifier.sh -i query_genome -d database.msh - m <model> mptp or bptp (bPTP by deafult)	                                                                                   #
# The working directory must contain the mash database (.msh) and the query genome in fasta format                 #
####################################################################################################################

# Rational: Compare a query_genome vs a curated MASH database, select the nearest phylogenetic neighbors;          #
# Estimate the ANI of the query vs the references,                                                                 #
# JolyTree estimate a rapid and robust genome based phylogeny 
# The tree is subjected to speciation hypothesis testing under Poisson Tree Processes Model                                                    #
                                                                                                            #
# THIS GENOME CLASSIFIER deals with the "Phylophenetic Species Concept" 
# AND WITH "Molecular Species Delimitation" by testing three working hypotheses:                            #
     # The Genomic Coherence measured through the genomic distance of Mash and the ANI                      #
     # The Phylogenetic Hypothesis of monophyly
     # The speciation or coalescence hypothesis under the Poisson tree processes model                      #
     # This makes it possible to evaluate phenetic, genomic and evolutionary-molecular 
     # elements in a corpus of speciation                                                                   #
                                                                                                            #
#  Ayixon Sánchez-Reyes                                                                   ayixon@gmail.com  #
#  Computational Microbiology                                                                               #
#  Microbiological Observatory                                                                              #
#  Institute of Biotechnology, UNAM, Cuernavaca, MEXICO                                                     #
#############################################################################################################
#                                                                                                           #
# ============                                                                                              #
# = VERSIONS =      VERSION=1.0. Written by Ayixon Sánchez Reyes                                            # 
#   VERSION=1.0                                                                                             #                       
# ============                                                                                              #
#############################################################################################################
#############################################################################################################

# Use getopts to parse the input options

start=`date +%s.%N`

function display_usage {
  echo "Usage: $0 -i <input_file> -d <database_file> -m <model>"
  echo
  echo "Options:"
  echo "  -i <input_file>    Input fasta, fna, or fa archive"
  echo "  -d <database_file> Database file in .msh format"
  echo "  -m <model>         Select between two models: mptp or bptp"
  echo "  -h                 Display this help message"
}

input_file=""
database_file=""
model=""

while getopts "i:d:m:h" opt; do
  case $opt in
    i)
      input_file="$OPTARG"
      ;;
    d)
      database_file="$OPTARG"
      ;;
    m)
      model="$OPTARG"
      ;;
    h)
      display_usage
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z "$model" ]; then
  model="bptp"
fi

if [ "$model" == "mptp" ]; then
  # Execute mptp command
  echo "Executing mptp command"
else
  echo "bptp selected by default"
fi

# Estimate Mash distance with the input files
  mash dist $input_file $database_file -p 10| sort -n -k3 > output.mash.txt; 
 echo "" 
 echo -e "	\e[0;32mMASH analysis finished. These are the top hits   \e[0m " 

 echo "	Query	Reference	D	p_value	shared hashed"

head output.mash.txt


cut -f3 output.mash.txt |head -n 100 > distancias 


uniq distancias > dist.uniq 

 echo "" 

 echo -e " \e[0;32mPrinting the distance table (Mash D)   \e[0m " 

 echo "" 

echo -e "	\e[0;32mThese are the closest genomes to your query and their approximate ANI   \e[0m " 

# Approximate ANI estimation as 1-Distance
echo "" 

for i in $(fmt dist.uniq)
do  sort -n -k3 output.mash.txt| grep -m1 -F "$i" 
echo "" 
var1=$(calc 1-"$i")
echo -e " \e[0;32mComputing approximate ANI   \e[0m "
echo $var1
echo "---------"  
if ((`bc <<< "$var1>=0.95"`)) ; then
  echo -e '\e[0;33mGenomic Coherence Detected   \e[0m '
  echo "" 
  echo "--------------------------" 
   
fi
  done 
echo "--------" 
echo "" 

#Selection and download of the neighboring genomes with the smallest genomic distance

echo -e " \e[0;32mDownloading the reference genomes from NCBI   \e[0m " 
echo "#####################################" 

for z in $(fmt dist.uniq)
do grep -F -m 5 "$z" output.mash.txt| awk '{print $2}'| sed 's/_[A-Z]/   / gI' | awk '{print $1}' >> Genome_accnumber.txt
  done
   echo ""


cat Genome_accnumber.txt | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget "$url/$fname" ;
        done ;
		
    done
echo "" 
echo -e " \e[0;32mDownload successful   \e[0m "
echo "" 
gzip -d *.gz
#Genomes are copied to a separate directory for phylogenetic analysis
#working directory cleanup operations

echo -e " \e[0;32mStarting Phylogenetic Analysis with JolyTree   \e[0m "
echo "##############################################" 

mkdir Mash_out
cp dist.uniq distancias Mash_out/ 

mkdir JolyTree_in 
cp *.fna JolyTree_in/
cp *.fasta JolyTree_in/
cp *.fa JolyTree_in/

JolyTree.sh -i JolyTree_in/ -b out_tree -t 10 
	
mkdir JolyTree_out
cp out_tree* JolyTree_out/
echo ""
echo -e " \e[0;32m#Phylogenetic analysis finished, results are in the directory JolyTree_out  \e[0m" 
echo ""

echo -e " \e[0;32m#Computing fastANI  \e[0m "
echo "--------------------"
echo  Genome_Query: $input_file 

rm *.fna
rm out_tree*

cp Genome_accnumber.txt output.mash.txt Mash_out/ 
rm dist.uniq distancias output.mash.txt Genome_accnumber.txt

cd JolyTree_in/
for ani in *.fna
do fastANI -q $input_file -r $ani -t 4  --fragLen 1000 --minFraction 0.1 -o $ani.txt
done
echo "--------------------"
echo -e " \e[0;32m#fastANI done \e[0m "
echo "--------------------"
awk '{print $1,$2,$3}' *.txt
#awk '{print substr(FILENAME, 3)" \""$0"\" "}' ./*.txt|sort -t\" -nrk2 | sed 's/_[A-Z]/   / gI'| awk '{print $1,"ANI=", $NF}'

cd ..


mkdir fastANI_out
cp JolyTree_in/*.txt fastANI_out/
echo "" 
rm JolyTree_in/*.txt
rm *.gz
echo ""

echo -e " \e[0;32m#Computing Species Delimitation under Markov Chain Monte Carlo  \e[0m "

mkdir sp_Results; \

cd JolyTree_out/; \
#------------------------------------------------------------------------------------------------
if [[ $model == mptp ]]; then
 mptp --mcmc 500000 --single --mcmc_sample 10000 --mcmc_burnin 100000 \
--tree_file out_tree.nwk --output_file PTP_sp --tree_show ; \

else
  echo -e " \e[0;32m#mptp was not selected, it is ok; lets try bptp algorithm  \e[0m "
  
fi

if [[ $model == bptp ]]; then
 bPTP.py -t out_tree.nwk  -o PTP_sp  -s 1234 -r -i 1000000
fi
#------------------------------------------------------------------------------------------------
cd .. ; \

cp JolyTree_out/PTP*  sp_Results/ ; \

echo -e " \e[0;32mSpeciation test done, look at the sp_Results directory  \e[0m "

echo ""
	echo -e " \e[0;32mFast_Genome_Classifier pipeline done  \e[0m "
echo ""
echo -e " \e[0;33mWe recommend that you test the phylophenetic hypothesis about *Diagnostic Characters*  \e[0m "

echo -e " \e[0;33mYou may want to explore the following phenotypic determinants:  \e[0m "

echo -e " \e[0;33mMenaquinone biosynthesis including the menaquinone reductase gene (MenJ)  \e[0m "

echo -e " \e[0;33mFatty acid synthesis cycle FASI, FASII \e[0m "

echo -e " \e[0;33m...any other phenotypic determinant of interest  \e[0m "

echo ""
	echo "Thank you for your support, enjoy it"
	echo ""

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc )
echo Time $runtime seconds
