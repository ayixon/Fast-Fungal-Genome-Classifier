# Fast-Fungal-Genome-Classifier
A tool for quickly and accurately classifying fungal genomes

This soft extends the phylophenetic concept of species for its application in the classification of fungal genomes. It evaluates several working hypotheses:

Genomic coherence, monophyly, and the speciation or coalescence hypothesis under the Bayesian model of Poisson tree processes (bPTP). This allows the evaluation
of phenetic, genomic, phylogenetic and evolutionary-molecular elements in a corpus of speciation.

    Usage: ./fclassifier.sh -i <input_file> -d <database_file> -m <model>

    Options:
    -i <input_file>    Input fasta, fna, or fa archive
    -d <database_file> Database file in .msh format
    -m <model>         Select between two species delimitation methods mptp or bptp
    -h                 Display this help message
   More information on the species delimitation methods here [mPTP](https://academic.oup.com/bioinformatics/article/33/11/1630/2929345) or [bPTP](https://academic.oup.com/bioinformatics/article/29/22/2869/314968)
    
Before you begin, install the following dependencies: 

> **you might want to try mamba instead of conda in the following commands**

	sudo apt install apcalc
	sudo apt install ncbi-entrez-direct
	sudo apt install ncbi-blast+ 
	sudo apt-get install -y mptp
	conda install -c bioconda fastani
	conda install -c bioconda jolytree
	conda install -c bfurneaux bptp
	
    
 > **Download preconfigured Fungal** database here [Fungal Mash DB](https://drive.google.com/file/d/1C3nl1MJjJytGJ9_F0A-fCSJRoLwVCcJM/view)

Contact: Ayixon Sánchez-Reyes  ayixon@gmail.com 
Computational Microbiology                                                                              
Microbiological Observatory                                                                             
Institute of Biotechnology, UNAM, Cuernavaca, MEXICO                                                    
                                                                                                     
============                                                                                              
= VERSIONS =      VERSION=1.0. Written by Ayixon Sánchez Reyes                      
