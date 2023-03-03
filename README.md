# Fast-Fungal-Genome-Classifier
A tool for quickly and accurately classifying fungal genomes

This soft extends the phyllophenetic concept of species for its application in the classification of fungal genomes. It evaluates several working hypotheses:

Genomic coherence, monophyly, and the speciation or coalescence hypothesis under the Bayesian model of Poisson tree processes (bPTP). This allows the evaluation
of phenetic, genomic, phylogenetic and evolutionary-molecular elements in a corpus of speciation.

    Usage: ./fclassifier.sh -i <input_file> -d <database_file> -m <model>

    Options:
    -i <input_file>    Input fasta, fna, or fa archive
    -d <database_file> Database file in .msh format
    -m <model>         Select between two models: mptp or bptp
    -h                 Display this help message
