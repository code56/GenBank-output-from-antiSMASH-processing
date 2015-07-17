# the Synprio.py is the main program that is calling the antismash_processing.py and MGB_processing.py

The antismash_process.py is taking the GenBank output files from the from antiSMASH software (http://antismash.secondarymetabolites.org/) 
which outputs GenBank cluster files. 

The Synprio also receives the computational outputs of the MGB_processing.py file. 

The aim of the Synprio.py is to classify and rank BioSynthetic Gene Clusters according to their features including GC content,
the length of the genes, the length distance of the largest gene from the user's preference, the total numbe of codons in the cluster
the synonymous codon usage of the cluster and the distance from the E.coli's SCU. 

For more information on the aims and how to run the Synprio.py please refer to the manpage.


