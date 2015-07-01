__author__ = 'Evanthia_Kaimaklioti'
#this is the script manipulating the antiSMASH output GenBank cluster files
#this pipeline was constructed for the purposes of the RP1 project within the MSc in Bioinformatics and Systems Biology
#course of The University of Manchester, UK
# this is a Python script

def count_codons(cds):                #cds is the iteration variable; cds is the sequence from which the codons will be found and then counted
    usage = {}                        # it's how you initiate a dictionary - Check though. #it only applies within the loop.
                                     # local variable
    for i in range(0, len(cds), 3):
        codon = cds[i:i+3]
        if usage.has_key(codon):            #if the codon is already in the dictionary, add one to it.
            usage[codon] += 1               # Is this the right place that I should compute the [(this_codon_total_number/total_number_of_codons?? )]
        else:                              # and the other tricky part is that I need to calculate the total number of codons of this cluster...
            usage[codon] = 1
    return usage


def get_keys(dict1):
    for key, value in dict1.iteritems():
        return key

def get_values(dict):
    for key, value in dict.iteritems():
        return value

def geneProcess(inputFilename,outputFilename, outputFilename2, ideal_GC_usage, ideal_len):
    from Bio import SeqIO, SeqFeature, SeqRecord
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    from Bio.SeqUtils import GC
    import sys
    import argparse
    import numpy
    import operator
    import re
    from operator import itemgetter
    from collections import OrderedDict

    size_list = []
    record = SeqIO.read(inputFilename, "gb")
    output_handle = open(outputFilename, "w")

    count = 0
    total_seq = 0    #keep the count of the total sequence length (keep it outside the loop)
    total_gc = 0
    largest_seq = 0
    largest_score = 0
    gc_dict = {}
    cds_list = []
    joined_list = []


    for feature in record.features:
        if feature.type == "CDS":
            count = count + 1

            feature_name = feature.qualifiers                   #instead of feature.qualifiers, I can write the product name
                                                                # or the actual coordinates or complement()
                                                                ## in the origin    #sequence #or also it could be the locus_tag="
            feature_seq = feature.extract(record.seq)           # the key here is the feature.extract, it extracts the CDS sequence/

            sequence_cds = feature_seq                          #created this because I got the AtrributeError: 'Seq' object has no attribute 'append'
            cds_list.extend(feature_seq)                        #list to store all the cds sequences in this cluster
            joined_list=''.join(cds_list)                       #gives the cds sequences of the cluster as one big cds in a list, a big long string of nucleotides.
            total_codons = len(joined_list)                     # total number of codons in this cluster

            length_seq = len(feature_seq)                       #computes the length of this sequence
            gc_seq = GC(feature_seq)                            #finds the GC_content (in percentage %) of the sequence. this should automatically cope with mixed case sequences and the ambiguous
                                                                #and mixed case nucleotide S which means G or C.

            sequence = feature_seq                              #check later if this is necessary to create another variable.
            for sequence in range(0, 20000):                    #large enough value
                if length_seq > largest_seq:                    #an if test to find the largest in length CDS in the cluster
                    largest_seq = length_seq                    #lengthiest sequence in the cluster
            biggest_seq = feature_seq                           # this is to print the sequence of this largest in length sequence, but I can remove this "biggest_seq" variable might not be necessary to see it.


            prot_seq = feature_seq.translate(table=11, to_stop=True, cds=True)              # the sequence is translated into protein sequence using the table 11. notice: when using the to_stop argument, the stop
                                                                #codon itself is not translated-and the stop symbol is not included at the end of your protein sequence.
            

            seqGC = ((gc_seq*length_seq)/100)                   #seqGC: finds the NUMBER of GCs of each sequence is the cluster file
                                                                #do I need to make it as int , so it isn't 130.0?

            total_gc = seqGC + total_gc                         #the total GC of the cluster file, adds GCs as it reads the CDSs in the cluster (GenBank file)


            total_seq = length_seq + total_seq                  #the total length of the coding sequences in the cluster file


            cluster_gcoverall = ((total_gc/total_seq)*100)      #the overall GC content of the cluster

           # Simple FASTA output without line wrapping:
            output_handle.write(">{0}\n{1}\n".format(str(feature_name), str(feature_seq) + "\n" + "translation:" +
                                                     str(prot_seq)))
            output_handle.write("The length of this sequence is" + " " + str(length_seq) + " " + "nucleotides" +
                                ", and its GC content is" + " " + str(gc_seq) + "%" + " " + "the number of GCs is"
                                 + " " + str(seqGC) + " " + " " + "total length of sequences is" + " " + str(total_seq) +
                                " " + "The overall GC number of the cluster SO FAR  is" + " " + str(total_gc) +
                                "\n" + "\n")


    output_handle.close()                                       #close/stop reading this GenBank file


    gc_values = str(gc_seq)
    sorted_gc = sorted(gc_values)

    FreqCodon = count_codons(joined_list)                  #list variable that stores the dictionary in a

    SynonymousCodons = {
        'CYS': ['TGT', 'TGC'],
        'ASP': ['GAT', 'GAC'],
        'SER':['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
        'GLN': ['CAA', 'CAG'],
        'MET': ['ATG'],
        'ASN':['AAC', 'AAT'],
        'PRO':['CCT', 'CCG', 'CCA', 'CCC'],
        'LYS': ['AAG', 'AAA'],
        'STOP':['TAG', 'TGA', 'TAA'],
        'THR':['ACC', 'ACA', 'ACG', 'ACT'],
        'PHE':['TTT','TTC'],
        'ALA': ['GCA', 'GCC', 'GCG', 'GCT'],
        'GLY':['GGT', 'GGG', 'GGA', 'GGC'],
        'ILE':['ATC', 'ATA', 'ATT'],
        'LEU':['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
        'HIS':['CAT', 'CAC'],
        'ARG':['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
        'TRP':['TGG'],
        'VAL':['GTA', 'GTC','GTG', 'GTT'],
        'GLU':['GAG', 'GAA'],
        'TYR':['TAT', 'TAC']

    }

    FracCodon1 = {
            'CYS': {'TGT':0.00, 'TGC':0.00},
            'ASP': {'GAT':0.00, 'GAC':0.00},
            'SER': {'TCT':0.00, 'TCG':0.00, 'TCA': 0.00, 'TCC':0.00, 'AGC':0.00, 'AGT':0.00},
            'GLN': {'CAA':0.00, 'CAG':0.00},
            'MET': {'ATG':0.00},
            'ASN': {'AAC':0.00, 'AAT':0.00},
            'PRO': {'CCT':0.00, 'CCG':0.00, 'CCA': 0.00, 'CCC':0.00},
            'LYS': {'AAG':0.00, 'AAA':0.00},
            'STOP':{'TAG':0.00, 'TGA':0.00, 'TAA': 0.00},
            'THR': {'ACC':0.00, 'ACA':0.00, 'ACG': 0.00, 'ACT':0.00},
            'PHE': {'TTT':0.00,'TTC':0.00},
            'ALA': {'GCA':0.00, 'GCC':0.00, 'GCG': 0.00, 'GCT':0.00},
            'GLY': {'GGT':0.00, 'GGG':0.00, 'GGA': 0.00, 'GGC':0.00},
            'ILE': {'ATC':0.00, 'ATA':0.00, 'ATT': 0.00},
            'LEU': {'TTA':0.00, 'TTG':0.00, 'CTC': 0.00, 'CTT':0.00, 'CTG':0.00, 'CTA':0.00},
            'HIS': {'CAT':0.00, 'CAC':0.00},
            'ARG': {'CGA':0.00, 'CGC':0.00, 'CGG': 0.00, 'CGT':0.00, 'AGG':0.00, 'AGA':0.00},
            'TRP': {'TGG':0.00},
            'VAL': {'GTA':0.00, 'GTC':0.00,'GTG': 0.00, 'GTT':0.00},
            'GLU': {'GAG':0.00, 'GAA':0.00},
            'TYR': {'TAT':0.00, 'TAC':0.00}
    }

    for keyS in SynonymousCodons:                                                           #for key (e.g. 'CYS') in SynonymousCodons
        totalCodon = 0                                                                      #each time you go to each key (e.g. 'CYS' you tear (make zero) the additions
        for codonS in SynonymousCodons[keyS]:                                               # 'x' is the value(e.g.'TGT') of this key 'CYS'
            for keyF in FreqCodon:
                if keyF == codonS:
                    FracCodon1[keyS][codonS] = FreqCodon[keyF]
                    totalCodon = totalCodon + FreqCodon[keyF]                               #computes the total number of codons (i.e. the synonymous codons) that code for this key in the sequence
        for codonS in SynonymousCodons[keyS]:                                               
            FracCodon1[keyS][codonS] = float(FracCodon1[keyS][codonS])/float(totalCodon) 
#

    e_coli_usage = {
        'CYS': {'TGT':0.45, 'TGC':0.55},
        'ASP': {'GAT':0.63, 'GAC':0.37},
        'SER': {'TCT':0.15, 'TCG':0.15, 'TCA': 0.12, 'TCC':0.15, 'AGC':0.28, 'AGT':0.15},
        'GLN': {'CAA':0.35, 'CAG':0.65},
        'MET': {'ATG':1.00},
        'ASN': {'AAC':0.55, 'AAT':0.45},
        'PRO': {'CCT':0.16, 'CCG':0.52, 'CCA': 0.19, 'CCC':0.12},
        'LYS': {'AAG':0.23, 'AAA':0.77},
        'STOP':{'TAG':0.07, 'TGA':0.29, 'TAA': 0.64},
        'THR': {'ACC':0.44, 'ACA':0.13, 'ACG': 0.27, 'ACT':0.17},
        'PHE': {'TTT':0.57,'TTC':0.43},
        'ALA': {'GCA':0.21, 'GCC':0.27, 'GCG': 0.36, 'GCT':0.16},
        'GLY': {'GGT':0.34, 'GGG':0.15, 'GGA': 0.11, 'GGC':0.40},
        'ILE': {'ATC':0.42, 'ATA':0.07, 'ATT': 0.51},
        'LEU': {'TTA':0.13, 'TTG':0.13, 'CTC': 0.10, 'CTT':0.10, 'CTG':0.50, 'CTA':0.04},
        'HIS': {'CAT':0.57, 'CAC':0.43},
        'ARG': {'CGA':0.06, 'CGC':0.40, 'CGG': 0.10, 'CGT':0.38, 'AGG':0.02, 'AGA':0.04},
        'TRP': {'TGG':1.00},
        'VAL': {'GTA':0.15, 'GTC':0.22,'GTG': 0.37, 'GTT':0.26},
        'GLU': {'GAG':0.31, 'GAA':0.69},
        'TYR': {'TAT':0.57, 'TAC':0.43}
    }

    difference_dict = {}

    for key in e_coli_usage:                                                                     #for key (e.g. 'CYS') in e_coli_usage
        for x in e_coli_usage[key]:                                                             # 'x' is the value (e.g. 'TGT' and then 'TGC' of this key 'CYS'
            for y in FracCodon1[key]:                                                             #'y' is the key (e.g. 'CYS' in the e_coli_usage
                if (x==y):                                                                      #if the 'x' (e.g. 'TGT') is found as a key in the FreqCodon then move on: #so put an else x=!y make the same equal to 0
                    difference_dict[x] = abs(e_coli_usage[key][x] - FracCodon1[key][y])          #creating the double dictionary


    diff_score = sum(difference_dict.values())                                                      ## now here is the code to sum the absolute values of the differences.


    text_file = open(outputFilename, "r")
    out_file = open(outputFilename2, "w")      #23.03.15 it seems that it is not writing new files with this

    for line in text_file:
        if line[0] == '>':
            pass                    #ignore lines which begin with ">"
        elif line[0] == '\n':       #ignore blank lines
            pass
        elif line[0] == 't':        #ignore lines that start with (translation)
            pass
        elif line[0] == 'T':        #ignore lines that start with "The"
            pass
        elif line[-2] == '%':       #ignore the lines whose 2nd from the last character is %
            pass  #also can put something to ignore anything that doesn't start with the above things.in case the previous file created something unnecessary.
        else:
            out_file.write(str(count_codons(line.rstrip('\n'))) + "\n" + "\n")    #rstip. to remove the newline character from the end of the gene sequence

    text_file.close()

    gc_diff = ((float(ideal_GC_usage)-cluster_gcoverall)**2)   #the difference of GC content to the power of 2


    gc_absdiff = math.sqrt(gc_diff)              #the absolute difference of GC content
    gc_absdiff1 = math.floor(gc_absdiff)         # converting the absolute difference of GC content from a float number to an integer
    cluster_name = inputFilename.split(".")      # splitting the cluster name (as it is given as the path to the file and splitting it every time there is a "."
    cluster_name1 = str(cluster_name[-2])        #take only the "clusterxxx" of each cluster

    gc_dict = {cluster_name1: gc_absdiff}

    distance_from_ideal_len = abs(ideal_len - largest_seq)
    size_distance = {}
    size_distance[cluster_name1] = distance_from_ideal_len

    return (get_keys(gc_dict),get_values(gc_dict),largest_seq, diff_score, get_values(size_distance))

#Further work
#Also rank and sort according to how many genes have sequence length closest to the desired sequence length. so that e.g. not to limit ourselves for one bottleneck sequence.
#also rank according to the base pairs
#also rank according to the maximal number of CDSs in the cluster. again do it by distance.
