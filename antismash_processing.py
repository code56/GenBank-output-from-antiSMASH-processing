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


__author__ = 'Evanthia_Kaimaklioti'
#this is the script manipulating the MultiGeneBlast output files
#this pipeline was constructed for the purposes of the RP1 project within the MSc in Bioinformatics and Systems Biology
#course of the university of Manchester


from Bio import Entrez, SeqIO
import os
import re
import linecache
import operator
from operator import itemgetter
from collections import OrderedDict


def get_keys(dict1):
    for key, value in dict1.iteritems():
        return key

def get_values(dict):
    for key, value in dict.iteritems():
        return value

def get_taxonomy(NCBI_id):
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=NCBI_id, retmode="text")     #handle to fetch the NCBI record with that ID, retmode: return mode
    record = SeqIO.read(handle, "genbank")                                   #handle to read the NCBI GenBank file
    for key in record.annotations.keys():                                    # go the annotations section of the GenBank file
        key = key.lower()                                                    #make all keys lowercase
        if key == 'taxonomy':
            string_taxonomy = record.annotations['taxonomy']
            return string_taxonomy

def multiGeneBlast_process(mgb_inputfile , ideal_mgb_score):
    Entrez.email = "evanthia.kaimaklioti@postgrad.manchester.ac.uk"

    count = 0
    string_taxonomy = []                    #dictionary that returns 1 when there is a match of the query taxonomy list to the one in hand and 0 when there is no match.
    line_wanted = ""

    accession_no_query = ""
    prevline1= ""
    hit_dict = {}
    query_taxonomy = []         # the list that hold the taxonomy information for the query

    with open(mgb_inputfile) as aFile:
        for line in aFile:
            if 'Significant' in line:                              # jump the file go to "Details" and then start reading from there again
                line_wanted = str(aFile.next())
                                                                 # but can I access this line_wanted from the loop below? The word significant is before the word Details, and the Details needs to be accessed first
                m = re.match(r'\d+\.\s*([a-zA-Z0-9_]+).*', line_wanted)         #match the line that holds the accession number of the hit
                accession_no_query = m.group(1)
                print "query accession is : ", accession_no_query
                sep = "_"                                     #this applies only if the separator is there. If it isn't it will return the whole of the oldacc which is O.K. also:
                id = accession_no_query.split(sep, 1)[0]                              #this is the ID of the NCBI file
                print id
                query_taxonomy = get_taxonomy(id)

    aFile.close()
    openFlag = 0
    comparison_taxonomy = []
    with open(mgb_inputfile) as aFile:
        prevline = " "                        # string to store the line that shows the accession number of the hit
        total_blast_score = 0                      #variable to store the total MGB_score (i.e. the sum of the MGB_scores that are above user (e.g. 40.0)

        for num, line in enumerate(aFile, 1):
            if 'Details' in line:                                                                   # jump the file go to "Details" and then start reading from there again
                for num, line in enumerate(aFile, 1):                                               # num will hold the number of the line
                    m = re.match(r'\d+\.\s*.*', line)                                               #match the line that holds the accession number of the hit
                    if m:
                        prevline = re.match(r'(\d+)\.\s*(.*)', line).group()                        # If there is any match found then this would store the whole line to the variable prevline.
                        accession_no = prevline[3:]                                                 #from 3rd position till the end of the string. Slicing the string. Necessary to produce the ID of the hit to search in Entrez

                    searchObj = re.search(r'(\bMultiGeneBlast\b\s+\w+\:\s+(\d+.\d+))', line)        #search object to capture the line "MultiGeneBlast score:  "
                    if searchObj:
                        score = float(searchObj.group(2))                                           # store the MGB score float , e.g. 56.0
                        if score > float(ideal_mgb_score):                                          #check if the MGB score is bigger than the user defined MGB score
                            count = count + 1                                                       # count the number of times this comparison is true                   # eventually the low number of hits are
                            openFlag = 1

                            print "the number of times we found hits above" , ideal_mgb_score,  "are", count       #printing the count to see

                            print "the MGB score is", score

                            #print prevline                                                          #the line that holds the accession number of the hit
                            print accession_no                                                       # print the accession_no that met the criterion where the MGB_score is > user input
                            oldacc = accession_no                                                    # store to a different variable to convert the accession number to the id.
                            total_blast_score = total_blast_score + score                            #keep the addition of the MGB_scores that are individually bigger than the user input

                            sep = "_"
                            id_rest = oldacc.split(sep, 1)[0]
                            print ("the ID of this NCBI record is " + " " + (str(id_rest)))

                            hit_taxonomy = get_taxonomy(id_rest)

                            for i in query_taxonomy:
                                if i in hit_taxonomy:
                                    comparison_taxonomy.append(1)
                                else:
                                    comparison_taxonomy.append(0)

    aFile.close()


    tax_match = sum(comparison_taxonomy)

    if count>0:
        average_cluster_tax_distance = tax_match/count
    else:
        average_cluster_tax_distance = 100

    cluster_name = mgb_inputfile.split("/")      # splitting the cluster name (as it is given as the path to the file and splitting it every time there is a "."

    cluster_name2 = str(cluster_name[-1])        #take only the "clusterxxx" of each cluster


    print "the cluster name1 is", cluster_name2
    hit_dict[cluster_name2] = float(total_blast_score)            # and another dictionary of clustername : total MGB score
    print "the total blast score so far is", total_blast_score
    print "hit_dict is:",  hit_dict, '\n'

    count_hit_dict = {}                    # clustername : number of hits that scored a MGB score above the user defined
    count_hit_dict[cluster_name2] = count
    print "count_hit_dict is:",  count_hit_dict, '\n'

    tax_dict1 = {}

    tax_dict1[cluster_name2] = average_cluster_tax_distance
    print "this is the avergae tax distance for this cluster", tax_dict1

    return (cluster_name2,total_blast_score,get_values(count_hit_dict), get_values(tax_dict1))

#############################################################################################################################################################################################################################################
__author__ = 'Evanthia_Kaimaklioti'
#this is the script manipulating the MultiGeneBlast output files
#this pipeline was constructed for the purposes of the RP1 project within the MSc in Bioinformatics and Systems Biology
#course of the university of Manchester


from Bio import Entrez, SeqIO
import os
import re
import linecache
import operator
from operator import itemgetter
from collections import OrderedDict


def get_keys(dict1):
    for key, value in dict1.iteritems():
        return key

def get_values(dict):
    for key, value in dict.iteritems():
        return value

def get_taxonomy(NCBI_id):
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id=NCBI_id, retmode="text")     #handle to fetch the NCBI record with that ID, retmode: return mode
    record = SeqIO.read(handle, "genbank")                                   #handle to read the NCBI GenBank file
    for key in record.annotations.keys():                                    # go the annotations section of the GenBank file
        key = key.lower()                                                    #make all keys lowercase
        if key == 'taxonomy':
            string_taxonomy = record.annotations['taxonomy']
            return string_taxonomy

def multiGeneBlast_process(mgb_inputfile , ideal_mgb_score):
    Entrez.email = "evanthia.kaimaklioti@postgrad.manchester.ac.uk"

    count = 0
    string_taxonomy = []                    #dictionary that returns 1 when there is a match of the query taxonomy list to the one in hand and 0 when there is no match.
    line_wanted = ""

    accession_no_query = ""
    prevline1= ""
    hit_dict = {}
    query_taxonomy = []         # the list that hold the taxonomy information for the query

    with open(mgb_inputfile) as aFile:
        for line in aFile:
            if 'Significant' in line:                              # jump the file go to "Details" and then start reading from there again
                line_wanted = str(aFile.next())
                                                                 # but can I access this line_wanted from the loop below? The word significant is before the word Details, and the Details needs to be accessed first
                m = re.match(r'\d+\.\s*([a-zA-Z0-9_]+).*', line_wanted)         #match the line that holds the accession number of the hit
                accession_no_query = m.group(1)
                print "query accession is : ", accession_no_query
                sep = "_"                                     #this applies only if the separator is there. If it isn't it will return the whole of the oldacc which is O.K. also:
                id = accession_no_query.split(sep, 1)[0]                              #this is the ID of the NCBI file
                print id
                query_taxonomy = get_taxonomy(id)

    aFile.close()
    openFlag = 0
    comparison_taxonomy = []
    with open(mgb_inputfile) as aFile:
        prevline = " "                        # string to store the line that shows the accession number of the hit
        total_blast_score = 0                      #variable to store the total MGB_score (i.e. the sum of the MGB_scores that are above user (e.g. 40.0)

        for num, line in enumerate(aFile, 1):
            if 'Details' in line:                                                                   # jump the file go to "Details" and then start reading from there again
                for num, line in enumerate(aFile, 1):                                               # num will hold the number of the line
                    m = re.match(r'\d+\.\s*.*', line)                                               #match the line that holds the accession number of the hit
                    if m:
                        prevline = re.match(r'(\d+)\.\s*(.*)', line).group()                        # If there is any match found then this would store the whole line to the variable prevline.
                        accession_no = prevline[3:]                                                 #from 3rd position till the end of the string. Slicing the string. Necessary to produce the ID of the hit to search in Entrez

                    searchObj = re.search(r'(\bMultiGeneBlast\b\s+\w+\:\s+(\d+.\d+))', line)        #search object to capture the line "MultiGeneBlast score:  "
                    if searchObj:
                        score = float(searchObj.group(2))                                           # store the MGB score float , e.g. 56.0
                        if score > float(ideal_mgb_score):                                          #check if the MGB score is bigger than the user defined MGB score
                            count = count + 1                                                       # count the number of times this comparison is true                   # eventually the low number of hits are
                            openFlag = 1

                            print "the number of times we found hits above" , ideal_mgb_score,  "are", count       #printing the count to see

                            print "the MGB score is", score

                            #print prevline                                                          #the line that holds the accession number of the hit
                            print accession_no                                                       # print the accession_no that met the criterion where the MGB_score is > user input
                            oldacc = accession_no                                                    # store to a different variable to convert the accession number to the id.
                            total_blast_score = total_blast_score + score                            #keep the addition of the MGB_scores that are individually bigger than the user input

                            sep = "_"
                            id_rest = oldacc.split(sep, 1)[0]
                            print ("the ID of this NCBI record is " + " " + (str(id_rest)))

                            hit_taxonomy = get_taxonomy(id_rest)

                            for i in query_taxonomy:
                                if i in hit_taxonomy:
                                    comparison_taxonomy.append(1)
                                else:
                                    comparison_taxonomy.append(0)

    aFile.close()


    tax_match = sum(comparison_taxonomy)

    if count>0:
        average_cluster_tax_distance = tax_match/count
    else:
        average_cluster_tax_distance = 100

    cluster_name = mgb_inputfile.split("/")      # splitting the cluster name (as it is given as the path to the file and splitting it every time there is a "."

    cluster_name2 = str(cluster_name[-1])        #take only the "clusterxxx" of each cluster


    print "the cluster name1 is", cluster_name2
    hit_dict[cluster_name2] = float(total_blast_score)            # and another dictionary of clustername : total MGB score
    print "the total blast score so far is", total_blast_score
    print "hit_dict is:",  hit_dict, '\n'

    count_hit_dict = {}                    # clustername : number of hits that scored a MGB score above the user defined
    count_hit_dict[cluster_name2] = count
    print "count_hit_dict is:",  count_hit_dict, '\n'

    tax_dict1 = {}

    tax_dict1[cluster_name2] = average_cluster_tax_distance
    print "this is the avergae tax distance for this cluster", tax_dict1

    return (cluster_name2,total_blast_score,get_values(count_hit_dict), get_values(tax_dict1))
###############################################################################################################################################################################################################################################

__author__ = 'Evanthia_Kaimaklioti'
#this is the main programme of the SynPrio pipeline
#this pipeline was constructed for the purposes of the RP1 project within the MSc in Bioinformatics and Systems Biology
#course of the University of Manchester

import sys
import os
import csv
from os import path
from os import system
from operator import itemgetter
import operator
import argparse
from argparse import ArgumentParser
from os import walk
from os import listdir
from os.path import isfile, join
from pprint import pprint


#e.g. Run this code by typing in the command line:
# python SynPrio1.py -inputfile NC_003888.gbk -gcu 72 -wgc 0.2 -wl 0.13 -wcu 0.18 -il 12000 -wil 0.34
#  -mgb_inputfile cluster -mgb_score 50.0 -wms 0.21 -whc 0.24 -wtx 0.28


def get_asteric(string):					#function that checks whether there is an asteric in front of the weight value defined by the CLAs
    character = list(string)
    wt_string = ""
    dictionary_option = {}
    if character[0] == "*":                                 # if there is a star
       	order_flag1 = True				    # raise the ordering flag to be True
        just_weight1 = character[1:]			    # what remains by the string is the actual weight value
        wt_string = ''.join(just_weight1)			
    else:						    # otherwise if there is no star
        order_flag1 = False				    #raise the ordering flag to be False 
        wt_string = ''.join(character)				

    return order_flag1, float(wt_string)

def make_weighted_mean(weights):			   #function to compute the weighted means
    if sum(weights) == 0:			           # if the sum of ranks is 0 then make all the feature weights equal to 1.0
        weights = 1.0
    wt = sum(weights)					   # then return the sum of weights

    def wm(seq):					   
        return sum(w * v for w, v in zip(weights, seq))/wt    # return the weighted mean 
    return wm


def get_ranks(sd): 					# function that gets the rank of the clusters in each of the feature lists, which are alredy sorted
    k, val = sd[0]  
    r = 1
    rank = {k: r}
    delta = 1						# delta holds the number of times there is an equal ranking 

    for k, v in sd[1:]:
        if v == val:					# if there is an equal ranking
            delta += 1					
        else:
            val = v					# make the value in the dictionary equal to the value to compare everything else which 
            r += delta					# to the rank is added the delta
            delta = 1					# make the delta back to 1 
        rank[k] = r					#put the rank of the cluster in the dictionary

    return rank						


parser = argparse.ArgumentParser()
parser.add_argument('inputfile', nargs='?', default=argparse.SUPPRESS)             

parser.add_argument('-inputfile',  dest='inputfile', default='NO FILE SUPPLIED')

parser.add_argument('GC_usage', nargs='?', default=argparse.SUPPRESS)               #the ideal GC usage for the cluster command line argument (CLA)
parser.add_argument('-gcu', dest='GC_usage', default=63)

parser.add_argument('GC_weight', nargs='?', default=argparse.SUPPRESS)              #the weight applied to the GC usage
parser.add_argument('-wgc', dest='GC_weight', default=0)

parser.add_argument('weight_length', nargs='?', default=argparse.SUPPRESS)        # this is the weight for the ranking according to the lengthiest sequence in the cluster.
parser.add_argument('-wl', dest='weight_length', default=0)

parser.add_argument('codon_usage', nargs='?', default=argparse.SUPPRESS)       # this is the wt for the ranking according to the codon usage difference by E.coli of the cluster.
parser.add_argument('-wcu', dest='codon_usage', default=0)

parser.add_argument('ideal_length', nargs='?', default=argparse.SUPPRESS)       # the ideal length of the largest sequence of the cluster
parser.add_argument('-il', dest='ideal_length', default=1000)

parser.add_argument('weight_ideal_length', nargs='?', default=argparse.SUPPRESS)       # this is the wt for the ranking according to the length distance
parser.add_argument('-wil', dest='weight_ideal_length', default=0)

parser.add_argument('mgb_inputfile', nargs='?', default=argparse.SUPPRESS)          		# the name of the MGB folder containing the MGB output files
parser.add_argument('-mgb_inputfile',  dest='mgb_inputfile', default='NO FILE SUPPLIED')

parser.add_argument('mgb_score', nargs='?', default=argparse.SUPPRESS)              # CLA for the MGB score threshold above which we want to count scores.
parser.add_argument('-mgb_score', dest='mgb_score', default=20.0)

parser.add_argument('weight_mgb_score', nargs='?', default=argparse.SUPPRESS)       # this is the wt for the ranking according to the total MGB_score for the cluster
parser.add_argument('-wms', dest='weight_msg_score', default=0)

parser.add_argument('weight_hits_count', nargs='?', default=argparse.SUPPRESS)       # this is the wt for the ranking according to the count (number) of hits that scored a MGB_score above 
										                                            #the user defined MGB score (threshold)
parser.add_argument('-whc', dest='weight_hits_count', default=0)

parser.add_argument('weight_taxonomy', nargs='?', default=argparse.SUPPRESS)       # this is the wt for the ranking according to the average taxonomic distance of the hits that had the MGB_score
                                                                                   # #higher than the user defined
parser.add_argument('-wtx', dest='weight_taxonomy', default=0)

filename = parser.parse_args().inputfile					# the name of the folder containing the antiSMASH output files. It is named always by the name of the genome 
fileWithoutExt = filename.split(".")						#accession number that fed in the antiSMASH software. e.g. "NC_003888".

path_anti_folder = os.path.dirname(os.path.abspath(fileWithoutExt[0]))		# find the path to the antiSMASH output files folder (which is the current directory) 

antiSMASH_folder_path = path_anti_folder + "/" + fileWithoutExt[0] + "/"			# go to the path of the antiSMASH output file folder, i.e. access the folder.

onlyfiles = [f for f in listdir(antiSMASH_folder_path) if isfile(join(antiSMASH_folder_path, f))] #the list that holds the number of .gbk files in the folder

ideal_GC_usage = parser.parse_args().GC_usage                                #ideal GC_usage value
ideal_len = int(parser.parse_args().ideal_length)                            #ideal length of the largest sequence

wt_GC = parser.parse_args().GC_weight                                 # weight value for the GC_usage comparison
wt_GC1 = get_asteric(wt_GC)

wt_cuE = parser.parse_args().codon_usage                        #weight value for the codon usage difference from E.coli
wt_cuE1 = get_asteric(wt_cuE)

wt_il = parser.parse_args().weight_ideal_length         # weight for the the length distance from the user defined distance.
wt_il1 = get_asteric(wt_il)

wt_mgb = parser.parse_args().weight_msg_score                       #weight value for the total MGB score for the cluster
wt_mgb1 = get_asteric(wt_mgb)

wt_len = parser.parse_args().weight_length        #weight value for the size of the lengthiest sequence in the cluster
wt_len1 = get_asteric(wt_len)

wt_mhc = parser.parse_args().weight_hits_count            #weight value for the total number of hits above the threshold MGB score for the cluster
wt_mhc1 = get_asteric(wt_mhc)

wt_tax = parser.parse_args().weight_taxonomy	# weight for the taxonomy distance from the total taxonomic distance (TD) of the hits with a MGB score above the threshold. 
wt_tax1 = get_asteric(wt_tax)

mgb_folder = parser.parse_args().mgb_inputfile		
ideal_mgb_score = parser.parse_args().mgb_score		# ideal total MGB score 


weights = [float(wt_GC1[1]),float(wt_len1[1]), float(wt_cuE1[1]), float(wt_il1[1]),float(wt_tax1[1]),float(wt_mgb1[1]),float(wt_mhc1[1])]  #put the weights for the seven features in a list 

#if float(sum(weights)) == 0:					#if the sum of weights is equal to 0
#    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]		# then put 1 for the weight values of the cluster features


inputstring = path_anti_folder + "/" + fileWithoutExt[0]+"/" + fileWithoutExt[0]+".3.cluster"   # the path to the antiSMASH output cluster files. Reasoning for the concatenation: because the files 													#provided by antiSMASH have NC_003888.3.clusterxxx.gbk name format

outputstring = path_anti_folder + "/" + "fasta_files" + "/" + fileWithoutExt[0]+ ".3.cluster "  #the path to create the files created after the antiSMASH output handling- creates in fasta_files  													#folder and the files are named automatically retaining their parent cluster name.
outputstring1 = path_anti_folder + "/cds_files/" + fileWithoutExt[0]+ ".3.cluster "		#save the synonymous codon usage of all sequences in the cluster in a separate "cds_file" folder	
												# retain the name of the input antiSMASH cluster file 


import antiSMASH_processing1					# run the antiSMASH_processing.py script through the SynPrio.py script

GC_distance = {}						# declare the different dictionaries 
sizeDict = {}
codon_usaged = {}
sizeDistDict = {}

for i in range(1, (len(onlyfiles)+1)):				#for all the .gbk files in the antiSMASH output files containing folder, take them one at a time
   #create the path for some of the input arguments of the antiSMASH_processing.py
    inputFilename = inputstring + '%0.3d'%i + '.gbk'                               #NC_003888.3.clusterxxx.gbk
    outputFilename = outputstring + '%0.3d'%i + '_cds.fasta'                       #NC_003888.3.clusterxxx_cds.fasta
    outputFilename2 = outputstring1 + '%0.3d'%i + 'codonFreqs.txt'                 #NC_003888.3.clusterxxxcodonFreqs.txt

    x, y, z, xi, xii = antiSMASH_processing1.geneProcess(inputFilename, outputFilename, outputFilename2, ideal_GC_usage,ideal_len) # LHS: output variables of antiSMASH_processing; RHS: input arguments 															          #of antiSMASH_processing.py

    GC_distance[str(x)] = float(y)                                               # a dictionary to hold the clusterxxx : GC distance_from_ideal value.
    sorted_GC_distance = sorted(GC_distance.items(), key=operator.itemgetter(1), reverse = wt_GC1[0])   #sort the dictionary according to the CLA-whether the order flag is True/False, if True: sort in 
													# descending order for the GC_distance feature

    sizeDict[x] = z                                                                         #dictionary of length of the largest sequence 
    sorted_lenDict = sorted(sizeDict.items(), key=operator.itemgetter(1), reverse=wt_len1[0]) # sort the dictionary according to to the CLA-whether the order flag is True/False, if False: sort in 
												# ascending order for the largest sequence's length feature

    codon_usaged[x] = xi                         # dictionary to store the synonymous codon usage (SCU) distance from the SCU of E.coli
    sorted_codon_usage = sorted(codon_usaged.items(), key=operator.itemgetter(1), reverse = wt_cuE1[0])   #sorting according to ... (please see above)

    sizeDistDict[str(x)] = int(xii)              #dictionary that stores the size distance from the user input ideal max length
    sorted_sizeDictDict = sorted(sizeDistDict.items(), key=operator.itemgetter(1), reverse = wt_il1[0])


import mgb_processing


path_mgb_folder = os.path.dirname(os.path.abspath(mgb_folder))		# find the path to the folder containing the MGB output files. The folder is in the current directory 

mgb_folder_path = path_mgb_folder + "/" + mgb_folder + "/"	# go to the path of the MGB output file folder, i.e. access the folder.

onlyfiles1 = [f for f in listdir(mgb_folder_path) if isfile(join(mgb_folder_path, f))] #the list that holds the number of MGB output files in the folder

inputstring_mgb = path_mgb_folder + "/" + mgb_folder + "/" + mgb_folder


mgb_Dict = {}			# declare the dictionaries that hold the distance of a cluster's feature value from the ideal user-defined value of the feature
count_hits_Dict = {}
tax_distance_Dict = {}
sorted_mgb_Dict = {}
sorted_count_hits_Dict={}
sorted_mgb_Dict1 = {}

for index in range(1, (len(onlyfiles1)+1)):
    mgb_inputfile = inputstring_mgb + '%0.3d'%index         # access the MGB output cluster files which have a name numbering of 3 digits, e.g. cluster001

    clust,sum_score, counthits,tax_distance = mgb_processing.multiGeneBlast_process(mgb_inputfile, ideal_mgb_score)

    mgb_Dict[clust] = sum_score				    
    sorted_mgb_Dict = sorted(mgb_Dict.items(), key=operator.itemgetter(1), reverse= wt_mgb1[0])    

    count_hits_Dict[clust] = counthits
    sorted_count_hits_Dict = sorted(count_hits_Dict.items(), key=operator.itemgetter(1), reverse= wt_mhc1[0])    #so that the clusters with the highest number of hits with  (total) cumulative 														 # MGB_score are ranked highest

    tax_distance_Dict[clust] = tax_distance
    sorted_tax_distance_Dict = sorted(tax_distance_Dict.items(), key=operator.itemgetter(1), reverse = wt_tax1[0])            



weighted_mean = make_weighted_mean(weights)

sorted_dicts = [sorted_GC_distance,sorted_lenDict, sorted_codon_usage, sorted_sizeDictDict, sorted_tax_distance_Dict, sorted_mgb_Dict, sorted_count_hits_Dict]      #the order of the weights on weights 																				#variable is the same as that on the sorted_dicts

print "Sorted dicts:sorted_GC_distance, sorted_lenDict, sorted_codon_usage, sorted_sizeDictDict, sorted_tax_distance_dictm sorted_mgb_dict, sorted_count_hits_dict"
pprint(sorted_dicts, indent=4)

all_ranks=[get_ranks(sd) for sd in sorted_dicts]
print '\nAll ranks: sorted_GC_distance, sorted_lenDict, sorted_codon_usage, sorted_sizeDictDict, sorted_tax_distance, sorted_mgb_Dict, sorted_count_hits_Dict'
pprint(all_ranks, indent=4)

keys = sorted(k for k, v in sorted_dicts[0])

means = [(k, weighted_mean([ranks[k] for ranks in all_ranks])) for k in keys]
print "\nWeighted means:"
pprint(means, indent=4)

sorted_means = sorted(means, key=lambda tup: tup[1], reverse=True)   #sort the ranking list of the clusters according to their weighted mean of ranks score (WMRS) in ascending order
print "\nSorted weighted means:"
pprint(sorted_means, indent=4)

with open('weightedMeans.txt', 'w') as file:			 # save this cluster ranking list to a tab delimited text in the current directory 
    for key, value in means:
        file.writelines([key])
        file.writelines('\t')
        file.writelines(str(value))
        file.writelines('\n')










