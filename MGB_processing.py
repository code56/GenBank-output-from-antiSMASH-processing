
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

