
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
