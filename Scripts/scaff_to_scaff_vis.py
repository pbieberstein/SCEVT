import pandas
import pickle
import re
import matplotlib.pyplot as plt
import random
from Bio import SeqIO



class ScaffoldBNG:
    '''
    This is a class that has all the relevant information for each scaffold.
    '''
    def __init__(self, seqIO_record):
        self.id = seqIO_record.id # Scaffold Name from fasta file
        self.length = len(seqIO_record.seq) # Scaffold length
        self.sequence = seqIO_record.seq
        self.gaps = self.N_region_finder() # this will be a dictionary with all the gaps! {1:[4,200], 2:[455, 7676]}
        self.record = seqIO_record

    def __repr__(self):
        return "<Scaffold Object> with ID: " + self.id

    def N_region_finder(self):
        '''
        Find the regions of NNNNs in the sequence and records the indices of starting and ending points in two lists.
        Return a tuple of two lists. If the sequence starts of known BasePairs, then the 'N-region_end' list will have index 0.
        :param seq_record:
        :return: returns a tuple of two lists. (N_region_start, N_region_end) which each contain the coordinates of either starting or
        ending 'N' regions (gaps)
        '''
        # Initialize N_region_start & N_region_end lists
        # Initialize the vertical position of the pointer for plotting (move ~50 px upwards)
        sequence = str(self.sequence).upper()
        N_region = False
        N_region_start = []
        N_region_end = []
        gaps = {}
        index = 0

        for letter in sequence:
            # Base Case: Check if N or other
            # If N, set N_region to True & include index (0) as marker in "N_region_start" list
            # Else, set N region to False & include index (0) as marker in "N_region_end" list
            if index == 0:
                if letter == "N":
                    N_region = True
                    N_region_start.append(0)
                else:
                    N_region = False
                index += 1
            # Then for all other bases, check if the same as before
            # If theres a change in "N_region":
            # include current index as marker in N_region_start if the N_region is now True
            else:
                if letter == "N" and N_region == False:
                    N_region_start.append(index)
                    N_region = True
                # Else include current index as marker in N_region_end if the N_region is now False
                elif letter != "N" and N_region == True:
                    N_region_end.append(index)
                    N_region = False
                index += 1
        for i in range(len(N_region_start)):
            gaps[i] = [N_region_start.pop(0), N_region_end.pop(0)]
        return gaps


def generate_scaffold_objects(list_of_scaffold_names, path_to_input_fasta):
    dict_of_scaffold_objects = {}
    for seq_record in SeqIO.parse(path_to_input_fasta, "fasta"):
        if seq_record.id in list_of_scaffold_names:
            index = list_of_scaffold_names.index(seq_record.id) # get the index of the scaffold ID in the list
            list_of_scaffold_names.pop(index) # Delete that item from the scaffold list, so at the end we know what we didn't find!
            scaffold_object = ScaffoldBNG(seq_record)
            dict_of_scaffold_objects[scaffold_object.id] = scaffold_object # Add the seq Record to the dictionary

    for item in list_of_scaffold_names: # Go through the original list and print out all the scaffold names that weren't found in the fasta file...
        print item, " Was not found in the sequence file"
    return dict_of_scaffold_objects


def create_temp_fasta_files(list_of_scaffold_names,path_to_input_fasta, output_file_name):
    '''
    This function extracts the sequences from specific scaffolds from a big fasta file (containing the entire genome)
    :param list_of_scaffold_names:
    :return:
    '''
    with open(output_file_name, "w") as temp_file:
        for seq_record in SeqIO.parse(path_to_input_fasta, "fasta"):
            if seq_record.id in list_of_scaffold_names:
                index = list_of_scaffold_names.index(seq_record.id) # get the index of the scaffold ID in the list
                list_of_scaffold_names.pop(index) # Delete that item from the scaffold list, so at the end we know what we didn't find!
                SeqIO.write(seq_record, temp_file, "fasta")
    for item in list_of_scaffold_names:
        print item, " Was not found in the sequence file"
    temp_file.close()
    return None





region_fasta_file = "60444_07april_cmd2_region.fasta"

# TME3 cmd2 scaffold list
cmd2_scaffolds_list = ['Super-Scaffold_44', 'Super-Scaffold_1951', 'Super-Scaffold_730', 'Super-Scaffold_1022', '001839F', '003097F', '002893F', '006625F', '007880F']

# FOR 60444
#cmd2_scaffolds_list = ['Super-Scaffold_111','Super-Scaffold_29', 'Super-Scaffold_885','006018F','Super-Scaffold_526','001529F','004409F','001787F','004969F','008344F','Super-Scaffold_572']

#create_temp_fasta_files(list_of_scaffold_names=cmd2_scaffolds_list, path_to_input_fasta="../sampleData/Genome_sequences/TME3_draft.fasta", output_file_name="../temp/60444_sub_selection.fasta")

dictionary =  generate_scaffold_objects(list_of_scaffold_names=cmd2_scaffolds_list, path_to_input_fasta="../sampleData/Genome_sequences/TME3_draft.fasta")

print dictionary['Super-Scaffold_44'].gaps