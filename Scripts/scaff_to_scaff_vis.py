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
    def __init__(self, seqIO_record, psl_file):
        self.id = seqIO_record.id # Scaffold Name from fasta file
        self.length = len(seqIO_record.seq) # Scaffold length
        self.sequence = seqIO_record.seq
        self.gaps = self.N_region_finder() # this will be a dictionary with all the gaps! {1:[4,200], 2:[455, 7676]}
        self.record = seqIO_record
        self.genes = self.get_gene_coordinates(psl_filename = psl_file) # This gets the the gene or SNP mappings (coordinates) and ties them to this scaffold object


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

    def get_gene_coordinates(self, psl_filename, quality_of_mapping= 0.95):
        '''
        The .psl file should look like this:
        matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName   qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts
        from this we want to extract all the gene names & their coordinates that lie on this specific scaffold

            The Gene name is in column: 9
        The scaffold that the gene was mapped to is in column: 13
        The coordinate where the match STARTS within the scaffold/contig is in column: 15
        The coordinate where the match ENDS within the scaffold/contig is in column: 16

        :param type:
        :param psl_file:
        :return:
        '''

        gene_dict = {}

        # This function should go through the gene mapping file, and extract all the gene names and their coordinates that lie on this scaffold!
        # Filter the .psl file so it only includes CMD2 scaffold information
        psl_gene_mappings = pandas.read_csv(psl_filename, sep='\t',
                                            header=None)  # load in the gene mapping .psl file
        psl_gene_mappings = psl_gene_mappings[
            psl_gene_mappings[13] == self.id]  # Filter so we only keep rows of relevant cmd2 scaffolds
        # Filtering out genes that don't have at least 90% matching bps of the mRNA
        psl_gene_mappings = psl_gene_mappings[psl_gene_mappings[0] > quality_of_mapping * (psl_gene_mappings[10])]

        for index, row in psl_gene_mappings.iterrows():
                if row[13] == self.id:
                    coordinate_range = [row[15], row[16]]
                    gene_name = row[9]
                    if gene_name in gene_dict:
                        gene_dict[gene_name].append(coordinate_range)
                    else:
                        gene_dict[gene_name] = [coordinate_range]
        return gene_dict


    def invert_scaffold(self):
        half = int(self.length/2)
        for key in self.gaps: # invert the gap coordinates
            gaps[key] = mirror_lists(half, gaps[key]) # give it the midpoint of the object and the list of coordinates it should invert on that object

        new_gene_coordinates = {}
        for key in self.genes: # {"gene1:[[22,532],[2343,64545]],"gene2:[[2,44]]...}
            for instance in genes[key]:
                if key in new_gene_coordinates:
                    new_gene_coordinates[key].append(mirror_lists(instance))
                else:
                    new_gene_coordinates[key] = [mirror_lists(instance)]
        self.genes = new_gene_coordinates
        print self.id + " was Inverted."
        return self.genes


    def __repr__(self):
        return "<Scaffold Object> with ID: " + self.id


def generate_scaffold_objects(list_of_scaffold_names, path_to_input_fasta, psl_file):
    dict_of_scaffold_objects = {}
    for seq_record in SeqIO.parse(path_to_input_fasta, "fasta"):
        if seq_record.id in list_of_scaffold_names:
            index = list_of_scaffold_names.index(seq_record.id) # get the index of the scaffold ID in the list
            list_of_scaffold_names.pop(index) # Delete that item from the scaffold list, so at the end we know what we didn't find!
            scaffold_object = ScaffoldBNG(seq_record, psl_file)
            dict_of_scaffold_objects[scaffold_object.id] = scaffold_object # Add the seq Record to the dictionary

    for item in list_of_scaffold_names: # Go through the original list and print out all the scaffold names that weren't found in the fasta file...
        print item, " Was not found in the sequence file"
    return dict_of_scaffold_objects


def mirror_lists(mid_point_coord, list_to_invert):
    # This function takes a list of coordinates and the coordinate of the halfway points
    # it then goes through the list and inverts it by:
    # 1. Move it to the left by half
    # 2. Multiply all coordinates by -1
    # 3. move to the right by half
    if isinstance(list_to_invert, list):
        for i, item in enumerate(list_to_invert):
            result = ((item - mid_point_coord)*-1)+mid_point_coord
            list_to_invert[i] = result

    elif isinstance(list_to_invert, int):
            item = list_to_invert
            result = ((item - mid_point_coord) * -1) + mid_point_coord
            list_to_invert = result

    return list_to_invert

def create_subset_fasta_files(list_of_scaffold_names,path_to_input_fasta, output_file_name): # NOT NEEDED ANYMORE
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


def create_temp_psl_mapping_file(psl_gene_mapping_file, scaffolds_list, output_prefix="temp", quality_of_mapping="0.95"): # NOT NEEDED ANYMORE
    '''
    This function pre-processes our initially large .psl mapping file (From BLAT output) (because they contain information about all genes on the entire genome)
    We only want rows with genes that mapped to our cmd2 specific scaffolds
    :param psl_gene_mapping_file:
    :param scaffolds_list: List of cmd2 related scaffolds
    :return: The filename of the output file
    '''
    # Filter the .psl file so it only includes CMD2 scaffold information
    psl_gene_mappings = pandas.read_csv(psl_gene_mapping_file, sep='\t', header=None) # load in the gene mapping .psl file
    psl_gene_mappings = psl_gene_mappings[psl_gene_mappings[13].isin(scaffolds_list)] # Filter so we only keep rows of relevant cmd2 scaffolds
    # Filtering out genes that don't have at least 90% matching bps of the mRNA
    psl_gene_mappings = psl_gene_mappings[psl_gene_mappings[0] > quality_of_mapping*(psl_gene_mappings[10])] # The match must include at least <quality_of_mapping>% of the total basepairs in that gene
    output_filename = output_prefix+"_filtered_gene_mapping.pickle"
    psl_gene_mappings.to_pickle(output_filename) # save it as a temporary pickle python file!
    return output_filename

def create_temp_reference_gff_file(ref_genome_gff_file, output_prefix="temp"): # NOT NEEDED ANYMORE
    '''
    This function pre-processes our initially large reference genome file Because we only want gene information from the reference file, not cDNA, mRNA, etc.

    Then the ref_gene_information file is the Cassava Reference Genome v6.1 file from Phytozome (Mesculenta_305_v6.1.gene.gff3
    With the following columns:
    Start coordinate in reference genome column: 3
    End coordinate in reference genome column: 4
    Gene_ID (includes other junk) column: 8
    ** We only care about rows that have "gene" in column 2, otherwise it also includes mRNA, CDS, etc. information - we only need gene coordinates


    :param ref_genome_gff_file: downloaded from phytozome (tested with Cassava v6.1)
    :return: the filename of the output file
    '''
    # and filter the .gff3 file so it includes only rows with gene coordinates (Filter out mRNA, cDNA, and other rows objects)
    ref_gene_info = pandas.read_csv(ref_genome_gff_file, sep='\t', header=None) # load in the reference gene information file
    ref_gene_info = ref_gene_info[ref_gene_info[2] == "gene"] # Filter so we only keep rows of gene coordinate information (not cDNA, RNA, etc...)
    output_filename = output_prefix+"_reference_info.pickle"
    ref_gene_info.to_pickle(output_filename) # save it as a temporary pickle python file!
    return output_filename


#####################################
## Running the script
#####################################

## PARAMTERS


ref_fasta_file = "TME3_07april_cmd2_region.fasta"

target_fasta_file = "TME3_07april_cmd2_region.fasta"

genome_fasta_file = "../sampleData/Genome_sequences/TME3_draft.fasta"

ref_psl_file = "../sampleData/Gene_BLAT_mappings/TME3_BNG_plus_notscaff.psl"

target_psl_file = "../sampleData/Gene_BLAT_mappings/TME3_BNG_plus_notscaff.psl"

# TME3 cmd2 scaffold list
ref_scaffold_list = ['Super-Scaffold_44'] #, 'Super-Scaffold_1951', 'Super-Scaffold_730', 'Super-Scaffold_1022', '001839F', '003097F', '002893F', '006625F', '007880F']

target_scaffold_list = ["Super-Scaffold_1951"]


ref_scaff = generate_scaffold_objects(list_of_scaffold_names=ref_scaffold_list, path_to_input_fasta=genome_fasta_file, psl_file=ref_psl_file)

target_scaff = generate_scaffold_objects(list_of_scaffold_names=target_scaffold_list, path_to_input_fasta=genome_fasta_file, psl_file=target_psl_file)

print ref_scaff


print target_scaff