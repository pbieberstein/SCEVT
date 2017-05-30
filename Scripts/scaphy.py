import pandas
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import re
from Bio import SeqIO
from numpy import sign, mean
import os
from datetime import datetime
'''
This script was written using Python 2.7
You'll need BioPython and matplotlib installed on your machine to make this script work
This is easiest done by downloading either miniconda or the entire anaconda package


To run this script, edit the part below this (Parameters) to what you like. Then execute this script by typing:

python scaphy.py

The output will be saved in the same directory as <timestamp>_scaphy_output.pdf

Written by:
Philipp v. Bieberstein
rphilipp@ethz.student.ch

'''

#####################################################################################################################
#########                                      || Parameters ||                                           ###########
#####################################################################################################################

## Fasta files were the scaffold sequences are included (Can be entire genome files, or just fasta files that contain the scaffold sequence)
# Needed to draw scaffolds with their gaps

reference_fasta_file = "../Data/Genome_sequences/TME3_draft.fasta"

## BLAT output files (.psl) that includes gene mappings onto the genome which includes the specified scaffolds
# The paths should be relative to where this script is... or absolute paths
ref_psl_file = "../Data/Gene_BLAT_mappings/TME3_BNG_plus_notscaff.psl"

# Reference genome information file (.gff3) from Phytozome (I deleted the first 3 header rows for simplicity)
reference_gff_file = "../Data/Reference_genome_info/Mesculenta_305_v6.1.gene.gff3"

# Parameters for Plotting

N_region_min = 100  # threshold how big NNNN regions have to be in order to plot them


# TME3 cmd2 scaffold list
# For example: 'Super-Scaffold_1951' or 'Super-Scaffold_730'

list_of_scaffolds = ['Super-Scaffold_1022', 'Super-Scaffold_44', '000979F', '002893F', 'Super-Scaffold_1206', '006535F', '004561F', 'Super-Scaffold_12475']

# Whether you want to invert the reference or target scaffold:
# simply list all scaffold names that you want inverted like: ['scaffold1', 'scaffold2', 'scaffold3']

scaffolds_to_invert = ['Super-Scaffold_44']

## Manual Plotting Settings
# For each scaffold, you need to specifiy a dictionary to declare the x_offset value, vertical_pointer value, the color of the lines
## colors:
# b: blue
# g: green
# r: red
# c: cyan
# m: magenta
# y: yellow
# k: black
# w: white

# If you want automatic plotting, simply comment this dictionary with the block comments: '''...'''
# if you want manual plotting, un-comment this and specify the x_offset, vertical_pointer and line_color for each scaffold
'''
manual_plotting = {'Super-Scaffold_1022': {"x_offset":1000000,
                                           "vertical_pointer": -100,
                                           "line_color": 'b'},
                   'Super-Scaffold_44': {"x_offset":1000000,
                                         "vertical_pointer": 100,
                                         "line_color": 'r'},
                   '000979F': {"x_offset":2000000,
                               "vertical_pointer": -200,
                               "line_color": 'b'}
                   }
'''

## Size of Reference Chromosome on Plot:
# pseudo_length = 31578400 # estimated chromosome 12 length from phytozome Cassava v6.1
reference_chromosome_length = 11000000  # only part of it so we can see everything better, but if genes lie beyond this region, then make this longer


## What reference chromosome to focus on (this needs be exact same string as is used in the reference annotation file)
reference_chromosome = 'Chromosome12'


#####################################################################################################################
#########                                   || End of Parameters ||                                       ###########
#####################################################################################################################





class OverallPlotter:
    def __init__(self, scaffold_dictionary, reference_vertical_pointer=0):
        self.vertical_pointer = 100
        self.line_color = 'k'
        self.current_chr_gene_color = 'y'
        self.remote_gene_color = 'c'
        self.scaffolds = scaffold_dictionary
        # Dictionary looks like: {'scaffold1':scaffold_object, 'scaffold2':scaffold_object}
        # each scaffold object has the gene_mapping dictionary like this:
        # {1: {'scaffold_coordinate': 229627, 'reference_chromosome': 'Chromosome12', 'reference_coordinate': 7495447.0, 'name': 'Manes.12G073900'}, 2:...}}

        self.reference_vertical_pointer = reference_vertical_pointer

    def auto_plot(self):
        '''
        This method should create the entire plot automatically without looking for manual tweaks
        It will first plot the reference chromosome
        Then plot the scaffolds
        Label the scaffolds with scaffold names and also gene hit proportions
        Then plot the genes on the scaffolds:
            - Unique genes in 1 color
            - mutual genes in another color and with connects to the reference chromosome with marks in it as well

        '''
        self.plot_reference_chromosome()
        self.plot_scaffolds_and_genes()

        return None

    def plot_reference_chromosome(self):
        # the reference chromosome will by default sit on y=0 but it can also be manually adjusted
        plt.plot([0, reference_chromosome_length], [self.reference_vertical_pointer, self.reference_vertical_pointer],
                 color='c', linestyle='-', linewidth=4)

        return None

    def plot_scaffolds_and_genes(self):
        # draw the objects with the given parameters
        for scaffold in self.scaffolds:
            object = self.scaffolds[scaffold]
            object.plot_scaffold(vertical_pointer=self.vertical_pointer, x_offset=object.x_offset)
            object.plot_all_genes(vertical_pointer=self.vertical_pointer, x_offset=object.x_offset, mutual_color=self.current_chr_gene_color , unique_color=self.remote_gene_color, line_color=object.line_color)
            self.vertical_pointer = (abs(self.vertical_pointer) + 50) * (sign(self.vertical_pointer) * -1)

        return None

    def manual_scaffold_plotter(self, manual_plotting):
        '''
        This script plots the scaffolds with manual arrangments
        :param manual_plotting:
        :return:
        '''
        # here all variables will have to be defined, and then it just draws the scaffolds exactly how we want it to
        # with vertical pointer and x_offsets set
        self.plot_reference_chromosome()
        # Give manual Plotting attributes
        for scaffold in manual_plotting: # only iterate trough the ones we have specified in manual dictionary
            object = self.scaffolds[scaffold]

            # saves the manual attributes to the objects themselves
            object.get_manual_attributes(manual_plotting)

            object.plot_scaffold(vertical_pointer=object.vertical_pointer, x_offset=object.x_offset)
            object.plot_all_genes(vertical_pointer=object.vertical_pointer, x_offset=object.x_offset, mutual_color=self.current_chr_gene_color , unique_color=self.remote_gene_color, line_color=object.line_color)

        return None



        return None

    def color_chooser(self):
        # set line_color to something new, but also increment so next call it'll be something new
        return None


class ScaffoldBNG:
    '''
    This is a class that has all the relevant information for each scaffold.

    Very importantly, it has a method that connects the reference gene location information to the scaffold gene locations...
    '''

    def __init__(self, seqIO_record, psl_file, line_color):
        self.id = seqIO_record.id  # Scaffold Name from fasta file
        self.length = len(seqIO_record.seq)  # Scaffold length
        self.sequence = seqIO_record.seq
        self.gaps = self.N_region_finder()  # this will be a dictionary with all the gaps! {1:[4,200], 2:[455, 7676]}
        self.record = seqIO_record
        self.genes = self.get_gene_coordinates(
            psl_filename=psl_file)  # This gets the the gene or SNP mappings (coordinates) and ties them to this scaffold object
        # self.genes => {'Manes.12G061900.1': [[1404095, 1405697]], 'Manes.12G059000.1': [[2613631, 2614138]]...}
        self.gene_mappings = None
        # {1: {'scaffold_coordinate': 229627, 'reference_chromosome': 'Chromosome12', 'reference_coordinate': 7495447.0, 'name': 'Manes.12G073900'}, 2:...}}
        self.line_color = line_color  # this specifies the line color for the gene matches connectors


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

    def get_gene_coordinates(self, psl_filename, quality_of_mapping=0.95):
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

        pattern = re.compile('.*.(?=.\d)')  # Pattern matching to gene underlying gene name instead of CDS names...

        for index, row in psl_gene_mappings.iterrows():
            if row[13] == self.id:
                coordinate_range = [row[15], row[16]]
                gene_name = row[9]
                gene_name = re.match(pattern,
                                     gene_name).group()  # This should erase the CDS specifications at the end of the gene names ".3"

                almost_same = False
                if gene_name in gene_dict:
                    for coordinates in gene_dict[gene_name]:  # coordinates will be a list of start & end [start, end]
                        # gene exists already in the gene dictionary, so we have to check that the match is different enough for us to care.. (not splicing)
                        if abs(coordinates[0] - coordinate_range[0]) < 300 or abs(
                                        coordinates[1] - coordinate_range[1]) < 300:
                            almost_same = True  # probably just difference in splicing, so we can ignore it
                        else:
                            None
                    if almost_same:
                        None
                    else:  # coordinates are not just different splicing, so we add them to the dictionary
                        gene_dict[gene_name].append(coordinate_range)
                else:  # we have no coordinates for this gene yet, so we add it freshly
                    gene_dict[gene_name] = [coordinate_range]  # that gene doesn't exist yet, so I add it freshly
        return gene_dict

    def calculate_automatic_offset_and_mapping_stats(self):
        # calculate averge gene location on the reference genome and move this scaffold by that much
        avg_coordinates = 0
        num_on_ref_chr = 0
        num_elsewhere = 0
        temp_for_avg = []

        for match in self.gene_mappings:
            gene_match = self.gene_mappings[match]
            if gene_match['reference_chromosome'] == reference_chromosome:
                temp_for_avg.append(gene_match['reference_coordinate'])  # add all the ref locations to the list
                num_on_ref_chr += 1
            else:
                num_elsewhere += 1
        if len(temp_for_avg) >= 1:
            avg_coordinates = mean(temp_for_avg)-(self.length/2) # because if there is no coordinate... then its an error
        else:
            print("No gene was on this, so x_offset will be 0")
            avg_coordinates = 0

        self.x_offset = avg_coordinates
        self.reference_chromosome_hits = num_on_ref_chr  # this is how many genes lie on the specified reference chromosome
        self.other_chromosome_hits = num_elsewhere  # this is how many genes lie on some other chromosome or another contig
        return None

    def plot_scaffold(self, vertical_pointer=100, x_offset=0):
        print("Plotting: " + self.id + "...")
        plt.plot([0 + x_offset, self.length + x_offset], [vertical_pointer, vertical_pointer], color='k',
                 linestyle='-', linewidth=4)
        text_location = vertical_pointer + 20 + len(list_of_scaffolds)
        # drawing the gaps
        for key in self.gaps:
            gap = self.gaps[key]  # returns a list witht the start of the gap and the end of that same gap
            start = gap[0]
            end = gap[1]
            if abs(end - start) > N_region_min:  # checks if gap is big enough for us to care about drawing it
                offsetter = -12  # silly hack to make the gap drawing fine enough for small gaps... otherwise the giant red blots just overlap way too much
                for i in range(1, 48):
                    offsetter += .5
                    plt.plot([start + x_offset, end + x_offset],
                             [vertical_pointer + offsetter, vertical_pointer + offsetter], color='r', linestyle='-',
                             linewidth=0.1, alpha=1)

        # Now add the labeling
        # the scaffold name
        plt.text(0 + x_offset, text_location, self.id, fontsize=5)

        text_location = vertical_pointer - 5
        # The gene location stats
        stats_text = str(self.reference_chromosome_hits) + "/" + str(
            self.reference_chromosome_hits + self.other_chromosome_hits) + " Genes Mapped Here"
        plt.text(100000 + x_offset + self.length, text_location, stats_text, fontsize=5)

    def plot_all_genes(self, vertical_pointer=100, x_offset=0, unique_color='c', mutual_color='y', line_color='c'):
        print("Plotting all genes for: ", self.id)
        # access all the genes, for each that is on the right chromosome, draw the connection, if not, just draw it with the other color
        for match in self.gene_mappings:
            gene_dictionary = self.gene_mappings[match]
            ref_chromosome = gene_dictionary[
                'reference_chromosome']  # Accesses that match's chromosome name where the gene is located
            gene_id = gene_dictionary['name']  # Accesses the gene name
            scaffold_coordinate = gene_dictionary['scaffold_coordinate']  # gets the coordinate for the scaffold
            reference_coordinate = gene_dictionary[
                'reference_coordinate']  # gets the coordinate for the reference chromosome
            # If the gene is on chromosome 12, draw it on the scaffold, on the reference, and a line between them
            if ref_chromosome == reference_chromosome:
                plot_thin_marker(scaffold_coordinate + x_offset, vertical_pointer,
                                 color=mutual_color)  # draw it in green
                plot_thin_marker(reference_coordinate, 0, color=mutual_color)  # draw it in green
                plt.plot([scaffold_coordinate + x_offset, reference_coordinate], [vertical_pointer, 0], linestyle='-',
                         color=line_color,
                         linewidth=0.1, alpha=0.8)
            else:
                plot_thin_marker(scaffold_coordinate + x_offset, vertical_pointer,
                                 color=unique_color)  # draw it in green
        return None

    def invert_scaffold(self):
        '''
        This function aims to invert this scaffold including its NNN regions, and its gene locations
        :return:
        '''
        half = int(self.length / 2)
        for key in self.gaps:  # invert the gap coordinates
            self.gaps[key] = self.mirror_lists(half, self.gaps[
                key])  # give it the midpoint of the object and the list of coordinates it should invert on that object

        new_gene_coordinates = {}
        for key in self.genes:  # {"gene1:[[22,532],[2343,64545]],"gene2:[[2,44]]...}
            for instance in self.genes[key]:
                if key in new_gene_coordinates:
                    new_gene_coordinates[key].append(self.mirror_lists(half, instance))
                else:
                    new_gene_coordinates[key] = [self.mirror_lists(half, instance)]
        self.genes = new_gene_coordinates

        for instance in self.gene_mappings:
            match_dictionary = self.gene_mappings[instance]
            self.gene_mappings[instance]['scaffold_coordinate'] = self.mirror_lists(half, match_dictionary[
                'scaffold_coordinate'])

        print(self.id + " was Inverted.")
        return self.genes

    def mirror_lists(self, mid_point_coord, list_to_invert):
        # Helper function for the invert method
        # This function takes a list of coordinates and the coordinate of the halfway points
        # it then goes through the list and inverts it by:
        # 1. Move it to the left by half
        # 2. Multiply all coordinates by -1
        # 3. move to the right by half
        if isinstance(list_to_invert, list):
            for i, item in enumerate(list_to_invert):
                result = ((item - mid_point_coord) * -1) + mid_point_coord
                list_to_invert[i] = result

        elif isinstance(list_to_invert, int):
            item = list_to_invert
            result = ((item - mid_point_coord) * -1) + mid_point_coord
            list_to_invert = result

        return list_to_invert

    def get_reference_mappings(self, reference_gff_file):
        # Pre-filter the phytozome reference genome gff file so that the following operations are faster...
        filtered_gff = pre_filter_phytozome_gff_file(
            reference_gff_file)  # this is now a pandas object with all the reference genome information

        # For this scaffold object, go through all the genes (.genes attribute) and find the location of them in the reference genome
        # Then we generate a new attribute called .gene_mappings which looks like: {1:{"name":Manes.12G059000.1, "scaffold_coordinate": 23213, "reference_chromosome": 12, "reference_coordinate":12432}, 2:...}
        # NOTE - we only care about starting coordinates, because the graphics are not detailed enough for specific start and end coordinates anyways...

        # It is important that we match all same genes on scaffold with all same genes on the reference genome... (not just pick one of them to connect)
        # self.genes => {'Manes.12G061900.1': [[1404095, 1405697]], 'Manes.12G059000.1': [[2613631, 2614138],[12,233]]...}

        gene_mappings = {}

        ref_genome_loc = get_gene_locations_from_reference(
            filtered_gff)  # puts all reference genome locations in one big dictionary... hopefully fast to search through :/
        match_num = 0
        for gene_name in self.genes:
            for scaff_instance in self.genes[
                gene_name]:  # Now scaff_instance is a specific coordinate of that gene on the scaffold
                scaff_instance_coordinate = scaff_instance[
                    0]  # Now we extracted only the starting coordinate of that instance (makes it easier later b/c we usually only need starting coordinate)
                # print "ref_genome: ",ref_genome_loc
                # print "Gene name:", gene_name
                for ref_instance in ref_genome_loc[
                    gene_name]:  # accesses all instances of that gene on the reference genome
                    match_num += 1
                    gene_mappings[match_num] = {"name": gene_name, "reference_chromosome": ref_instance['chromosome'],
                                                "reference_coordinate": ref_instance['start_coordinate'],
                                                "scaffold_coordinate": scaff_instance_coordinate}
        # {1:{"name":Manes.12G059000.1, "scaffold_coordinate": 23213, "reference_chromosome": 12, "reference_coordinate":12432}, 2:...

        self.gene_mappings = gene_mappings
        self.calculate_automatic_offset_and_mapping_stats()  # just calling it at the end of this function because now we have enough information
        return None

    def get_manual_attributes(self, manual_dicationary):
        # This  gives the scaffold object the manual x_offset and vertical pointer values from the manual specification dic
        attributes = manual_dicationary[self.id]
        try:
            self.x_offset = attributes['x_offset']
        except:
            None
        self.vertical_pointer = attributes['vertical_pointer']
        self.line_color = attributes['line_color']

    def __repr__(self):
        return "<Scaffold Object> with ID: " + self.id


def get_new_color(current_color_index):
    RGBA_tuple = colors[current_color_index]
    current_color_index += 1
    return RGBA_tuple, current_color_index

def plot_thin_marker(x, y, color):
    # This function draws very thin vertical lines by doing it with horizontal lines because vertical bars are too fat!
    # plt.plot(scaffold_coordinate + x_offset, vertical_pointer, '|', color=mutual_color)  # draw it in green
    # A hack to get the lines of the above command thinner...
    offsetter = -12
    for i in range(1, 60):
        offsetter += .4
        plt.plot([x, x + 1000], [y + offsetter, y + offsetter], color=color, linestyle='-', linewidth=0.1, alpha=1)
    return None

def get_gene_locations_from_reference(gff_pandas_object):
    # This function generates a ref_genome_loc dictionary which contains the coordinates for each gene in the reference genome
    # {"Manes.12G059000.1":[{'chromsome':chr12, 'start_coordinate': 12, 'end_coordinate': 421}, {'chromsome':chr11, 'start_coordinate': 1231, 'end_coordinate': 5331}], ...}
    # NOTE - dictionary in list in dictionary...
    ref_genome_loc = {}
    pattern = re.compile('ID.*.(?=.v6.1;)')  # for pattern matching in regular expressions when extracting gene_id

    for index, row in gff_pandas_object.iterrows():
        chr_name = row[0]
        start_coord = float(row[3])
        end_coord = float(row[4])
        gene_name = row[8]
        # orignally looks like this 'ID=Manes.01G000600.v6.1;Name=Manes.01G000600;ancestorIdentifier=cassava4.1_022534m.g.v4.1'
        # So I need to extract just the gene ID which comes right before ".v6.1;"...
        gene_name = re.match(pattern, gene_name).group()  # Extracts the gene name from the big string....
        gene_name = gene_name[3:]  # To delete the beginning 'ID=' part
        if gene_name in ref_genome_loc:  # The key is already in there so we need to append the additional coordinates
            ref_genome_loc[gene_name].append(
                {'chromsome': chr_name, 'start_coordinate': start_coord, 'end_coordinate': end_coord})
        else:  # The key doesn't exist yet, so we can initiallize the key and add the first coordinates
            ref_genome_loc[gene_name] = [
                {'chromosome': chr_name, 'start_coordinate': start_coord, 'end_coordinate': end_coord}]

    return ref_genome_loc

def pre_filter_phytozome_gff_file(gff_file):
    '''
    This function pre-processes our initially large reference genome file Because we only want gene information from the reference file, not cDNA, mRNA, etc.

    Then the ref_gene_information file is the Cassava Reference Genome v6.1 file from Phytozome (Mesculenta_305_v6.1.gene.gff3
    With the following columns:
    Start coordinate in reference genome column: 3
    End coordinate in reference genome column: 4
    Gene_ID (includes other junk) column: 8
    ** We only care about rows that have "gene" in column 2, otherwise it also includes mRNA, CDS, etc. information - we only need gene coordinates


    :param ref_genome_gff_file: downloaded from phytozome (tested with Cassava v6.1)
    :return: a pandas object that contains the reference gene information
    '''
    # filter the .gff3 file so it includes only rows with gene coordinates (Filter out mRNA, cDNA, and other rows objects)
    ref_gene_info = pandas.read_csv(gff_file, sep='\t', header=None)  # load in the reference gene information file
    ref_gene_info = ref_gene_info[ref_gene_info[
                                      2] == "gene"]  # Filter so we only keep rows of gene coordinate information (not cDNA, RNA, etc...)
    return ref_gene_info

def generate_scaffold_objects(list_of_scaffold_names, path_to_input_fasta, psl_file):
    # Generates a dictionary with multiple scaffold objects as specified.
    # {scaffold_1: Object, scaffold2: Object...}
    dict_of_scaffold_objects = {}
    list_of_found_scaffolds = []
    current_color_index = 0  # to help us cycle through the color map

    for seq_record in SeqIO.parse(path_to_input_fasta, "fasta"):
        if seq_record.id in list_of_scaffold_names:
            line_color = colors[current_color_index]  # needed to get new colors for the match lines
            current_color_index += 1  # needed to get new colors for the match lines
            index = list_of_scaffold_names.index(seq_record.id)  # get the index of the scaffold ID in the list
            name = list_of_scaffold_names[
                index]  # Delete that item from the scaffold list, so at the end we know what we didn't find!
            list_of_found_scaffolds.append(name)
            dict_of_scaffold_objects[seq_record.id] = ScaffoldBNG(seq_record, psl_file,
                                                                  line_color)  # Add the seq Record to the dictionary

    not_found = list(set(list_of_scaffold_names) - set(list_of_found_scaffolds))

    for item in not_found:  # Go through the original list and print out all the scaffold names that weren't found in the fasta file...
        print(item, " Was not found in the sequence file")
    return dict_of_scaffold_objects


#####################################
## Running the script
#####################################

cm = plt.get_cmap('Paired')
num_colors = len(list_of_scaffolds)  # Decides how many different colors we generate
## Generate enough RGBA Color tuples to have new line colors for each scaffold
colors = []
for i in range(num_colors):
    colors.append(cm(1.0 * i / num_colors))

scaff_dict = generate_scaffold_objects(list_of_scaffold_names=list_of_scaffolds,
                                       path_to_input_fasta=reference_fasta_file,
                                       psl_file=ref_psl_file)

## Get mappings to the reference genome annotations
for key in scaff_dict:
    scaff_dict[key].get_reference_mappings(reference_gff_file)


## Invert the scaffold objects if the user chose to do so:
for key in scaff_dict:
    if scaff_dict[key].id in scaffolds_to_invert:
        scaff_dict[key].invert_scaffold()

## PLOTTING
Plotter = OverallPlotter(scaff_dict)

try: # if there is a manual dictionary, then we will use it to make the plot
    Plotter.manual_scaffold_plotter(manual_plotting)
except:
    Plotter.auto_plot()

## This is just to create a standard sized output plot
plt.plot([-1000, -1000], [800, -800], linestyle='-', linewidth=0.0, color='r')
plt.plot([4500000, 4500000], [800, -800], linestyle='-', linewidth=0.0, color='r')


## Show the plot
#plt.show()

### Naming our plots uniquely so we don't keep overrwiting the old ones
cwd = os.getcwd()
plot_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
output_name = "../Output/scaphy_%s.pdf" % plot_time


# Saves the plot
plt.savefig(output_name, dpi=300, figsize=(400, 100))  # Switch between tme3 or 60444

print("You can find your plot at: %s/%s" % (cwd, output_name))

print("Have a nice day now and don't forget to take a break every now and then ;)")
print("Cheers! \n -your friendly SCEVT staff")
