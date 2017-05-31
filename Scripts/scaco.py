import pandas
import matplotlib.pyplot as plt
from Bio import SeqIO
import os
from datetime import datetime
import sys

'''
You'll need BioPython and matplotlib installed on your machine to make this script work
I recommend you get it via miniconda or anaconda.

To run this script, edit the part below this (Parameters) to what you like. Then execute this script by typing:

python scaco.py

The output will be saved in the same directory as <timestamp>_scaco_output.pdf

Written by:
Philipp v. Bieberstein
rphilipp@ethz.student.ch


'''

#####################################################################################################################
#########                                      || Parameters ||                                           ###########
#####################################################################################################################

# Fasta files were the scaffold sequences are included (Can be entire genome files, or just fasta files that contain the scaffold sequence)

reference_fasta_file = "../Data/Genome_sequences/60444_draft.fasta"

target_fasta_file = "../Data/Genome_sequences/TME3_draft.fasta"


# BLAT output files (.psl) that includes gene mappings onto the genomes which includes the specified scaffolds
# The paths should be relative to where this script is... or absolute paths
ref_psl_file = "../Data/Gene_BLAT_mappings/60444_BNG_plus_notscaff.psl"

target_psl_file = "../Data/Gene_BLAT_mappings/TME3_BNG_plus_notscaff.psl"

# Parameters for Plotting

N_region_min = 100 # threshold how big NNNN regions have to be in order to plot them

ref_scaff_x_offset = 0 # moving the reference scaffold to the left or  for nicer plots

target_scaff_x_offset = 300000 # moving the target scaffold to the left or right for nicer plots


# TME3 cmd2 scaffold list
#For example: 'Super-Scaffold_1951' or 'Super-Scaffold_730'

ref_scaffold = '004409F'

target_scaffold = "Super-Scaffold_1022"

# Whether you want to invert the reference or target scaffold:

invert_reference = False

invert_target = True

#####################################################################################################################
#########                                   || End of Parameters ||                                       ###########
#####################################################################################################################





class MatchObject:
    def __init__(self, ref_scaffold_object, target_scaffold_object):
        self.ref_object = ref_scaffold_object
        self.target_object = target_scaffold_object
        (unique_ref_genes, unique_target_genes, mutual_gene_connectors, total_ref_genes, total_target_genes,
         total_mutual_genes) = self.get_comparison_stats()

        # New awesome characteristics of this class
        self.unique_ref = unique_ref_genes
        self.unique_target = unique_target_genes
        self.mutual_genes = mutual_gene_connectors # looks like {1:{ref_coord:232, target_coord:5342, name:"masd222"}, 2:...}
        self.total_ref_genes = total_ref_genes
        self.total_target_genes = total_target_genes
        self.total_mutual_genes = total_mutual_genes


    def get_comparison_stats(self):
        '''
        This function takes two scaffold objects that we want to compare and generates 3 lists:
        ref_unique_genes
        target_unique_genes
        mutual_genes
        :param ref_scaffold:
        :param target_scaffold:
        :return:
        '''
        # Scaff_Object.genes => {'Manes.12G061900.1': [[1404095, 1405697]], 'Manes.12G059000.1': [[2613631, 2614138]]...}
        # Go through each gene in reference Genome,
        ref_scaffold = self.ref_object
        target_scaffold = self.target_object

        mutual_gene_connectors = {}  # Looks like this: {1:{ref_coord:232, target_coord:5342, name:"masd222"}, 2:...}
        unique_ref_genes = {} # Looks like this: {23:{"name": ref_gene, "ref_coord": ref_instance[0]}, 44:{"name": ref_gene, "ref_coord": ref_instance[0]}}
        unique_target_genes = {}
        match_number = 0
        total_ref_genes = 0
        total_target_genes = 0
        total_mutual_genes = 0
        for ref_gene in ref_scaffold.genes: # now ref_gene is the name of each gene on the ref_scaffold
            total_ref_genes += 1
            if ref_gene in target_scaffold.genes: # checks if this gene also exists on the target scaffold
                total_mutual_genes += 1
                for ref_instance in ref_scaffold.genes[ref_gene]: # If yes, get each coordinate instance of that gene from the ref_scaffold
                    for target_instance in target_scaffold.genes[ref_gene]: # for each instance on the ref_scaffold, get each instance of that gene on the target Scaffold

                        mutual_gene_connectors[match_number] = {"name": ref_gene,"ref_coordinate":ref_instance[0], "target_coordinate": target_instance[0]} # build up mutual gene dictionary
                        match_number += 1 # increase counter for mutual genes
            else: # if the gene does not exist on target scaffold, add it to a dictionary that only contains genes that are unique to the ref_scaffold
                for ref_instance in ref_scaffold.genes[ref_gene]: #
                    unique_ref_genes[match_number] = {"name": ref_gene, "coordinate": ref_instance[0]}  # build up unique to reference scaffold gene dictionary
                    match_number += 1

        for target_gene in target_scaffold.genes: # iterate through all target scaffold genes
            total_target_genes += 1
            if target_gene not in ref_scaffold.genes: # if the gene does not exist on the reference scaffold
                for target_instance in target_scaffold.genes[target_gene]: # add each instance of that gene to the unique to target dictionary
                    unique_target_genes[match_number] = {"name": target_gene, "coordinate": target_instance[0]}
                    match_number += 1

        # Calculate how many genes matched out of all genes


        return unique_ref_genes, unique_target_genes, mutual_gene_connectors, total_ref_genes, total_target_genes, total_mutual_genes




    def plot_scaffold(self, vertical_pointer=100,x_offset=0, scaffold_object=None):
        print("Plotting Scaffold: ", scaffold_object.id)
        plt.plot([0 + x_offset, scaffold_object.length + x_offset], [vertical_pointer, vertical_pointer], color='k',
                 linestyle='-', linewidth=4)
        text_location = vertical_pointer + 40
        # drawing the gaps
        for key in scaffold_object.gaps:
            gap = scaffold_object.gaps[key] # returns a list with the start of the gap and the end of that same gap
            start = gap[0]
            end = gap[1]
            if abs(end-start) > N_region_min: # checks if gap is big enough for us to care about drawing it
                offsetter = -12 # silly hack to make the gap drawing fine enough for small gaps... otherwise the giant red blots just overlap way too much
                for i in range(1,48):
                    offsetter += .5
                    plt.plot([start+x_offset, end+x_offset], [vertical_pointer + offsetter, vertical_pointer + offsetter], color='r', linestyle='-', linewidth=0.1, alpha=1)

        # Now add the labeling
        # the scaffold name
        plt.text(0 + x_offset, text_location, scaffold_object.id, fontsize=5)

    def plot_statistics(self, ref_vertical_pointer, target_vertical_pointer, ref_scaff_x_offset=0, target_scaff_x_offset=0):
        # draw statistics for reference scaffold
        text_location = ref_vertical_pointer - 10
        x_offset = ref_scaff_x_offset
        reference_stat = str(self.total_mutual_genes) + "/" + str(self.total_ref_genes) + " Genes Match"
        plt.text(self.ref_object.length + x_offset + 100000, text_location, reference_stat, fontsize=5)

        # draw statistics for target scaffold
        text_location = target_vertical_pointer - 10
        x_offset = target_scaff_x_offset
        target_stat = str(self.total_mutual_genes) + "/" + str(self.total_target_genes) + " Genes Match"
        plt.text(self.target_object.length + x_offset + 100000, text_location, target_stat, fontsize=5)

    def plot_unique_genes(self, ref_vertical_pointer, target_vertical_pointer=0, ref_scaff_x_offset=0, target_scaff_x_offset=0):
        # mutual_gene_connectors = {}  # Looks like this: {1:{ref_coord:232, target_coord:5342, name:"masd222"}, 2:...}
        # unique_ref = {} # Looks like this: {23:{"name": ref_gene, "ref_coord": ref_instance[0]}, 44:{"name": ref_gene, "ref_coord": ref_instance[0]}}
        print("Plotting unique genes on Scaffolds")
        # First, draw unique genes on reference scaffold
        x_offset = ref_scaff_x_offset
        vertical_pointer = ref_vertical_pointer
        for key in self.unique_ref:
            #name = self.unique_ref[key]["name"]
            coord = self.unique_ref[key]["coordinate"]
            plot_thin_marker(coord + x_offset, vertical_pointer, color='g')  # draw it in green

        # Draw unique genes on target scaffold
        x_offset = target_scaff_x_offset
        vertical_pointer = target_vertical_pointer
        for key in self.unique_target:
            #name = self.unique_target[key]["name"]
            coord = self.unique_target[key]["coordinate"]
            plot_thin_marker(coord + x_offset, vertical_pointer, color='g')  # draw it in green
        return None

    def plot_mutual_genes(self, ref_vertical_pointer, target_vertical_pointer, ref_scaff_x_offset=0, target_scaff_x_offset=0):
        # mutual_gene_connectors = {}  # Looks like this: {1:{ref_coordinate:232, target_coordinate:5342, name:"masd222"}, 2:...}
        print("Plotting Mutual Gene connections...")
        # get the coordinates for each match
        for key in self.mutual_genes:
            ref_coord = self.mutual_genes[key]["ref_coordinate"]
            target_coord = self.mutual_genes[key]["target_coordinate"]
            # Draw gene location of reference scaffold
            plot_thin_marker(ref_coord + ref_scaff_x_offset, ref_vertical_pointer, color='m')  # draw it in green

            # Draw gene location on target scaffold
            plot_thin_marker(target_coord + target_scaff_x_offset, target_vertical_pointer, color='m')  # draw it in green

            # Draw connection between scaffolds
            plt.plot([ref_coord + ref_scaff_x_offset, target_coord + target_scaff_x_offset],
                     [ref_vertical_pointer, target_vertical_pointer], linestyle='-', color='m',
                     linewidth=0.1, alpha=0.5)
        return None



    def plot_match_scaffolds(self):
        ref_vertical_pointer = 200
        target_vertical_pointer = -200
        # Plot Reference Scaffold
        self.plot_scaffold(vertical_pointer=ref_vertical_pointer, scaffold_object=self.ref_object, x_offset=ref_scaff_x_offset)

        # Plot Target Scaffold
        self.plot_scaffold(vertical_pointer=target_vertical_pointer, scaffold_object=self.target_object, x_offset=target_scaff_x_offset)

        # Plot Stats (# genes out of # are mutual)
        self.plot_statistics(ref_vertical_pointer, target_vertical_pointer, ref_scaff_x_offset, target_scaff_x_offset)

        # Plot the Unique genes
        self.plot_unique_genes(ref_vertical_pointer, target_vertical_pointer, ref_scaff_x_offset, target_scaff_x_offset)

        # Plot the mutual genes
        self.plot_mutual_genes(ref_vertical_pointer, target_vertical_pointer, ref_scaff_x_offset, target_scaff_x_offset)




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
        # self.genes => {'Manes.12G061900.1': [[1404095, 1405697]], 'Manes.12G059000.1': [[2613631, 2614138]]...}

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
                    almost_same = False
                    if gene_name in gene_dict:
                        for coordinates in gene_dict[gene_name]: # coordinates will be a list of start & end [start, end]
                            # gene exists already in the gene dictionary, so we have to check that the match is different enough for us to care.. (not splicing)
                            if abs(coordinates[0] - coordinate_range[0]) < 300 or abs(coordinates[1] - coordinate_range[1]) < 300:
                                almost_same = True # probably just difference in splicing, so we can ignore it
                            else:
                                None
                        if almost_same:
                            None
                        else: # coordinates are not just different splicing, so we add them to the dictionary
                            gene_dict[gene_name].append(coordinate_range)
                    else: # we have no coordinates for this gene yet, so we add it freshly
                        gene_dict[gene_name] = [coordinate_range] # that gene doesn't exist yet, so I add it freshly
        return gene_dict


    def invert_scaffold(self):
        '''
        This function aims to invert this scaffold including its NNN regions, and its gene locations
        :return:
        '''
        half = int(self.length/2)
        for key in self.gaps: # invert the gap coordinates
            self.gaps[key] = self.mirror_lists(half, self.gaps[key]) # give it the midpoint of the object and the list of coordinates it should invert on that object

        new_gene_coordinates = {}
        for key in self.genes: # {"gene1:[[22,532],[2343,64545]],"gene2:[[2,44]]...}
            for instance in self.genes[key]:
                if key in new_gene_coordinates:
                    new_gene_coordinates[key].append(self.mirror_lists(half, instance))
                else:
                    new_gene_coordinates[key] = [self.mirror_lists(half, instance)]
        self.genes = new_gene_coordinates
        print self.id + " was Inverted."
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


    def __repr__(self):
        return "<Scaffold Object> with ID: " + self.id


def generate_scaffold_object(scaffold_name, path_to_input_fasta, psl_file):
    # Generates one specific scaffold object by scanning through the fasta file, finding the scaffold ID, and creating a scaffold object from that
    for seq_record in SeqIO.parse(path_to_input_fasta, "fasta"):
        if seq_record.id == scaffold_name:
            print(seq_record.id, " Found")
            scaffold_object = ScaffoldBNG(seq_record, psl_file)
            print("Length: ", scaffold_object.length)
            return scaffold_object
    print(scaffold_name, " Could not be found in the fasta file -> ERROR")
    print('Try with a different Scaffold, or make sure you gave me the right FASTA file')
    raise SystemExit(0)

def plot_thin_marker(x, y, color):
    # This function draws very thin vertical lines by doing it with horizontal lines because vertical bars are too fat!
    # plt.plot(scaffold_coordinate + x_offset, vertical_pointer, '|', color=mutual_color)  # draw it in green
    # A hack to get the lines of the above command thinner...
    offsetter = -12
    for i in range(1, 60):
        offsetter += .4
        plt.plot([x, x+1000], [y+offsetter, y+offsetter], color=color, linestyle='-', linewidth=0.1, alpha=1)
    return None





#####################################
## Running the script
#####################################
print(sys.version)


ref_scaff = generate_scaffold_object(scaffold_name=ref_scaffold, path_to_input_fasta=reference_fasta_file,
                                     psl_file=ref_psl_file)

target_scaff = generate_scaffold_object(scaffold_name=target_scaffold, path_to_input_fasta=target_fasta_file,
                                        psl_file=target_psl_file)


# Invert the scaffold object if the user chose to do so:
if invert_reference:
    ref_scaff.invert_scaffold()
if invert_target:
    target_scaff.invert_scaffold()

# create a match instance
match = MatchObject(ref_scaff, target_scaff)

# PLOTTING

# This is just to create a standard sized output plot
plt.plot([-1000, -1000], [800, -800], linestyle='-', linewidth=0.0, color='r')
plt.plot([4500000, 4500000], [800, -800], linestyle='-', linewidth=0.0, color='r')

# Plots the scaffold and genes
match.plot_match_scaffolds()




cwd = os.getcwd()
plot_time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
output_name = "../Output/scaco_%s.pdf" % plot_time


# Saves the plot
plt.savefig(output_name, dpi=300, figsize=(400, 100))  # Switch between tme3 or 60444

# Show the plot
# plt.show()

print("You can find your plot at: %s/%s" % (cwd, output_name))

print "Have a nice day now and don't forget to take a break every now and then ;)"
print "Cheers! \n -your SCEVT staff"



