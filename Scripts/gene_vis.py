import pandas
import pickle
import re
import matplotlib.pyplot as plt
import random
from Bio import SeqIO
from numpy import sign, mean

random.seed(1)

def get_rand_color(index):
    colors = "bbrrbbrrcmbrcmbgrcm"
    return colors[index]


def create_temp_files(psl_gene_mapping_file, ref_gene_info_gff_file, cmd2_scaffolds_list):
    '''
    This function pre-processes our initially large files (because they contain information about all genes on the entire genome)
    We only want genes that mapped to our cmd2 specific scaffolds && we only want gene information from the reference file, not cDNA, mRNA, etc.
    :param psl_gene_mapping_file:
    :param ref_gene_info_gff_file:
    :param cmd2_scaffolds_list: List of cmd2 related scaffolds
    :return:
    '''
    # Filter the .psl file so it only includes CMD2 scaffold information
    psl_gene_mappings = pandas.read_csv(psl_gene_mapping_file, sep='\t', header=None) # load in the gene mapping .psl file
    psl_gene_mappings = psl_gene_mappings[psl_gene_mappings[13].isin(cmd2_scaffolds_list)] # Filter so we only keep rows of relevant cmd2 scaffolds
    # Filtering out genes that don't have at least 90% matching bps of the mRNA
    #print psl_gene_mappings[10]
    psl_gene_mappings = psl_gene_mappings[psl_gene_mappings[0] > 0.95*(psl_gene_mappings[10])]
    psl_gene_mappings.to_pickle("temp_filtered_gene_mapping.pickle") # save it as a temporary pickle python file!

    # and filter the .gff3 file so it includes only rows with gene coordinates and maybe only rows which describe genes that are included in the .psl file after filtering
    ref_gene_info = pandas.read_csv(ref_gene_info_gff_file, sep='\t', header=None) # load in the reference gene information file
    ref_gene_info = ref_gene_info[ref_gene_info[2] == "gene"] # Filter so we only keep rows of gene coordinate information (not cDNA, RNA, etc...)
    ref_gene_info.to_pickle("temp_reference_info.pickle") # save it as a temporary pickle python file!
    return None


def N_region_finder(seq_record):
    '''
    Find the regions of NNNNs in the sequence and records the indices of starting and ending points in two lists.
    Return a tuple of two lists. If the sequence starts of known BasePairs, then the 'N-region_end' list will have index 0.
    :param seq_record:
    :return: returns a tuple of two lists. (N_region_start, N_region_end) which each contain the coordinates of either starting or
    ending 'N' regions (gaps)
    '''
    # Initialize N_region_start & N_region_end lists
    # Initialize the vertical position of the pointer for plotting (move ~50 px upwards)
    N_region = False
    N_region_start = []
    N_region_end = []

    sequence = str(seq_record.seq).upper()
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

    return (N_region_start, N_region_end)

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


def invert_scaffold(scaffold_id, seq_len,N_region_start, N_region_end, gene_ids, gene_coord):
    # This function takes all relevant coordinate lists for a scaffold and then
    # Inverts all relevant ones so that when we plot it, the scaffold and all features are simply inverted...
    # gene_ids = {"scaffold1":[Manes.12G070700, Manes.12G070200]}
    # gene_coord {"scaffold1_Manes.12G070700": [[300, 500], [649, 859]], "scaffold1_Manes.12G079900": [[0, 220]]}
    # ref_genome_loc {"Manes.12G070700": [[ChrXII, 30, 500], [ChrXII, 900, 1040]]}
    half = seq_len/2
    N_region_start = mirror_lists(half,N_region_start)
    N_region_end = mirror_lists(half,N_region_end)
    # This for loop will invert each of the lists of a gene (because there could me multiple copies of 1 gene)
    for gene in gene_ids[scaffold_id]:
        scaffold_gene_id = scaffold_id + "_" + gene
        for i, instance in enumerate(gene_coord[scaffold_gene_id]):
            gene_coord[scaffold_gene_id][i] = mirror_lists(half, instance)

    return (N_region_start, N_region_end, gene_coord)


def get_gene_coordinates_from_pickle(temp_gene_mapping_filename, temp_ref_gene_info_filename):
    '''
    This function gets the coordinate information from all the mapped genes on our scaffolds. The gene_mapping_filename is a .psl file which looks like this:
    matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName   qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts
    The Gene name is in column: 9
    The scaffold that the gene was mapped to is in column: 13
    The coordinate where the match STARTS within the scaffold/contig is in column: 15
    The coordinate where the match ENDS within the scaffold/contig is in column: 16

    Then the ref_gene_information file is the Cassava Reference Genome v6.1 file from Phytozome (Mesculenta_305_v6.1.gene.gff3
    With the following columns:
    Start coordinate in reference genome column: 3
    End coordinate in reference genome column: 4
    Gene_ID (includes other junk) column: 8
    ** We only care about rows that have "gene" in column 2, otherwise it also includes mRNA, CDS, etc. information - we only need gene coordinates

    :param filename: name of the csv file (*Actually its a .psl file) and another file which was downloaded from Phytozome gene.gff3 file (deleted the first 3 rows... header information)
    :return:
    Returns Panda dataforms
    returns three dictionaries -
    gene_id = one that contains the scaffold_ID with all gene_ids that mapped to it {"scaffold1":[Manes.12G070700.1,Manes.12G070200.1,snp19]}
    gene_cord = One that contains the physical coordinates of each gene specific to scaffold {"scaffold1_Manes.12G070700.1":[[300,500], [649,859]], "scaffold1_Manes.12G079900.1":[[0,220]] }
    (Can have multiple copies of the same gene on the same scaffold)
    ref_genome_loc = A dictionary that contains the gene_ID locations on the reference genome {"Manes.12G070700.1":[[ChrXII, 30,500], [ChrXII, 900, 1040]]
    (Can occur multiple times on the reference genome && also on different chromosomes!
    '''

    gene_ids = {}
    gene_coord = {}
    ref_genome_loc = {}
    gene_list = []

    ####### read in the pickled files from pre-processing
    gene_mapping = pandas.read_pickle(temp_gene_mapping_filename)
    ref_info = pandas.read_pickle(temp_ref_gene_info_filename)

    pattern = re.compile('.*.(?=.\d)') # Pattern matching to gene underlying gene name instead of CDS names...
    # Filling the gene_id dictionary
    for index, row in gene_mapping.iterrows():
        scaff_name = row[13]
        gene_name = row[9]
        gene_name = re.match(pattern, gene_name).group() # This should erase the CDS specifications at the end of the gene names ".3"

        start_coord = float(row[15])
        end_coord = float(row[16])
        scaffold_gene_id = scaff_name + "_" + gene_name
        if scaff_name not in gene_ids:
            gene_list.append(gene_name)
            gene_ids[scaff_name] = [gene_name]
            gene_coord[scaffold_gene_id] = [[start_coord, end_coord ]]
            #print "Key added: ", scaff_name+"_"+gene_name
        else:
            #print "Adding this key: ", scaff_name+"_"+gene_name
            if gene_name not in gene_list: # Don't add duplicate genes
                gene_list.append(gene_name)
            gene_ids[scaff_name].append(gene_name)
            almost_same = False  # boolean to signal if the new gene coordinate is actually on a different location on scaffold, or if its just different splicing!
            if scaffold_gene_id in gene_coord: # Need to check this, because there could be multiple copies of the same gene on the same scaffold!!
                for coordinates in gene_coord[scaffold_gene_id]: # Do some naive checks to try not to add different splicing matches... because the coordinates will be almost the same and create unnecessary lines in plot
                    if abs(coordinates[0] - start_coord) < 300 or abs(coordinates[1] - end_coord) < 300:
                        almost_same = True
                    else:
                        None
                if almost_same:
                    None
                else:
                    gene_coord[scaffold_gene_id].append([start_coord, end_coord])
            else:
                gene_coord[scaffold_gene_id] = [[start_coord, end_coord]]

    pattern = re.compile('ID.*.(?=.v6.1)')  # for pattern matching in regular expressions when extracting gene_id
    for index, row in ref_info.iterrows():
        chr_name = row[0]
        start_coord = float(row[3])
        end_coord = float(row[4])
        gene_name = row[8]
        # orignally looks like this 'ID=Manes.01G000600.v6.1;Name=Manes.01G000600;ancestorIdentifier=cassava4.1_022534m.g.v4.1'
        # So I need to extract just the gene ID which comes right before ".v6.1;"...
        gene_name = re.match(pattern, gene_name).group() # Extracts the gene name from the big string....
        gene_name = gene_name[3:] # To delete the beginning 'ID=' part
        if gene_name in ref_genome_loc: # The key is already in there so we need to append the additional coordinates
            ref_genome_loc[gene_name].append([chr_name,start_coord, end_coord])
        else: # The key doesn't exist yet, so we can initiallize the key and add the first coordinates
            ref_genome_loc[gene_name] = [[chr_name,start_coord, end_coord]]

    return (gene_ids, gene_coord, ref_genome_loc, gene_list)

def plot_graphic(region_fasta_file, gene_ids, gene_coord, ref_genome_loc, gene_list, output_name):
    # Plotting the pseudo molecule
    plt.plot([0, pseudo_length], [0, 0], color='c', linestyle='-', linewidth=4)
    vertical_pointer = 1  # Just because if its 0, then the numpy sign is '0' and it breaks my technique for spreading scaffolds around the y=0 axis
    color_index = 0 # for the color index choosing
    pseudo_chromosome_loc = 0 # this is aritifically spreading out the x coordinates from genes that match to other chromosomes besides chr12

    for seq_record in SeqIO.parse(region_fasta_file, "fasta"):
        print color_index
        color = get_rand_color(color_index) # For the color choosing
        color_archive = color # Just so we can switch inbetween but always revert back to the main color for this scaffold
        seq_name = seq_record.id
        print seq_name
        (N_region_start, N_region_end) = N_region_finder(seq_record)
        print N_region_start
        print N_region_end
        # Get sequence length and store in seq_len
        seq_len = len(seq_record.seq)
        print seq_len

        vertical_pointer = (abs(vertical_pointer) + 100) * (sign(vertical_pointer) * -1)  # This is so we get alternating + and - coordinates to draw around the y=0 exis

        ####### Invert specific scaffolds
        try:
            if seq_name in invert_these_scaffolds:
                (N_region_start, N_region_end, gene_coord) = invert_scaffold(seq_name, seq_len,N_region_start, N_region_end, gene_ids, gene_coord)
        except:
            pass
        print "After Inversion!"
        print N_region_start
        print N_region_end
        manual_horizontal_moving = 0
        # Checking if we manually want to set the vertical pointer for each scaffold, or if we do it manually

        try:
            if seq_name in scaffold_y_ordering:
                vertical_pointer = scaffold_y_ordering[seq_name]
            else:
                continue
        except:
            pass
        # Same as above except for manual right and left moving now
        try:
            if seq_name in scaffold_x_ordering:
                manual_horizontal_moving = scaffold_x_ordering[seq_name]
            else:
                None
        except:
            pass



        # Draw entire chromosome in black
        # x_offset will make the scaffold plots more readable by offsetting the scaffolds to the right proportionaly to where the SNPs match


        # gene_ids = {"scaffold1":[Manes.12G070700, Manes.12G070200]}
        # gene_coord {"scaffold1_Manes.12G070700": [[300, 500], [649, 859]], "scaffold1_Manes.12G079900": [[0, 220]]}
        # ref_genome_loc {"Manes.12G070700": [[ChrXII, 30, 500], [ChrXII, 900, 1040]]}
        # gene_list = ['Manes.S103300', 'Manes.04G035500', 'Manes.04G106800', 'Manes.04G106900', 'Manes.S094100']
        # calculating offset
        num_matched_genes = str(len(gene_ids[seq_name]))
        avg_coordinates = 0
        temp_for_avg = []
        num_on_chr12 = 0
        num_elsewhere = 0
        duplicated_genes = 0
        for gene in gene_ids[seq_name]: # now 'gene' is a gene name that sits on this scaffold
            for item in ref_genome_loc[gene]: # now 'item' is a list of locations where this gene is in the ref genome # actually these genes exist only once on the ref genome!
                if len(ref_genome_loc[gene]) >1:
                    print "THIS GENE EXISTS FOR THAN ONCE!!!!!"+gene
                    print item
                    #  actually these genes exist only once on the ref genome! (Tested for current cmd2 locus on 60444)

                if item[0] == "Chromosome12":  # check if its on chromosome 12 on reference genome... otherwise print out a message
                    temp_for_avg.append(item[1]) # adding beginning coordinate of all matched genes to this list
                    avg_coordinates = mean(temp_for_avg)
                    num_on_chr12 += 1 # count how many genes are on chromosome 12
                else:
                    num_elsewhere += 1

        x_offset = avg_coordinates
        x_offset += manual_horizontal_moving # Allowing for manual moving of the scaffolds to the left or right

        plt.plot([0+x_offset, seq_len+x_offset], [vertical_pointer, vertical_pointer], color='k', linestyle='-', linewidth=4)

        # Label the scaffold name above the line
        text_location = vertical_pointer + 40
        plt.text(0+x_offset, text_location, seq_name, fontsize=5)

        # Plot how many genes matched on the specific scaffold
        # So we can write the number of genes that matched for each scaffold
        # for information to the right of scaffold
        chr12_gene_ratio_text= str(num_on_chr12)+"/"+str(num_elsewhere+num_on_chr12) + " genes on Chr12"
        plt.text(seq_len+x_offset+100000, text_location - 45, chr12_gene_ratio_text, fontsize=5)


        # Now draw N regions in red
        # base case for normal behavior of starting with known bases

        for index in range(len(N_region_start)):
            start = N_region_start.pop(0)
            end = N_region_end.pop(0)
            if abs(end - start) > N_region_min:
                offsetter = -12

                for i in range(1,48):
                    offsetter += .5
                    plt.plot([start+x_offset, end+x_offset], [vertical_pointer + offsetter, vertical_pointer + offsetter], color='r', linestyle='-', linewidth=0.1, alpha=1)



        # gene_ids = {"scaffold1":[Manes.12G070700, Manes.12G070200]}
        # gene_coord {"scaffold1_Manes.12G070700": [[300, 500], [649, 859]], "scaffold1_Manes.12G079900": [[0, 220]]}
        # ref_genome_loc {"Manes.12G070700": [[ChrXII, 30, 500], [ChrXII, 900, 1040]]}
        # gene_list = ['Manes.S103300', 'Manes.04G035500', 'Manes.04G106800', 'Manes.04G106900', 'Manes.S094100']
        if seq_name in gene_ids:
            color_index += 1  # For the color of SNP location mapping lines
            marker_text_location = vertical_pointer + 40
            for gene in gene_ids[seq_name]: # now 'gene' is a gene name that is located on this scaffold
                scaffold_gene_id = seq_name + "_" + gene # generates the scaffold_name+gene_name combination to find the specific coordinates
                phys_location = gene_coord[scaffold_gene_id]  # take the coordinate where the gene starts... length probably doesn't matter
                ref_location = []
                phys_coordinate = []
                alternative_references = []
                for gene_instance in phys_location:
                    phys_coordinate.append(gene_instance[0]) # Adds starting coordinate of every instance of that gene on that scaffold
                for item in ref_genome_loc[gene]:
                    if item[0] == "Chromosome12": # check if its on chromosome 12 on reference genome... otherwise print out a message
                        ref_location.append(item[1]) # Add all instances of where the gene is on the reference genome
                    else:
                        print "This gene does not sit on Chromosome 12 in Ref Genome: ", item
                        if item[0] == 'Scaffold03391': # This draws a black line where the S111500 gene is!
                            #plt.plot([phys_coordinate[0] +x_offset, phys_coordinate[0] + x_offset],
                            #         [vertical_pointer + offsetter, vertical_pointer + offsetter+60], color='k',
                            #         linestyle='-', linewidth=0.1, alpha=1)
                            None
                        alternative_references.append(item[1])
                #print ref_location
                #print phys_location
                # PLOTTING THE GENE LOCATION LINES ONTO THE PSEUDO-MOLECULE
                alpha_value = .4

                if len(phys_coordinate) > 1:
                    print "THIS GENE EXISTS MULTIPLE TIMES ON THE SCAFFOLD!!!!"
                    print gene
                    print phys_coordinate
                    duplicated_genes += 1
                    color = 'k' # Make the gene black if it exists more than once
                else:
                    color = color_archive # Revert to original color for this scaffold if the gene only exists once
                if len(ref_location) > 1:
                    print "this gene exists multiple times on the reference genome!!!!!!", ref_location
                for phys_instance in phys_coordinate:  # connect all lines from all instances of that gene in the scaffold
                    for instance in ref_location: # connect lines to all instances of that gene on the reference genome
                        #print "Drawing new line!!!"
                        print "Drawing a new line from:", phys_instance, " to: ", instance, "For gene: ", gene
                        plt.plot(phys_instance + x_offset, vertical_pointer, '|', color='y')
                        plt.plot([phys_instance+x_offset, instance], [vertical_pointer, 0], linestyle='-', color=color, linewidth=0.1, alpha=alpha_value)

                    # This marks in green on the scaffold when there is a gene that doesn't map to chr 12
                    for instance in alternative_references:
                        # Plots the molecules that are not on Chr2 to the top of the graph
                        plt.plot(phys_instance + x_offset, vertical_pointer, '|', color='g')
                        #plt.plot([phys_instance+x_offset, pseudo_chromosome_loc], [vertical_pointer, 1500], linestyle='-', color='g', linewidth=0.1, alpha=0.3)
                        #pseudo_chromosome_loc += 10000 # to artificially spread out green lines
                marker_text_location += 2  # so they don't overlay completely
            # Add text for how many genes are duplicated on that scaffold:
            #plt.text(seq_len + x_offset + 100000, text_location - 75, str(duplicated_genes) +" genes duplicated", fontsize=5)

    plt.plot([-1000,-1000],[800,-800], linestyle='-', linewidth=0.0, color='r')
    plt.plot([2500000,2500000],[800,-800], linestyle='-', linewidth=0.0, color='r')
    plt.savefig(output_name, dpi=300, figsize=(400, 100))  # Switch between tme3 or 60444

    plt.show()


############################
## MAIN SCRIPT
############################


pre_processing = True
psl_filename = 'TME3_BNG_plus_notscaff.psl'
#psl_filename = '60444_BNG_plus_notscaff.psl'

gff3_filename = 'Mesculenta_305_v6.1.gene.gff3'

# TME3 cmd2 scaffold list
cmd2_scaffolds_list_tme3 = ['Super-Scaffold_44', 'Super-Scaffold_1951', 'Super-Scaffold_730', 'Super-Scaffold_1022', '001839F', '003097F', '002893F', '006625F', '007880F']

cmd2_scaffolds_list_60444 = ['Super-Scaffold_111','Super-Scaffold_29', 'Super-Scaffold_885','006018F','Super-Scaffold_526','001529F','004409F','001787F','004969F','008344F','Super-Scaffold_572']


#pseudo_length = 31578400 # estimated chromosome 12 length from phytozome
pseudo_length = 11000000 # only part of it so we can see everything better


if pre_processing == True:
    # Filter the .psl file so it only includes CMD2 scaffold information
    create_temp_files(psl_filename, gff3_filename, cmd2_scaffolds_list_tme3)
    # and the .gff3 file so it includes only rows with gene coordinates and maybe only rows which describe genes that are included in the .psl file after filtering

(gene_ids, gene_coord, ref_genome_loc, gene_list) = get_gene_coordinates_from_pickle("temp_filtered_gene_mapping.pickle", "temp_reference_info.pickle")

# Manually decide where on the y axis the specific scaffolds should lie: { "scaffold_1":100, "scaffold_234":-100}

#scaffold_y_ordering = {"Super-Scaffold_111":200, "Super-Scaffold_526":300, "Super-Scaffold_29":-200,"Super-Scaffold_885":-300, "001787F":400}
# The scaffolds will the moved to the right or left by this amount

scaffold_y_ordering = {"Super-Scaffold_44":200, "Super-Scaffold_1022":300, "Super-Scaffold_1951":-200,"Super-Scaffold_730":-300}#, "001787F":400}

#scaffold_x_ordering = {"Super-Scaffold_111":-1000000, "Super-Scaffold_526":0, "Super-Scaffold_29":-1000000,"Super-Scaffold_885":0}

scaffold_x_ordering = {"Super-Scaffold_44":-1500000, "Super-Scaffold_1022":-300000, "Super-Scaffold_1951":-900000,"Super-Scaffold_730":-300000}

## Manually invert specific scaffolds:
invert_these_scaffolds = ["Super-Scaffold_111", "Super-Scaffold_885", "Super-Scaffold_526", "001787F", "Super-Scaffold_44", "001839F"]

# Specify what the minimum NNN region size needs to be in order for the script to mark it as a gap:
N_region_min = 100

#print gene_ids
#print gene_coord
#print ref_genome_loc
#print gene_list
output_name = 'tme3_graphic.pdf'
region_fasta_file = "tme3_07april_cmd2_region.fasta"


plot_graphic(region_fasta_file, gene_ids, gene_coord, ref_genome_loc, gene_list, output_name)


# DONE: Fix bug that prevents NNNN regions from showing correctly when flipping the scaffolds!!
# DONE: plot how many genes are there fully 0-XX compared to how many genes are just partially there...
# DONE: Invert scaffolds manually if we want!
# DONE: mark how many genes from chr12 vs. others
# DONE: Filter gaps so that only big ones get shown... not the tiny ones
# DONE: Check the double gene copies on each scaffolds