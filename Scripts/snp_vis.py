'''
This script should generate a matplot lib plot with pixels in black when we have a real sequences and gaps in the dots
when we have 'N's in the sequence... This way we have an image with proportional gaps etc.
'''

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import matplotlib.pyplot as plt
import csv
from numpy import sign
import random


def get_rand_color(index):
    colors = "bbggrrkkbbggrrcmykbgrcmykbgrcmyk"
    return colors[index]

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
                N_region_end.append(0)
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



def get_SNP_coordinates_from_csv(filename):
    '''
    This function gets the coordinate information form the mapped SNPs. The csv file should look like this:
    SNP_ID, ChromosomeID, cM_distance, Scaffold_ID, Physical_Distance
    :param filename: name of the csv file
    :return: returns three dictionaries, one that contains the SNP ids {"scaffold1":[snp1,snp5,snp19]}
    One that contains the physical distances of each SNP specific to scaffold {"scaffold1_snp2":2412}
    And one that contains the cM distances as above {"scaffold1_snp2":31.2}
    '''
    snp_ids = {}  # Dictionary to hold for each scaffold all the SNP IDs that matched to it! {"scaffold1":[snp1,snp5,snp19]}
    snp_phys_dic = {}  # Dictionary which holds for each SNP_id+scaffold match, the physical distance on the scaffold {"scaffold1_snp2":2412}
    snp_cM_dic = {}  # Same as above but for cM distances {"scaffold1_snp2":31.2}
    ####### read in csv file
    with open(filename, 'rb') as csv_file:  # Switch between tme3 and 60444
        # with open('60444_SNP_hits.csv', 'rb') as csv_file:
        snp_file = csv.reader(csv_file, delimiter=',')

        for row in snp_file:
            print row
            if row[3] not in snp_ids:  # if its not in there, initiate it in the dictionaries
                snp_ids[row[3]] = [row[0]]
                snp_phys_dic[row[3] + "_" + row[0]] = float(row[4])
                snp_cM_dic[row[3] + "_" + row[0]] = float(row[2])

            else:  # Otherwise just append the values to the exisitng dictionary values... in list to make sure the order stays the same!
                snp_phys_dic[row[3] + "_" + row[0]] = float(row[4])
                snp_cM_dic[row[3] + "_" + row[0]] = float(row[2])
                snp_ids[row[3]].append(row[0])
        print snp_ids
        print snp_cM_dic
    return (snp_ids, snp_phys_dic, snp_cM_dic)


def get_gene_coordinates_from_csv(gene_mapping_filename, ref_gene_info_filename):
    '''
    This function gets the coordinate information from all the mapped genes on our scaffolds. The gene_mapping_filename is a .psl file which looks like this:
    matches	misMatches	repMatches	nCount	qNumInsert	qBaseInsert	tNumInsert	tBaseInsert	strand	qName   qSize	qStart	qEnd	tName	tSize	tStart	tEnd	blockCount	blockSizes	qStarts	tStarts
    The Gene name is in column: 9
    The scaffold that the gene was mapped to is in column: 15
    The coordinate where the match STARTS within the scaffold/contig is in column: 17
    The coordinate where the match ENDS within the scaffold/contig is in column: 18

    Then the ref_gene_information file is the Cassava Reference Genome v6.1 file from Phytozome (Mesculenta_305_v6.1.gene.gff3
    With the following columns:
    Start coordinate in reference genome column: 3
    End coordinate in reference genome column: 4
    Gene_ID (includes other junk) column: 8
    ** We only care about rows that have "gene" in column 2, otherwise it also includes mRNA, CDS, etc. information - we only need gene coordinates

    :param filename: name of the csv file (*Actually its a .psl file) and another file which was downloaded from Phytozome gene.gff3 file (deleted the first 3 rows... header information)
    :return:
    Returns Panda dataforms
    returns three dictionaries - one that contains the scaffold_ID with all gene_ids that mapped to it {"scaffold1":[Manes.12G070700.1,Manes.12G070200.1,snp19]}
    One that contains the physical distances of each SNP specific to scaffold {"scaffold1_snp2":2412}
    And one that contains the cM distances as above {"scaffold1_snp2":31.2}
    '''
    snp_ids = {}  # Dictionary to hold for each scaffold all the SNP IDs that matched to it! {"scaffold1":[snp1,snp5,snp19]}
    snp_phys_dic = {}  # Dictionary which holds for each SNP_id+scaffold match, the physical distance on the scaffold {"scaffold1_snp2":2412}
    snp_cM_dic = {}  # Same as above but for cM distances {"scaffold1_snp2":31.2}
    ####### read in csv file
    with open(filename, 'rb') as csv_file:  # Switch between tme3 and 60444
        # with open('60444_SNP_hits.csv', 'rb') as csv_file:
        snp_file = csv.reader(csv_file, delimiter=',')

        for row in snp_file:
            print row
            if row[3] not in snp_ids:  # if its not in there, initiate it in the dictionaries
                snp_ids[row[3]] = [row[0]]
                snp_phys_dic[row[3] + "_" + row[0]] = float(row[4])
                snp_cM_dic[row[3] + "_" + row[0]] = float(row[2])

            else:  # Otherwise just append the values to the exisitng dictionary values... in list to make sure the order stays the same!
                snp_phys_dic[row[3] + "_" + row[0]] = float(row[4])
                snp_cM_dic[row[3] + "_" + row[0]] = float(row[2])
                snp_ids[row[3]].append(row[0])
        print snp_ids
        print snp_cM_dic
    return (snp_ids, snp_phys_dic, snp_cM_dic)

def plot_graphic(region_fasta_file, snp_ids, snp_phys_dic, snp_cM_dic, output_name):
    # Plotting the pseudo molecule
    plt.plot([offset_pseudo_molecule, pseudo_length + offset_pseudo_molecule], [0, 0], color='c', linestyle='-',
             linewidth=4)
    # Need to offeset to the right because the molecule doesn't start at 0cM!!!
    vertical_pointer = 1  # Just because if its 0, then the numpy sign is '0' and it breaks my technique for spreading scaffolds around the y=0 axis
    color_index = 0 # for the color index choosing

    for seq_record in SeqIO.parse(region_fasta_file, "fasta"):  # Switch between tme3 and 60444
        color_index += 1 # For the color of SNP location mapping lines
        print color_index
        color = get_rand_color(color_index) # For the color choosing
        # for seq_record in SeqIO.parse("60444_07april_cmd2_region.fasta", "fasta"):
        seq_name = seq_record.id
        print seq_name
        (N_region_start, N_region_end) = N_region_finder(seq_record)
        print N_region_start
        print N_region_end
        # Get sequence length and store in seq_len
        seq_len = len(seq_record.seq)
        print seq_len
        vertical_pointer = (abs(vertical_pointer) + 100) * (sign(vertical_pointer) * -1)  # This is so we get alternating + and - coordinates to draw around the y=0 exis
        try:
            if seq_name in scaffold_y_ordering:
                vertical_pointer = scaffold_y_ordering[seq_name]
            else:
                continue
        except:
            pass



        # Draw entire chromosome in black
        # x_offset will make the scaffold plots more readable by offsetting the scaffolds to the right proportionaly to where the SNPs match
        # snp_ids {"scaffold1": [snp1, snp5, snp19]}
        # snp_phys_dic {"scaffold1_snp2": 2412}
        # snp_cM_dic {"scaffold1_snp2": 31.2}
        min_cM_location = 9999

        for snp in snp_ids[seq_name]:
            cM_location = snp_cM_dic[seq_name + "_" + snp]
            offset_phys_location = snp_phys_dic[seq_name + "_" + snp]
            if cM_location < min_cM_location:
                min_cM_location = cM_location
                min_phys_location = offset_phys_location
        x_offset = min_cM_location*bp_to_cM_ratio # This is the offset that all information should move to the right for better readability
        x_offset = x_offset - min_phys_location


        plt.plot([0+x_offset, seq_len+x_offset], [vertical_pointer, vertical_pointer], color='k', linestyle='-', linewidth=4)
        seq_len = len(seq_record.seq)

        # Label the scaffold name above the line
        text_location = vertical_pointer + 40
        plt.text(0+x_offset, text_location, seq_name, fontsize=5)

        # Now draw N regions in red
        # base case for normal behavior of starting with known bases
        if N_region_end[0] == 0:
            N_region_end.pop(0)
        for index in range(len(N_region_start)):
            start = N_region_start.pop(0)
            end = N_region_end.pop(0)
            plt.plot([start+x_offset, end+x_offset], [vertical_pointer, vertical_pointer], color='r', linestyle='-', linewidth=4,
                     alpha=0.3)

        # Now draw in the SNP locations && Identities!!!
        # snp_ids = {"scaffold1":[snp1,snp5,snp19]}
        # snp_phys_dic = {"scaffold1_snp2":2412}
        # snp_cM_dic = same as above but cM
        if seq_name in snp_ids:
            marker_text_location = vertical_pointer + 1
            cM_text_location = vertical_pointer - 1
            for snp in snp_ids[seq_name]:
                phys_location = snp_phys_dic[
                    seq_name + "_" + snp]  # thats the specific dictionary key (seq_name+"_"+snp)
                cM_location = snp_cM_dic[seq_name + "_" + snp]  # thats the specific dictionary key (seq_name+"_"+snp)
                plt.plot(phys_location+x_offset, vertical_pointer, '|', color='y')
                plt.text(phys_location+x_offset, marker_text_location, snp, fontsize=0.5, color='r')
                plt.text(phys_location+x_offset, cM_text_location, cM_location, fontsize=0.5, color='r')
                # PLOTTING THE SNP LOCATION LINES ONTO THE PSEUDO-MOLECULE

                plt.plot([phys_location+x_offset, (cM_location * bp_to_cM_ratio)], [vertical_pointer, 0], linestyle='-',
                         color=color, linewidth=0.1)

                marker_text_location += 2  # so they don't overlow completely
                cM_text_location -= 2
    plt.plot([-1000000,-1000000],[1000,-1000], linestyle='-', linewidth=0.0, color='r')
    plt.plot([2500000,2500000],[1000,-1000], linestyle='-', linewidth=0.0, color='r')

    plt.savefig(output_name, dpi=300, figsize=(400, 100))  # Switch between tme3 or 60444

    plt.show()






################################
# PARAMETERS
################################

################################
# PSEUDO MOLECULE LENGTH

pseudo_length = 1106536.68 # Pseudo length for 60444 cmd2 region
#pseudo_length = 1154528.91 # Pseudo length for tme3 cmd2 region

################################
# bp TO cM RATION

bp_to_cM_ratio = 36268 # Ratio for 60444
#bp_to_cM_ratio = 37841 # Ratio for tme3

################################
# OFFESET FOR THE PSEUDO MOLECULE (This is the distance of the start from the CMD2 locus to the beginning of chr 12... because CMD2 locaus doesn't start at 0 cM)

offset_pseudo_molecule = 797896 # For 60444
#offset_pseudo_molecule = 832502 # For tme3

################################
# SNP FILENAME
snp_filename = '60444_SNP_hits.csv'
#snp_filename = 'tme3_SNP_hits.csv'

################################
# PDF OUTPUT FILENAME
#output_name = 'tme3_cmd2_region_snps.pdf'
output_name = '60444_cmd2_region_snps.pdf'

################################
# CMD2 REGION FASTA FILE
#region_fasta_file = "tme3_07april_cmd2_region.fasta"
region_fasta_file = "60444_07april_cmd2_region.fasta"


################################
############## MAIN ############
################################

# Manually decide where on the y axis the specific scaffolds should lie: { "scaffold_1":100, "scaffold_234":-100}
#scaffold_y_ordering = {"Super-Scaffold_111":400,"Super-Scaffold_29":-400}#, "Super-Scaffold_526":200, ,"Super-Scaffold_885":-200}
scaffold_y_ordering = {"Super-Scaffold_526":400, "Super-Scaffold_885":-400}

# Getting the SNP information
(snp_ids, snp_phys_dic, snp_cM_dic) = get_SNP_coordinates_from_csv(snp_filename)

# Plotting and saving the graph
plot_graphic(region_fasta_file, snp_ids, snp_phys_dic, snp_cM_dic, output_name)




