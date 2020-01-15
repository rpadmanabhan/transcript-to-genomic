import os
import sys

'''
Convert Transcript Coordinates to Genomic Coordinates
'''

def parse_cigar(cigar):
    '''
    '''

def translate_coordinates(genomic_start_pos, cigar):
    '''
    :param int genomic_start_pos: Start position of the transcript on the Genome
    :param str cigar: Cigar String representing the alignment of the transcript to the genome
    '''
    cigar_pattern = re.compile("([0-9]+)([A-Z])")
    
    for num_bases, cigar_char in re.findall(cigar_pattern, cigar):
        if cigar_
        
        if cigar_char == "D" or cigar_char == "N":
            genome_pos += int(num_bases)
        elif cigar_char == "I":
            transcript_pos += int(num_bases)


def run():
    '''
    '''

    transcript_coordinates = {}

    with open(infile1, "r") as IN:
        for line in IN:
            transcript, chrom, genomic_start_pos, cigar = line.strip("\n").split("\t")
            interval_tree_obj = translate_coordinates(genomic_start_pos, cigar)
            assert transcript not in transcript_coordinates, "Duplicate Transcript encountered"
            transcript_coordinates[transcript]  = [chrom, interval_tree_obj]

    with open(infile2, "r") as IN:
        for line in IN:
            transcript, transcript_pos = line.strip("\n").split("\t")
    
    

if __name__ == '__main__':
    run()
