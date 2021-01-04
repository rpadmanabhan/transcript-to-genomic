# stdlib
import datetime
import os
import re
import sys

# 3rd Party
import intervaltree

'''
Convert Transcript Coordinates to Genomic Coordinates
'''

def index_cigar(genomic_start_pos, cigar):
    '''
    :param int genomic_start_pos: Start position of the transcript on the Genome
    :param str cigar: Cigar String representing the alignment of the transcript to the genome

    :rtype: intervaltree.IntervalTree obj
    :returns Interval Tree object which can be queried with transcript coordinates to get offset for calculating cooresponding genomic coordinate
    '''
    pattern_str = "([0-9]+)([A-Z])"
    cigar_pattern = re.compile(pattern_str)
    i = 0
    start  = 0
    stop   = 0
    offset = genomic_start_pos

    # interval tree obj
    cigar_index = intervaltree.IntervalTree()

    matches = re.findall(cigar_pattern, cigar)
    if not matches:
        raise Exception("Could not parse Cigar: {}. Expected patten: {}".format(cigar, pattern_str))

    for num_bases, cigar_char in matches:
        num_bases = int(num_bases)
        if cigar_char == "M" or cigar_char == "X" or cigar_char == "=":
            stop += num_bases
            cigar_index[start:stop] = (offset, False)
            start = stop

        elif cigar_char == "D" or cigar_char == "N":
            offset += num_bases

        elif cigar_char == "I" or cigar_char == "S":
            stop += num_bases
            cigar_index[start:stop] = (offset, True)
            offset -= num_bases
            start = stop

        else: # Unsupported cigar ops from SAM format: P and H (Padded and Hardclip)
            raise Exception("Unsupported Cigar Operation: {}. Only support M, X, =, D, N, I and S".format(cigar_char))

    return cigar_index


def query_cigar_index(transcript_pos, cigar_index):
    ''' Query cigar interval tree for transcript position
    :param int transcript_pos:
    :param intervaltree.IntervalTree cigar_index

    :rtype tuple
    :returns (interval_start, offset_for_genomic_position, is_on_a_insertion)
    '''
    result = cigar_index[transcript_pos]
    if len(result) == 0:
        raise Exception("Transcript coordinate is out-of-bounds. {}\t{}".format(chrom, line))
    elif len(result) > 1:
        raise Exception("Invalid cigar index. Transcript Coordinate is in two interval blocks. {}\t{}".format(line, res))

    res = next(iter(result))
    offset, is_ins = res.data
    start = res.begin
    return (start, offset, is_ins)

def get_insertion_representation(transcript_pos, start, offset):
    ''' Return genomic position spanning an insertion in the transcript
    :param int transcript_pos
    :param int start
    :param int offset

    :rtype str
    :returns  representation of the genomic position around the insertion context
    '''
    ins_start = max(0, start + offset - 1)
    return "I{}+{}".format(ins_start, transcript_pos - start + 1)


def run(infile1, infile2):
    ''' Main wrapper func
    :param str: infile1: Inputfile1 containing cigars for each transcript
    :param str: infile2: Inputfile2 containing coordinates to query
    '''
    print("\n{}: Started".format(datetime.datetime.now()))
    outfile = "outfile.txt"

    transcript_dict = {} # transcript -> (cigar_index, chrom)
    print("\n{}: Indexing Cigars".format(datetime.datetime.now()))
    with open(infile1, "r") as IN:
        for line in IN:
            transcript, chrom, genomic_start_pos, cigar = line.strip("\n").split("\t")
            assert transcript not in transcript_dict, "Duplicate Transcript encountered"
            cigar_index = index_cigar(int(genomic_start_pos), cigar)
            transcript_dict[transcript]  = (chrom, cigar_index)


    print("\n{}: Translating transcript coordinates to genomic. Writing to file: {}".format(datetime.datetime.now(), outfile))
    with open(infile2, "r") as IN, open(outfile, "w") as OUT:
        for line in IN:
            line = line.strip("\n")
            transcript, transcript_pos = line.split("\t")
            chrom, cigar_index = transcript_dict[transcript]
            transcript_pos = int(transcript_pos)

            start, offset, is_ins = query_cigar_index(transcript_pos, cigar_index)

            if is_ins:
                ins_query_representation = get_insertion_representation(transcript_pos, start, offset)
                OUT.write("{}\t{}\t{}\t{}\n".format(transcript, transcript_pos, chrom, ins_query_representation))
            else:
                OUT.write("{}\t{}\t{}\t{}\n".format(transcript, transcript_pos, chrom, transcript_pos + offset))
    
    print("\n{}: Finished\n".format(datetime.datetime.now()))

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
