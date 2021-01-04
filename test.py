# stdlib
import unittest
import os
import subprocess

# this project
import translate_coordinates

# 3rd party
from intervaltree import IntervalTree, Interval

class TestCmdlineExample(unittest.TestCase):
    ''' Test with a mock command line run
    '''
    def setUp(self):
        '''
        '''
        pass

    def test_cmd_line(self):
        '''
        '''
        subprocess.check_output("python3 translate_coordinates.py infile1.txt infile2.txt", shell=True)
        subprocess.check_output("diff outfile.txt outfile.expected.txt", shell=True)


class TestCigarParsing(unittest.TestCase):
    ''' Test cigar parsing
    '''
    def setUp(self):
        '''
        '''
        pass

    def test_cigar_parsing(self):
        ''' Check if we can parse cigars appropriately
        '''
        # simple parse with supported cigar operators
        cigar = "10S20M30M40D20I40M50D100X10M100N50M"
        cigar_index = translate_coordinates.index_cigar(3, cigar)
        self.assertEqual(cigar_index is not None, True)

        # character not supported - should fail
        cigar = "10M20S30M40P"
        cigar_index = None
        try:
            cigar_index = translate_coordinates.index_cigar(3, cigar)
        except Exception as e:
            pass
        self.assertEqual(cigar_index, None)

    def test_parsed_interval(self):
        ''' Check parsed interval from cigar
        '''
        cigar = "4M8I5M"
        #         0123456        789 11
        # Genome: ---ACGT--------ACGTG
        # Trans:     ACGTGGGGGGGGACGTG
        #            0123456789      16
        cigar_index = translate_coordinates.index_cigar(3, cigar) # start pos 3
        self.assertEqual(cigar_index, IntervalTree([Interval(0, 4, (3, False)), Interval(4, 12, (3, True)), Interval(12, 17, (-5, False))]))

        cigar = "2I3D4M1D" # leading insertion and end with deletion
        #           0123456
        # Genome: --GTGACGT--
        # Trans:  AC---ACGT--
        #         01   2345
        cigar_index = translate_coordinates.index_cigar(0, cigar) # start pos 0
        self.assertEqual(cigar_index, IntervalTree([Interval(0, 2, (0, True)), Interval(2, 6, (1, False))]))

        cigar = "3D2I3D4I" # all transcript intervals are on Insertions
        #         012  345
        # Genome: GTG--ACG---
        # Trans:  ---GG---ACGT
        #            01   2345
        cigar_index = translate_coordinates.index_cigar(0, cigar) # start pos 0
        self.assertEqual(cigar_index, IntervalTree([Interval(0, 2, (3, True)), Interval(2, 6, (4, True))]))


class TestIntervalQuery(unittest.TestCase):
    def setUp(self):
        '''
        '''
        self.transcript_pos = None
        self.genome_pos     = None
        self.cigar_index    = None
        self.ins_repr       = None

    def test_insertion_query(self):
        ''' test transcript pos on insertion
        '''
        cigar = "3D2I3D4I" # all transcript intervals are on Insertions
        #         012  345
        # Genome: GTG--ACG---
        # Trans:  ---GG---ACGT
        #            01   2345
        self.cigar_index = translate_coordinates.index_cigar(0, cigar)

        self.transcript_pos = 1
        self.ins_repr = "I2+2"
        self.helper_ins_query()

        self.transcript_pos = 5
        self.ins_repr       = "I5+4"
        self.helper_ins_query()

    def helper_ins_query(self):
        ''' Helper func for running test queries spanning insertion sites
        '''
        start, offset, is_ins = translate_coordinates.query_cigar_index(self.transcript_pos, self.cigar_index)
        self.assertEqual(is_ins, True)
        ins_representation = translate_coordinates.get_insertion_representation(self.transcript_pos, start, offset)
        self.assertEqual(ins_representation, self.ins_repr)
        
    def helper_non_ins_query(self):
        ''' Helper func for running test queries for non insertion case
        '''
        start, offset, is_ins = translate_coordinates.query_cigar_index(self.transcript_pos, self.cigar_index)
        self.assertEqual(is_ins, False)
        self.assertEqual(self.transcript_pos + offset, self.genome_pos)
        
    def test_transcript_ahead_of_genome(self):
        ''' Weird alignment - Technically still valid though 
        '''
        cigar = "4M8I5M"
        #         0123456        789 11
        # Genome: ---ACGT--------ACGTG
        # Trans:     ACGTGGGGGGGGACGTG
        #            0123456789      16
        self.cigar_index = translate_coordinates.index_cigar(3, cigar) # start pos 3
        self.transcript_pos = 16
        self.genome_pos     = 11
        self.helper_non_ins_query()

    def test_standard_transcript_alignment(self):
        ''' Some standard alignments
        '''
        cigar = "2M4D2M2I2M"
        #         01234567  89
        # Genome: ACGGGGTC--CA
        # Trans:  AC----TCGGCA
        #         01    234567
        self.cigar_index = translate_coordinates.index_cigar(0, cigar) # start pos 0

        self.transcript_pos = 7
        self.genome_pos     = 9
        self.helper_non_ins_query()

        self.transcript_pos = 3
        self.genome_pos     = 7
        self.helper_non_ins_query()

        cigar = "3M4D3M"
        #         012345678910
        # Genome:  ACGTTTTACT
        # Trans:   ACG----ACT
        #          012    345
        self.cigar_index = translate_coordinates.index_cigar(1, cigar) # start_pos 1
        
        self.transcript_pos = 5
        self.genome_pos     = 10
        self.helper_non_ins_query()

        self.transcript_pos = 3
        self.genome_pos     = 8
        self.helper_non_ins_query()

        self.transcript_pos = 2
        self.genome_pos     = 3
        self.helper_non_ins_query()

        self.transcript_pos = 0
        self.genome_pos     = 1
        self.helper_non_ins_query()


        cigar = "2I3D4M1D" # leading insertion and end with deletion
        #           0123456
        # Genome: --GTGACGT--
        # Trans:  AC---ACGT--
        #         01   2345
        self.cigar_index = translate_coordinates.index_cigar(0, cigar) # star pos 0

        self.transcript_pos = 5
        self.genome_pos = 6
        self.helper_non_ins_query()

        self.transcript_pos = 2
        self.genome_pos = 3
        self.helper_non_ins_query()

        self.transcript_pos = 0
        self.ins_repr = "I0+1"
        self.helper_ins_query()

        self.transcript_pos = 1
        self.ins_repr = "I0+2"
        self.helper_ins_query()



if __name__ == '__main__':
    unittest.main()
