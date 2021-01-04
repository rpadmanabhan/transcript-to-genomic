# transcript-to-genomic

##### An excercise question I solved for converting transcript coordinates to genomic coordinates. Run it as :
```
python3 translate_coordinates.py infile1.txt infile2.txt > outfile.txt
```

Where the input and output files are :

infile1.txt : A four column (tab-separated) file containing the transcripts. The first column is the transcript name, and the remaining three columns indicate it’s genomic mapping: chromosome name,
0-based starting position on the chromosome, and CIGAR string indicating the mapping.

infile2.txt : A two column (tab-separated) file indicating a set of queries. The first column is a transcript name, and the second column is a 0-based transcript coordinate.

outfile.txt : A four column tab separated file with one row for each of the input queries. The first two columns are exactly the two columns from the second input file, and the remaining two columns are the chromosome name and chromosome coordinate,
respectively.

##### Run tests as :
```
python3 test.py -v
```
