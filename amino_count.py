'''
Title: amino_count.py
Date: 2020-10-28
Author: Atte Oskari Räsänen

Description
    The programme takes an amino acid sequence file as an input and calculates
    the absolute abundancies of each amino acid in an output file.

List of functions:
    None.

List of "non standard" modules:
    None.

Procedure
    1. Create a dictionary for the 20 standard amino acids(AA), making the single
    letter code as the key and the value as the abundance (0 initially).
    2. Iterate through the given input file, selecting lines containing the AA sequences.
    Create an inner iteration corresponding for each AA letter in the sequence line.
    3. If the letter does not correspond to the ones found in the AA dictionary,
    name it as 'X'. If the letter is one of the standard AAs, update the occurence
    in the AA dictionary.
    4. Once all the lines have been covered, iterate over the updated AA and write
    the key and the value into the output file. value needs to be converted into
    a string before writing it into the file.
Usage
    amino_count.py amino.faa amino.list output_file.txt


'''


import sys


with open(sys.argv[1], "r") as aminofile, open(sys.argv[2], "w") as output:
    #first create a dictionary list of amino acids, the letters acting as keys
    #and values as abundancies
    AA_list={'A':0, 'R':0, 'N':0, 'D':0, 'B':0, 'C':0, 'Q':0, 'E':0, 'Z':0, 'G':0, 'H':0, 'I':0, 'L':0, 'K':0, 'M':0, 'F':0, 'P':0, 'S':0, 'T':0, 'W':0, 'Y':0, 'V':0, 'X':0}
    for lines in aminofile:             #iterate through the amino acids
        if not lines.startswith('>'):     #skip the ID lines
            #print(lines)
            for amino in lines.upper().strip():   #make sure the letters are in upper character
                if not amino in AA_list:
                    amino = 'X'
                AA_list[amino] += 1       #count the abundance of the certain amino acid

    #print(AA_list)
    for key,value in AA_list.items():       #write them into the output file
        output.write(key + ' ')
        output.write(str(value) + '\n')
