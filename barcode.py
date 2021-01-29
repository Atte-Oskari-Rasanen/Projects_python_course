'''
Title: barcode.py
Date: 2017-10-17
Author: Atte Oskari Räsänen

Description
    The programme takes a FASTQ file, extracts the ID headers and the sequences,
    writes the ones containing certain barcodes into one output and the rest
    of the IDs and sequences into a separate one.

List of functions:
    None.

List of "non standard" modules:
    None.

Procedure
    1. Initially, a general dictionary should be made which contains all sequence
    headers (as keys) and their corresponding DNA sequences. The rest of the
    content is ignored. This is done based on line numbers.
    2. Open the FASTQ file and read it line by line to make it into a list and
    save it as a variable.
    Iterate over the lines of the file, keeping count of the lines with a separate
    variable (). Every fourth line is the last line of the given sequence (quality
    score line). To get the header with the ID, subtract 3 from the number of lines
    counted (which is 4 for every sample) and for the sequence subtract 2. Read
    these into a general dictionary.
    3. Once the general dictionary has been generated, iterate over it, checking
    for the presence of barcodes at the start of each line and slice them off,
    after which save this into the given sample dictionary. There should be 3 sample
    dictionaries in total as there are three barcodes. If the sequence does not start
    with a barcode, save it into a separate dictionary (NoBarc in the code).
    4. Open a file in a write mode into which the ids+sequences containing barcodes
    will be written and another file into which the rest of the ids+sequences will be
    placed. Line changes should be included in each case.
Usage
    barcode.py input_filename output_file.txt

'''

import sys

with open(sys.argv[1], "r") as fastqfile:
    fastqfile=fastqfile.readlines()
    general_dict={}
    sample1_dict={}
    sample2_dict={}
    sample3_dict={}
    NoBarc={}
    Number_of_lines=0
    genseq_List=[]
    for lines in fastqfile:             #go over the fastqfile line by line
        Number_of_lines +=1         #keep count of the lines
        if Number_of_lines % 4 ==0:     #once you reach the 4th(quality score line, extract keys and values)
            ID_line= (Number_of_lines -1) -3   #find the line with the ID
            seq_line = (Number_of_lines -1) -2  #find the line with the sequence
            idkey=fastqfile[ID_line].rstrip()
            seq=fastqfile[seq_line].rstrip()
            if not ID_line in general_dict:
                general_dict[idkey]=seq
            else:
                print("Error: found two IDs that were similar")

    #now iterate through the general_dict
    for key,value in general_dict.items():
        if value.startswith('TATCCTCT'):
            seqmod1=value[8:]                      #get rid of the barcode, save it into the sample_dict
            sample1_dict[key]=seqmod1
        elif value.startswith('GTAAGGAG'):
            seqmod2=value[8:]
            sample2_dict[key]=seqmod2
        elif value.startswith('TCTCTCCG'):
            seqmod3=value[8:]
            sample3_dict[key]=seqmod3
        else:
             NoBarc[key]=value


    #write the sequences which contained barcodes into their own file and the rest of the sequences to a separate file
with open (sys.argv[2], "w") as output_samplefiles, open(sys.argv[3], 'w') as undetermined:
    for key,value in sample1_dict.items():
        output_samplefiles.write(key + '\n')
        output_samplefiles.write(value + '\n')
    for key,value in sample2_dict.items():
        output_samplefiles.write(key + '\n')
        output_samplefiles.write(value + '\n')
    for key,value in sample3_dict.items():
        output_samplefiles.write(key + '\n')
        output_samplefiles.write(value + '\n')
    for key,value in NoBarc.items():
        undetermined.write(key + '\n')
        undetermined.write(value + '\n')
