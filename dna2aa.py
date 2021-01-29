'''
Title: dna2aa.py
Date: 2017-10-17
Author: Atte Oskari Räsänen

Description
    The programme takes a multi fasta DNA sequence file as an input, converts the sequence
    into RNA (T->U). The RNA sequence is translated into a protein sequence
    based on the RNA codons which correspond to a single letter amino acid code.
    The subsequent sequence is saved into an output file.

List of functions:
    DNA(inputfile)
    RNA()

List of "non standard" modules:
    None.

Procedure
    1. Convert the T bases of the AA dictionary (keys=DNA codons, values=AA letters)
    by iterating over the keys and replacing the Ts with Us. Save into an RNA dictionary
    2. Create a function inside which you open the DNA sequence file and read it into
    a dictionary. Define idheader as None before iteration. When iterating, if the line
    starts with '>', save it as a the idheader and remove the '>'. If the line does
    not start with '>', append it into a list. To include the header and sequence,
    repeat the procedure for creating the ID key as mentioned above. Return the dictionary.
    3. Define the function's dictionary outside the function based on the input file.
    4. Create a function for creating RNA dictionary. Based on the DNA dictionary,
    create an RNA one by iterating through the DNA dictionary and replacing T bases
    with Us. Return RNA dictionary
    5. Define the RNA dictionary outside the function (in the code called RNA_table)
    6. Iterate through the keys of the RNA_table, saving the corresponding sequence
    into a variable. Then within this iteration, go through the current RNA sequence
    using a window of three letters. Save the windows as codons. Correspond the
    codon to the amino acid from the RNA dictionary created in step 1. Only include
    codons that are 3 letters long. Create a protein dictionary with IDs as keys and
    protein sequences as values.
    7. Open the outputfile in write mode, iterate over the protein dictionary
    created in the previous step, add "> rf 1 " to the key (define this as e.g. a header)
    and then write the header and the value into the file, including line changes.

Usage
    dna2aa.py DNA.faa output_file.txt



'''


import sys


#write a separate function for reading fasta file into a dictionary

#make tehe
AA_dict = []
# Aminoacid table
AA_table_withDNA = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

# Creating a RNA table, replacing T with U
RNA_dict={}
for old_key in AA_table_withDNA:
    if 'T' in old_key:
        new_key= old_key.replace('T','U')
        RNA_dict[new_key] = AA_table_withDNA[old_key]
    else:

        RNA_dict[old_key]=AA_table_withDNA[old_key]


#read the sequence file into a DNA dictionary
def DNA(inputfile):
    DNA_dict={}
    lista = []
    with open(inputfile) as file1:
        idheader = None
        for lines in file1:
            if lines.startswith('>'):
                if idheader:
                    DNA_dict[idheader]=''.join(lista) # adding the id to the dic, making the whole 3 lines into a one line, then deleting the lista content so we can start again
                    del lista[:]

                idheader = lines.strip().replace('>','').split()[0]
            else:
                lista.append(lines.strip())
                #print(lista)
        DNA_dict[idheader]=''.join(lista)       # need to add this once more
                                                #as the last id+seq ends with a seq
                                                #instead of the next id
        del lista[:]

    return DNA_dict

#make an RNA dictionary based on the DNA dictionary
def RNA(DNAseq_dict):
    RNA_dict={}
    for key,value in DNAseq_dict.items():
        #print(key,value)
        RNA_dict[key]=value.replace('T','U')
    #print(RNA_dict)
    return RNA_dict


DNAseq_dict=DNA(sys.argv[1])
RNA_table = RNA(DNAseq_dict)

Protein_dict ={}
codon_list=[]
AA_list=[]

#read RNA as codons and match to the AA from the list, making a dictionary for protein sequences
for singleID in RNA_table.keys():
    rnaseq=RNA_table[singleID]
    for index in range(0, len(rnaseq), 3):
        codon = rnaseq[index:index+3]
        AA = RNA_dict.get(codon)
        letters=len(codon)
        if letters == 3:
            AA_list.append(AA)

    AA=''.join(AA_list)
    Protein_dict[singleID]=AA
    AA_list = []


#iterate through the protein dictionary and write it into an output file
with open (sys.argv[2], "w") as output:
    for key,value in Protein_dict.items():
        Header = "> rf 1 " + key
        output.write(Header + '\n')
        output.write(value + '\n')
