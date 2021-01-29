


 """
Title
malaria.py
Date
2020-10-10
Author
Atte Oskari Räsänen
Description
The program creates a text file version of malaria FASTA file with protein names included
List of user defined functions
none
List of modules that are not explained in the course material
none
Procedure
1. Open the malaria blast file in read mode, create an empty dictionary called data. 
2. With for loop go over every line.  a receives the lines contents during every cycle. Every line is then split based on columns 
3. Select the first(ID) and 10th column (protein name). The file closes once all lines have been covered. 
4. Move on to the malaria fasta file by opening it the same way. We create two empty dictionaries: one for the contents of malaria and one for the header
5. Create an empty Header list. 
6. Iterate through all lines of the file, removing the line changes. Identify the lines with the ID(starts with >), remove >, split based on tabs and select the first column
   Save the result into the Header dictionary. Other results (sequences) save into the Fastadict.
7. Open malaria fasta file and an output file (this one in write-mode). Identify the ID lines just as before.
    If the ID matches the ID found in data dictionary, combine the ID with the protein name and append this into the Header2.
8. Iterate over the Header2 table, writing each line into the output file (along with line changes). 
    Remove >, split the lines and select the ID (first column). Find the equivalent from the Fastadict and
    write it into the output file. This is repeated until all the IDs have been covered. 
Usage
    python malaria.py malaria.fna malaria.blastx.tab output.txt
Bugs
    When I ran my code in the windows terminal, no output file was created but when I initially (before importing sys) ran
    the code on spyder, I got the text file along with correct answer. Taking care of this on Linux may have solved the issues. 
"""

#%%
import sys
fna = sys.argv[1]
blastx = sys.argv[2]
output = sys.argv[3]

with open(sys.argv[2], 'r') as malariablast:
    data={}
    for a in malariablast:
        lineList = a.split("\t")        
        if "null" not in lineList:             
            ids = lineList[0]           
            pname = lineList[9]         
            data[ids] = pname                           
with open(sys.argv[1], 'r') as malariafna:   
    Fastadict = {}  
    Header = {}
    Header2 = []                            
    for a in malariafna:
        a = a.replace('\n', '')        
        if a.startswith('>'):           
            fastaid = a.replace('>', '').split('\t')[0]
            Header[fastaid] = a 
        else:
            Fastadict[fastaid] = a
with open(sys.argv[1], 'r') as malariafna, open(sys.argv[3], 'w') as output:
    for b in malariafna:               
        if b.startswith('>'):           
            malaline=b.replace('>', '').replace('\n', '').split('\t')[0]                                       
                          
        if malaline in data:  
            h = Header[malaline]+'\tprotein='+data[malaline]
            Header2.append(h)                             
    for aline in Header2: 
        output.write(aline + '\n')
        aline=aline.replace('>', '')
        fileids=aline.split('\t')[0]
        output.write(Fastadict[fileids])


       


