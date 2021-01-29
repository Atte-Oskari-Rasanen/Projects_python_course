'''
Title: Project.py
Date: 2020-30-10
Author: Atte Oskari Räsänen

Description
    This programme takes two input files; Fasta file and a file with enzymes along 
    with their recognition sequences. The file has to be in the format described
    in the readme file. The programme asks for the enzyme with which to cleave the
    sequences and outputs a file with the Header ID and number of fragments on the
    first line and on the second the fragments. 
List of functions: 
    def Create_dict(fasta):
    def Choose_enzyme(selection):

List of "non standard" modules:
    None.

Procedure
    1. Create a function inside which you open the DNA sequence file and read it into
    a dictionary. Define idheader as None before iteration. When iterating, if the line
    starts with '>', save it as a the idheader. If the line does
    not start with '>', append it into a list. To include the header and sequence,
    repeat the procedure for creating the ID key as mentioned above. Return the dictionary.
    2. Define the function's dictionary outside the function based on the input file.
    3. Create a function in which the user can select the enzyme he/she wants to use.
    Extract the recognition site from the enzyme and return it.
    4. Find the index of the recognition site at which the cut happens. In the recognition
    site this is marked as '/', e.g AA/GC. Save this index into a variable (cut_site in the code).
    5. Remove recognition site's '/' and save into a variable (Recog_seq in the code).
    6. Iterate over the fasta dictionary created in the first function. Save the key
    as the ID variable. If recognition site is found in the value(sequence) of the dictionary,
    proceed. Define a variable pos(position) as -1. Add 0 into an empty list (List_index)
    which will later be filled with other values.
    7. While the statement of pos is true, add +1 to pos, making it 0 and redefine 
    pos as the site of the value(seq), making this being read as index 0(the starting 
    point recognition site=. Stop if the pos equals -1
    8. To the pos, add the index at which the cut happens and append this into the previously 
    defined list(List_index) and lastly add the total length of the value(sequence)
    in which the recognition site(s) was found. 
    9. Next, iterate over the created list. Define pos1 as the position i in the List_index
    and pos2 as the position next to this. Extract the sequence fragment from the 
    fasta dictionary value using these two pos' as the window. Append into a fragments
    list. 
    10. Make the fragments list into a string where the items are tab separated. use
    join method for this. Add the fragments string into the fragments dictionary where
    the ID is the fasta file's header. 
    11. Count the number of fragments in the fragments list and save it into a variable.
    12. Turn this variable into a string and add into a dictionary containing header ID
    as the key and fragment count as the value.
    13. Write into the output file where the first line contains the header and the 
    fragment count, the second line containing the fragments. 
Usage
    python Project.py regions.fna NE_Biolabs_Endonucleases.txt output_regions.txt
    python Project.py testfile.fna NE_Biolabs_Endonucleases.txt output_testfile.txt

Extra notes
    The example output files where generated using the enzyme BfaI from the NE_Biolabs_Endonucleases.txt
'''
import sys
with open(sys.argv[1], "r") as InputFile:
    def Create_dict(fasta):
        DNA_dict={}
        lista = []
        idheader = None
        for lines in fasta: 
            if lines.startswith('>'):
                if idheader:
                    DNA_dict[idheader]=''.join(lista)
                    #print(lista)
                    del lista[:]
    
                idheader = lines.strip().split()[0]  #saves the id. cut the lines based on where id starts, the the first
                #print(idheader)
            else:
                lista.append(lines.strip())   ##saves the sequences (multiline)
                
                #print(lista)
        DNA_dict[idheader]=''.join(lista)       #you need to add this once more as the last id+seq ends with a seq instead of the next id
        del lista[:]    
        return DNA_dict
    
    Fasta_dict= Create_dict(InputFile)
    
    def Choose_enzyme(selection):          #choose the enzyme from the given input file
        with open(selection) as enz_file:     #open the file

            enzymes=input("Please enter the enzyme name: ")  #enter the name 
            for row in enz_file:           #search the name from the file
                if enzymes in row:
                    row=row.split()         #split based on the space
                    Recog_site=row[0]           #take the first item which is the recognition site
                    Recog_site=str(Recog_site)  #make into a string

            return Recog_site
   
    
    Recog_seq=Choose_enzyme(sys.argv[2])    #call the function Choose_enzyme
            
    cut_site=int(Recog_seq.find('/'))        #get the cut site which is noted by / in the sequence
    Recog_seq=Recog_seq.replace('/','')    #remove '/' to be able to find the location from the target sequence


                
    Dict_fragments={}   #create an empty dictionary for the final fragments
    Fragments_list=[]   #used for storing the fragments
    List_index=[]
    ID_and_frag_count={}    #ID header acts as the key, frag_count as the number of fragments in a sequence

    for key,value in Fasta_dict.items():        #iterate over the fasta dict
        ID=key
        #print(ID)
        if Recog_seq in value:     #we want to split the sequence(value) in the spots where the recog.site is
                                #and split it based on the nth letter in these sites

            pos=-1              #define pos(position) initially as -1
            List_index=[0]         #0 into the list
            while True:                             #getting the indexes where each recog.site starts
                pos=value.find(Recog_seq, pos + 1)     #add +1 to pos, making it 0 
                #and redefine pos as the site of the value(seq), making this being read as index 0(the starting point recognition site)

                if pos ==-1:            #stop while loop if pos gets to -1
                    break
                #print(pos)
                pos+=cut_site               #to the pos, add the index at which the cut happens
                #print(pos)
                List_index.append(pos)              #append these indexes into a list
                #print(List_index)
            List_index.append(len(value))        #append the total length of the sequence
            #print(List_index)
            
        for i in range(len(List_index)-1):       #iterate over the aforementioned index
            pos1=List_index[i]              #pos1 is the item of the list at the i
            #print(pos1)
            pos2=List_index[i+1]            #pos2 is the item of the list at the position next to i
            #print(pos2)
            seq_frag=value[pos1:pos2]       #sequence fragment is between these
            Fragments_list.append(seq_frag)  #append this into a list


        Fragments='\t'.join(str(subfrag) for subfrag in Fragments_list).upper()  #make the list into a string, make sure all NTs are upper case
        Dict_fragments[ID]=Fragments        #add these strings tab separated strings into the dictionary
        frag_count=len(Fragments_list)     #caclulate the length of the fragments
        ID_and_frag_count[ID]= str(frag_count)  #add this into a dictionary containing the header IDs
        Fragments=''                     #empty the variables after a loop
        Fragments_list=[]
        
with open(sys.argv[3],"w") as out1:   #write the results into the outputfile

    for a,b in zip(Dict_fragments, ID_and_frag_count):      #perform parallel iteration
        out1.write(a + ' Number of fragments: ' + ID_and_frag_count[b] + '\n') #first line should contain the ID and the number of fragments
        out1.write(Dict_fragments[a] + '\n')        #second line should contain the fragments
        
