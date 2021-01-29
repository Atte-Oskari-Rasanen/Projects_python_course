'''
Title: MSA.py
Date: 2017-27-10
Author: Atte Oskari Räsänen

Description
    This code takes in a file with subject names as well mtDNA and Y chromosomal
    data. It parses through the file, extracting the mtDNA and Y chromosomal data
    into their own output files and subsequently performs a multiple sequence 
    alignment and calculates percentage identity on the output files, generating
    final files in the format name1 - name2 - percentage identity- score.
    
List of functions: 
    def gen_data_parser(fh)
    def scoring(seq1, seq2, header1, header2)
    def dict_score_loop(dict)
    def un_wrapper(my_list)
    def printer(my_list_of_lists, output_file)
    


List of "non standard" modules:
    None.

Procedure
    1. Parse the file into mtDNA and Ychrom dictionaries using gen_data_parser. 
    The function iterates over the input line(in this case GeneticData.txt), gets
    the name of the line (the subsequent lines belong to this subject), identifies
    whether the subject has mtDNA and/or Ychromosomal sequences based on the names.
    2. If a line found with mtDNA or Y chromosome name, jump to next line, which is
    the sequence, strip it from new line characters and save into the dictionary
    with the corresponding name. Return the dictionaries.
    3. Subsequent scoring function takes in the headers(names) and their sequences(seq)
    and iterates over them, applying parallel iteration. Score the corresponding
    nucleotides(or empty spaces) according to the values used and take into account
    transitions and transversions. Refer to the code for more details.
    4.At the end of this function, calculate the alignment score and percentage
    identity score. Return the scores.
    5. Dict_score_loop function is for saving the name pairs with the previous 
    function's scores by iterating over the dictionary and within each iteration iterates over the
    dictionary again (getting the name pairs). We call the scoring function here 
    and save the name1-name2 pair of the iteration with the corresponding 
    dictionary values that were  calculated in the previous function. This is 
    saved into a list 'Row', which is appended into the empty list. This is done
    for each name(key) pair. Return the empty list. 
    6. un_wrapper function eliminates duplicates of name1-name2 combinations
    by iterating over the rows of the list, seeing if the index 1(name1) corresponds
    to the index 2(name2). Refer to the code for more details
    7. In the last function, printer, input the lists with the duplicates removed,
    add the header line first with the required titles, loop over the list contents
    and then print them into the output file in the wanted order.
    
    
    
Usage
    Cannot enter a fasta file as the given file for the exercise (GeneticData.txt)
    is not in fasta format. The input file GeneticData.txt is already written 
    into the code, the user only needs to download the file. 


'''
def gen_data_parser(fh):
    with open(fh) as f:     
        mt_dict = {}
        y_dict = {}
        my_count = 0
        for line in f:      #go through the file line by line
            if my_count == 0:          #the first file contains the name
                name = line.strip()
                my_count+=1
                continue
            if line.startswith("mtDNA"):    #saving mtDNA seq into dict
                seq = next(f)           #the next line after mtDNA is the seq
                edit = seq.strip()      #strip from new line characters
                mt_dict[name] = edit      #save into dict with the corresponding name
                seq = ''          #empty the name and seq variable for the next loop
                edit = ''
                continue
            if line.startswith("Y chromosome"):    #repeat the same steps as when line starts with mtDNA
                seq = next(f)
                edit = seq.strip()
                y_dict[name] = edit
                seq = ''
                edit = ''
                continue
            if "hemophilia" in line:       #skip lines containing the word hemophilia
                continue
            if line.strip() != '':          #after stripping the file from new lines, if the line is not empty, this line has name
                name = line.strip().strip(">")
                continue
        return mt_dict, y_dict      #return the dictionaries for mtDNA and Y chromosome


def scoring(seq1, seq2, header1, header2): #to check the alignment score
    transition = ['AG', 'TC', 'GA', 'CT'] #the transition scores
    score_list = [] #a list to add all the scores when they are calculated
    identical_nt_list = []
    for a, b in zip(seq1, seq2): #this for loop goes through the nucleotides in seq1 and seq2 (a and b respectively) and compares them
        a = a.upper() #in case the sequence is lower case
        b = b.upper() #in case the sequence is upper case
        if a == b: #if the nucleotide in a is the same as in b
            if '-' in a and '-' in b: #if both are - the score is 0 and is added to the list score_list
                score = 0
                score_list.append(score)
            if '?' in a and '?' in b:
                score = -1
                score_list.append(score)
            else: #if they are not - they are nucleotides and the score is 1, which is added to the score_list
                score = 1
                score_list.append(score)
                identical_nt_list.append(score)
        else: #if the nucleotides are not the same
            if '-' in a or '-' in b: #if one sequence has a gap when aligned to the other one the score is -1 and is added to score_list
                score = -1
                score_list.append(score)
            if '?' in a or '?' in b:
                score = -1
                score_list.append(score)
            else:
                nt_sum = a + b
                if nt_sum in transition: #if they are in transition the score is -1 and it is appended to score_list
                    score = -1
                    score_list.append(score)
                else: #if it is not a transition it means that it is a transversion, which has a score of -2 that is added to score_list
                    score = -2
                    score_list.append(score)
    align_score = str(sum(score_list)) #then i do a sum to sum all the values in the score_list and assign it to a variable score_sum
    id_score = str(sum(identical_nt_list)/len(seq1)*100)
    return id_score, align_score #id and alignment must be returned


def dict_score_loop(dict):
    my_list = []
    for k1, v1 in dict.items():  #iterate over the dictionary
        for k2, v2 in dict.items():   #iterate over it within the first iteration
            id_score, align_score = scoring(v1, v2, k1, k2)  #call the scoring function
            row = [k1, k2, id_score, align_score]   #create list with the names and the scores
            my_list.append(list(row))   #append this into the list
    return my_list
mt_dict, y_dict = gen_data_parser("GeneticData.txt")    #
mt_score_out = dict_score_loop(mt_dict)    #mtDNA scores in a list
y_score_out = dict_score_loop(y_dict)    #Ychrom scores in a list
#print(mt_score_out)
#print(y_score_out)

#we send the contents from prev function to the un_wrapper
def un_wrapper(my_list):    #to get the rid of the rows with duplicates
    check_list = []
    new_list = []
    for row in my_list:    #loop over the rows that become the matrix in second part
        #print(row[0])
        if row[0] != row[1]:     #if the name1(index 1 of the row) does not equal to name2
            if row[1] not in check_list:    #if the name2 is not in checklist 
                new_list.append(row)  #append it into the new list. Thus, if there are duplicates, dont append
        check_list.append(row[0])   #append name1 to the check list


    return new_list

mt_unwrapped = un_wrapper(mt_score_out)  #call the functions
y_unwrapped = un_wrapper(y_score_out)

def printer(my_list_of_lists, output_file):    #write into the outputfile
    header = "SampleA\tSampleB\tIdentityScore\tScore"  #use as the top row, columns are tab separated
    with open(output_file, 'w') as f:       #open the output file in writemode
        print(f"{header}", file=f)
        for row in my_list_of_lists:    #go through the list, print out the values in columns
            print(f"{row[0]}\t{row[1]}\t{row[2]}%\t{row[3]}\t", file=f)
printer(mt_unwrapped, "mt_test.txt")
printer(y_unwrapped, "y_test.txt")
