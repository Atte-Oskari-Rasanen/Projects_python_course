# -*- coding: utf-8 -*-

#!/usr/bin/python3
'''
    Date: 2020-10-27
    Description:
        This program will take a user input and output either the full matrix, in case
        the user input is ALL, or the best match to the name given in the input.
    List of functions:
        None
    List of non-standard modules:
        None
    Procedure:
        In this script we use argparse to work from the command prompt. Then if the args.name
        is ALl the similarity input is considered as ALL and the full matrix is printed. If the
        args.name is not ALL it takes the name given, finds the best match (the highest number in my_split)
        and prints the name of the best match.
    Usage:
        score.py input_alignment_fasta_file output_file
'''

import argparse

parser = argparse.ArgumentParser(description="This is a description of what this program does")
parser.add_argument("input", help="name of input file")
parser.add_argument("name", help="the name to search")
args = parser.parse_args()

input_file = args.input
if args.name == 'ALL': #Checks if the name in args.name is ALL
    similarity_input = "ALL" #If it is the similarity input is ALL
else: #If it's not all
    similarity_input = args.name #It takes the name as similarity input

with open(input_file, 'r') as matrix:  # open the file

    if similarity_input == 'ALL': #If the input is all
        for element in matrix: #Prints the full matrix found in the input file
            print(element)
    else: #If the input is something other than all
        head_list = matrix.readline().strip().split("\t") #A head_list is created in which the first row is found
        for row in matrix: #It goes through the matrix row by row
            if similarity_input.capitalize() in row.capitalize(): #If the similarity input is in the row
                my_split = row.split("\t") #It splits the row through tabs and assigns it to my split
                my_max = max(my_split[1:]) #It takes my split from index one, finds the higest number in it and assigns it to my max
                my_index = my_split.index(my_max) #it finds the index of my max in my split and assigns it to my index
                best_name = head_list[my_index] #Finds the item in the same index as in my_index
                print(f"The best match is {best_name} with a score of {my_max}") #It prints the best name and the score
