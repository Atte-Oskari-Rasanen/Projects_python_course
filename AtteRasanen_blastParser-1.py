# -*- coding: utf-8 -*-
'''
Title: AtteRasanen_blastParser
Date: 26.1.2021
Author: Atte Räsänen

Description:
    Script for parsing a blast file and giving an output file file with the information
    about the queries and the related information including target, e-value, identities
    and score.

Usage:
    python AtteRasanen_blastParser.py yeast_vs_Paxillus.blastp
Arguments:
     yeast_vs_Paxillus.blastp: file containing the blastp results

List of functions:
    None

Exceptions:
    None

'''
import pandas as pd 
import sys

with open(sys.argv[1], "r") as file: #open the file
    parsingfile = file.readlines()   # read lines into a list

# create empty lists for sections we want to extract from the file
    querylist = []
    targetlist = []
    E_vallist = []
    identitylist = []
    scorelist = []
    
    
    for lines in parsingfile: #go over the file line by line

        if "Query=" in lines:   #get the query
         #   print(lines)
            querypart = lines.split(' ')[1]

        if lines.startswith(">"):     #get the target
             targetpart = lines.replace(">", "") #remove the > as well as line changes
             targetpart = targetpart.replace("\n","")

        if "***** No hits found" in lines:  #if no hits for the query, append empty space into lists
            querylist.append(querypart)
            targetlist.append(" ")
            E_vallist.append(" ")
            identitylist.append(" ")
            scorelist.append(" ")
            
        if lines.startswith(" Score"):   #get the score
            scorepart = lines.split(' ')[3]

        if "Expect = " in lines:    #get the e value
            E_valpart = lines.split(' ')[9]
            E_valpart = E_valpart.replace(",", "")

        if "Identities " in lines:   #get the identity score
            #this is the last value we need to find from a single query
            # so after this all values are appended to their correspondent lists
            identitypart = lines.split(' ')[4]
            identitypart = identitypart.replace("(", "")
            identitypart = identitypart.replace(")", "")
            identitypart = identitypart.replace(",", "")
        #appending values of a single query
            querylist.append(querypart)
            targetlist.append(targetpart)
            identitylist.append(identitypart)
            scorelist.append(scorepart)
            E_vallist.append(E_valpart)
            
    df = pd.DataFrame()  #create an empty dataframe and add the values
    df["#query"]=querylist
    df["target"]=targetlist
    df["e-value"]=E_vallist
    df["identity(%)"]=identitylist
    df["score"]=scorelist
    
with open("OutputFile.txt", "w") as out: #write the dataframe into an output file
#write the dataframe into an output file
    df.to_string(out, index=None)

        
        