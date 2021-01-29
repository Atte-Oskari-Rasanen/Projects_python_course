#!/usr/bin/python3
"""
    Description:
        This program will parse a file with a similarity matrix in squareform, then generate a dendrogram from it
        that can be interpreted as a phylogenetic tree.
    List of functions:
        main(): calls all functions for all files in list
        matrix_parser(): takes a file, makes a numpy array and returns it with a list of labels
        make_tree(): will take a file, send it to matrix_parser() and then generate a tree using dendrogram
        and output the tree as a file

"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist

inputs = ["mt_out_identity", "mt_out_align", "y_out_identity", "y_out_align"] #list of file names that are output from step2


def matrix_parser(fh):
    """
    :param fh: input file
    :return: list of headers and a numpy array
    """
    my_counter = 0
    m = []
    with open(fh) as f:
        for line in f:
            if my_counter == 0:
                names = list(line.strip("\n").split("\t"))
                my_counter += 1
            else:
                row = list(line.strip("\n").split("\t"))
                m.append(row[1:])
    matrix = np.array(m, dtype=float)
    return names, matrix


def make_tree(fh):
    """
    :param fh: input file
    :return: nothing, writes to files
    """
    names, matrix = matrix_parser(fh)
    distance_matrix = pdist(matrix) #dendrogram does not like squareform
    """
    Hierarchical clustering
    """
    linked = linkage(distance_matrix, 'single') #linkage analysis
    label_list = names[1:] #will slice the header to remove leading tab
    plt.figure(figsize=(25, 15)) #define a size of the figure to be plotted
    dendrogram(linked, #take results from linkage
               orientation='left', #root the tree to the left
               labels=label_list, #give leaves names from headers
               distance_sort='descending', #start by drawing most distant branches
               show_leaf_counts=True) #showing counts on leaves in complex trees
    thing = fh + ".png" #generate a friendly filename from the ingoing file
    plt.savefig(thing) #output to file


def main():
    """
    main function
    :return: one .png tree file per ingoing matrix
    """
    for f in inputs: #for all items in the list of files going in
        make_tree(f) #do everything else on all the files

main() #calling main
