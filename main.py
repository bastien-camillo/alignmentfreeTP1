#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@teacher: yoann dufresne
@students: camillo and izard
"""

from loading import load_directory
from kmers import stream_kmers
from alive_progress import alive_bar
import argparse

parser = argparse.ArgumentParser() # create an argument parser
parser.add_argument('--dir',   dest='data_directory', type=str, required=True, help='path to the data directory')
parser.add_argument('--k',   dest='k', type=int, required=True, help='kmer size')
args = parser.parse_args() # parse the arguments


def similarity(A : int, inter : int, B : int):
    """
    Calculate the similarity between two sequences from the
    number of kmers in each sequence and the number of kmers in common
    ------------
    parameters
    A : size of the first set
    inter : size of the intersection between the two sets
    B : size of the second set
    ------------
    outputs
    the coverage of the intersection on A
    the coverage of the intersection on B
    """    
    return inter/(A), inter/(B)

def jaccard(A : int, inter : int, B : int):
    """
    Calculate the jaccard index, the intersection over the union of two sets
    It quantify the coverage of the intersection over the union of two sets
    ------------
    parameters
    A : number of kmers in the first sequence
    inter : number of kmers in the intersection
    B : number of kmers in the second sequence
    ------------
    output
    jaccard index
    """
    return inter/(A - inter + B)

def intersection(kmers_A : list, kmers_B : list, bar) -> int:
    '''
    Return the number of kmers in common between A and B
    ------------
    parameters
    kmers_A : a list of kmers from sequence A
    kmers_B : a list of kmers from sequence B
    ------------
    outputs
    inter_A : the number of kmers in kmers_A that are also in kmers_B
    inter_B : the number of kmers in kmers_B that are also in kmers_A
    '''
    
    dico1, dico2 = {}, {}
    inter_kmers = set()
    inter, value = 0, 0

    # create a dictionary with the kmers as keys and the number of occurences as values
    for kmer in kmers_A:
        bar()
        if kmer not in dico1:
            dico1[kmer] = 1
        else:
            dico1[kmer] += 1
    
    # seek if there is a kmer in kmers_B that is also in kmers_A
    for kmer in kmers_B:
        bar()
        if kmer in dico1:
            inter_kmers.add(kmer)
            # if the kmer is also present add the number of occurences to the value
            if kmer not in dico2:
                dico2[kmer] = 1
            else:
                dico2[kmer] += 1
    
    # count the number of kmers in common between A and B
    for kmer in inter_kmers:
        value += min(dico1[kmer], dico2[kmer])
    inter = value

    return inter

def main( data_directory : str, k : int):
    '''
    Compare the similarity between all the sequences in a directory
    ------------
    parameters
    data_directory : path to the directory containing the sequences
    k : kmer size
    ------------
    output
    print the similarity between all the sequences
    '''

    files = load_directory(data_directory)    
    filenames = list(files.keys())

    for filename in filenames :
        with alive_bar(len(files[filename][0])-k+1, ctrl_c=True, title=f"{filename}\t") as bar: # create a progress bar
            files[filename] = stream_kmers(files[filename][0], k, bar)

    for i in range(len(files)):
        for j in range(i+1, len(files)):
            name1, name2 = '_'.join(filenames[i].split('_')[:2]), '_'.join(filenames[j].split('_')[:2])
            with alive_bar(len(files[filenames[i]]) + len(files[filenames[j]]), ctrl_c=True, title=f"{name1} x {name2}\t") as bar:
                A, inter, B = files[filenames[i]], intersection(files[filenames[i]], files[filenames[j]], bar), files[filenames[j]]
            
            print('\t' + name1 + ' : ' + str(round((similarity(len(A), inter, len(B))[0])*100,5)) + ' % covered\n\t' + name2 + ' : ' +  str(round((similarity(len(A), inter, len(B))[1])*100,5)) + ' % covered')
            print('\tJaccard index : ' + str(round(jaccard(len(A), inter, len(B)),5)))

if __name__ == "__main__":
    data_directory = args.data_directory
    k = args.k
    main( data_directory, k)