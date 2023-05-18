"""
PAD Project 2022: P3
Ekaterina Golubeva
Input: dict(tuple(int,int) -> tuple(string,string))
Output: list(list(float)
Outline :
    - Exception_ ACGT - input check helper function
    - Exception_Format - input check helper function
    - Extract_Pairs_Seqs - helper function extracts pairs of sequences from the input dictionary
    - Length - calculates the length of the alignment (without dashes)
    - P_Dist - calculates p-distance between two sequences
    - D_Dist - calculates d-distance between two sequences
    - ComputeDistMatrix - main function
"""
import numpy as np
def Exception_ACTG_(string):
    """Checks if the nucleotides string only contains allowed characters
    input : string, sequence of nucleotides with gaps
    output: True if there's malformed input, False if everything is ok"""
    counter_errors = 0
    for letter in string:
        if letter not in ['A', 'a', 'c', 'C', 'T', 't', 'g', 'G', '-'] or not type(letter) == str:
            counter_errors += 1
            return True
        else:
            counter_errors = 0
    if counter_errors == 0:
        return False


def Exception_Format(dict_tuples):
    """Checks if the input is of correct type and format
     input : dict{tuple(int,int):tuple(string,string)}
     output: True if there's malformed input, False if everything is ok"""
    Exception = False
    if type(dict_tuples) != dict:
        Exception = True
        print(Exception)
    else:
        for key in dict_tuples.keys():

            if type(key) != tuple or type(dict_tuples[key]) != tuple:
                Exception = True
            elif type(key[0]) != int or type(key[1]) != int:
                Exception = True
            elif type(dict_tuples[key][0])!= str or type(dict_tuples[key][1]) != str:
                Exception = True
            else:
                Exception_ACTG_(dict_tuples[key][0])
                Exception_ACTG_(dict_tuples[key][1])
    return Exception


def Extract_Pairs_Seqs(dict, i, j):
    """ Extracts pairs of sequences from the input dictionary
    input : dict, int, int
    output : the pair of sequences as a list of tuples
    """
    pairs = dict[(i, j)]
    seq1 = pairs[0]
    seq2 = pairs[1]
    return seq1, seq2


def Length(seq1, seq2):
    """ Calculates the length of the sequences alignment without dashes
     input : seq1,seq2
     output: length (in int type) of the alignment
    """
    l = len(seq1)
    for i in range(len(seq2)):
        if seq1[i] == '-' or seq2[i] == '-':
            l = l-1
    return l


def P_Dist(l, seq1, seq2):
    """Calculates the p distance between two sequences and stores it in a matrix
    input: length, seq,seq2
    output: p-distance between sequences
    """
    difference = 0
    for i in range(len(seq2)):
        if seq1[i] != seq2[i] and seq1[i] != '-' and seq2[i] != '-':
            difference += 1
    return difference/l


def D_Dist(p):
    """Calculates evolutionary distance between two sequences and stores in a matrix
    input : p-distance
    output: d-distance
    """
    if p >= 3/4:
        d = 30
    else:
        d = -3/4*np.log(1-4/3*p)
    return d


def ComputeDistMatrix(dict):
    """ Calculates evolutionary distance of aligned sequences.
    Input: dict(tuple(int,int) -> tuple(string,string))
    Output: list(list(float)"""
    if Exception_Format(dict):
        raise Exception("malformed input")
    l = len(dict)
    n = round(1/2 + ((1+8*l)**0.5)/2) # number of sequences, matrix dimension
    matrix = np.zeros((n, n))
    for i in range(0, n):
        for j in range(i+1, n):
            seq1, seq2 = Extract_Pairs_Seqs(dict, i, j)
            l = Length(seq1, seq2)
            p = P_Dist(l, seq1, seq2)
            d = D_Dist(p)
            if i == j:
                matrix[i, j] = 0
            else:
                matrix[i, j] = d
                matrix[j, i] = d
    return matrix

