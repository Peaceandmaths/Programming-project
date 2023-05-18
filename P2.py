"""
PAD Project 2022: P2
Ekaterina Golubeva

Input: list(tuple(string,string))
Output: dict(tuple(int,int) -> tuple(string,string))
Outline :
    - Exception_ ACGT - input check helper function
    - Exception_Format - input check helper function
    - Pairs_Combinations - helper function creating keys for the output dictionary
    - Forward_Phase, previous score - helper function calculating traceback matrix using Needleman-Wunsch algorithm
    - Backward_Phase - helper function creating alignments of sequences based on the traceback matrix
    - AlignByDP - main function

Inspired by :
https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
https://static-content.springer.com/esm/art%3A10.1038%2Fnbt0704-909/MediaObjects/41587_2004_BFnbt0704909_MOESM1_ESM.c
"""
import numpy as np
from itertools import combinations

# Scoring system
match = 5
mismatch = -2
indel = -6

def Exception_ACTG(string):
    """Checks if the nucleotides string only contains allowed characters
    input : string, sequence of nucleotides
    output: True if there's malformed input, False if everything is ok"""
    counter_errors = 0
    for letter in string:
        if letter not in ['A', 'a', 'c', 'C', 'T', 't', 'g', 'G', ' '] or not type(letter) == str:
            counter_errors += 1
            return True
        else:
            counter_errors = 0
    if counter_errors == 0:
        return False


def Exception_Format(list_tuples):
    """Checks if the input is of correct type
     input : list(tuple(string,string))
     output: True if there's malformed input, False if everything is ok"""
    Exception = False
    if type(list_tuples) != list:
        Exception = True
    else:
        for i in range(len(list_tuples)):
            if type(list_tuples[i]) != tuple:
                Exception = True
            else :
                for j in range(len(list_tuples[i])):
                    if type(list_tuples[i][j]) != str:
                        Exception = True
                    else :
                        Exception_ACTG(list_tuples[i][j])
    return Exception


def Pairs_Combinations(list):
    """ Creates list of all combinations of sequence indices
    input : list of tuples with labels and sequences
    output : list of tuples with combinations of alignments
    Inspired by https://www.geeksforgeeks.org/python-itertools-combinations-function/
    """
    output = []
    n = len(list)
    for i in combinations(range(n), 2):
        output.append(i)
    return output


def Forward_Phase(seq1, seq2):
    """ Fills up the matrix according to the algorithm
    input : seq1,seq2, a pair of sequences to align
    output : traceback matrix with the direction code
    diagonal = 1, upper = 2, left = 3
    """
    M, N = len(seq1), len(seq2)
    F = np.zeros((M + 1, N + 1))
    traceback_path = np.zeros((M + 1, N + 1))
    F[(0, 0)] = 0

    def previous_score(i, j):
        """ Calculates the score according to the scoring system :
        input: i, j coordinates of the cell in the matrix
        output: score in int type
        """
        if seq1[i-1] == seq2[j-1]:
            return match
        else:
            return mismatch

    for i in range(1, M+1):
        for j in range(1, N+1):
            F[(i, 0)] = i * indel
            F[(0, j)] = j * indel

            F1 = F[(i-1, j-1)] + previous_score(i, j)
            F2 = F[(i-1, j)] + indel
            F3 = F[(i, j-1)] + indel

            F[(i, j)] = max(F1, F2, F3)

            # Store the directions in the traceback path matrix
            if F[(i, j)] == F1:
                traceback_path[i, j] = 1  # diagonal score
            elif F[(i, j)] == F2:
                traceback_path[i, j] = 2  # upper score
            else:
                traceback_path[i, j] = 3  # left score
    return traceback_path


def Backward_Phase(traceback_path, seq1, seq2):
    """ Traceback using traceback_path, creates alignments
    input : Traceback matrix, seq1,seq2
    output : align1, align2
    """
    align1, align2 = list(), list()
    N, M = len(seq1), len(seq2)
    # We start from right-bottom cell and move up-left
    # through the traceback path till we hit hte first cell
    while N != 0 or M != 0:
        if traceback_path[N, M] == 1:
            # A diagonal direction represents a match,
            # so the letters of the column and the row align
            align1.append(seq1[N-1])
            align2.append(seq2[M-1])
            N -= 1
            M -= 1
        elif traceback_path[N,M] == 2:
            # Up direction represents a mismatch,
            # so a gap ("-") is aligned to the letter of the row
            align1.append(seq1[N-1])
            align2.append('-')
            N -= 1
            M = M
        else:
            # Left direction represents a mismatch,
            # so a gap ("-") is aligned to the letter of the column
            align1.append('-')
            align2.append(seq2[M-1])
            N = N
            M -= 1
    # The nucleotides were aligned in the backward order, so we have to reverse it
    align1.reverse()
    align2.reverse()
    align1 = ''.join(align1)
    align2 = ''.join(align2)
    return align1, align2


def AlignByDP(list):
    """ Takes a list of tuples and saves keys and aligned sequences in a dictionary
    input : list(tuple(string,string))
    output : dictionary {(i,j):(seq1_aligned,seq2_aligned)}"""
    if Exception_Format(list) == True:
        raise Exception("malformed input")
    aligned_seqs = {}
    for j in Pairs_Combinations(list):
        keys = j
        key1 = j[0]
        key2 = j[1]
        seq1 = list[key1][1]
        seq2 = list[key2][1]
        values = Backward_Phase(Forward_Phase(seq1, seq2), seq1, seq2)
        aligned_seqs[keys] = values
    return aligned_seqs
