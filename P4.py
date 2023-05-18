"""
PAD Project 2022: P4
Ekaterina Golubeva
Input: list(list(float)
Output: string
Outline :
    - Exception_Format - input check helper function
    - Dict_Dist_Matrix - helper function transfers the input matrix into a 2d dictionary
    - Find_Minimum_Indices - helper function finds minimum value in the dictionary
    - New_Dict - calculates new keys and values for new entries based on WGPMA algorithm
    - Delete_Old_Dict - reduces the dictionary
    - Cluster - main function

Inspired by : https://en.wikipedia.org/wiki/WPGMA
"""

from collections import defaultdict

def Exception_Format(matrix, labels):
    """ Checks that the input has the correct type
    input : list(list(float)), matrix with distance values
    output : True if there's malformed input, False if everything is ok
    """
    Exception = False
    if type(matrix) != list or type(labels) != list:
        Exception = True
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if type(matrix[i][j]) != float:
                Exception = True
            if type(labels[i]) != str:
                print(type(matrix[i][j]), type(labels[i]))
                Exception = True
    return Exception


def Dict_Dist_Matrix(matrix, labels):
    """ Represents the distance matrix as a 2D dictionary
    input : matrix of distances between sequences and list of labels
    output : 2d dictionary with distance values and labels as indices,
    representing the upper half of the matrix
    """
    d = defaultdict(dict)
    for i_i, i in enumerate(labels):
        for i_j, j in enumerate(labels):
            d[i][j] = matrix[i_i][i_j]

    return d


def Find_Minimum_Indices(dict_dist):
    """  Finds minimum pairwise value in a dictionary
    input : 2d dictionary that represents the upper half of the distance matrix
    output :  corresponding indices of the minimum pairwise value
   """
    # Initialize minimum value and corresponding indices arbitrarily
    min = float('inf')
    for i in dict_dist.keys():
        for j in dict_dist[i]:
            if i != j and dict_dist[i][j] < min:
                min = dict_dist[i][j]
                min_i = i
                min_j = j

    return min_i, min_j


def New_Dict(dict_dist, min_i, min_j):
    """ Recalculates the dictionary using WPGMA algorithm: replaces new row and column with average distances,
    input : min_i,min_j indices of the minimum value in the 2d distance dictionary
    output: newly recalculated dictionary
    """
    new_dict = defaultdict(dict)

    for i in dict_dist:
        if i == min_i:
            for k in dict_dist[i]:
                if min_i != k and min_j != k:
                    new_cluster_value = (dict_dist[min_i][k] + dict_dist[min_j][k]) / 2
                    new_dict[k] = new_cluster_value

    new_indices = "(" + min_i + "," + min_j + ")"
    dict_dist[new_indices] = new_dict
    new = "(" + min_i + "," + min_j + ")"
    new_cluster = '('.join(new)
    new_cluster = ')'.join(new_cluster)
    return dict_dist


def Delete_Old_Dict(dict_dist, min_i, min_j):
    """ Reduces the dictionary after joining clusters
    input : dict_dist, min_i, min_j
    output: dict_dist with reduced row, column corresponding to the minimum value
    """

    for i in dict_dist:
        if min_i in dict_dist[i]:
            del dict_dist[i][min_i]
    for j in dict_dist:
        if min_j in dict_dist[j]:
            del dict_dist[j][min_j]

    return dict_dist


def Cluster(matrix, labels):
    if Exception_Format(matrix, labels) == True:
       raise Exception("malformed input")
    else:
    dict_dist = Dict_Dist_Matrix(matrix, labels)
    for _ in range(len(matrix) - 1):
        min_i, min_j = Find_Minimum_Indices(dict_dist)
        new_dict = New_Dict(dict_dist, min_i, min_j)
        dict_dist_new = Delete_Old_Dict(new_dict, min_i, min_j)
        binary_tree = str(dict_dist_new.keys()).replace("dict_keys([", "").replace("])", "")
    return binary_tree

