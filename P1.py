"""
PAD Project 2022: P1
Ekaterina Golubeva

Input: string
Output: list(tuple(string,string))
Outline :
    - Exception_ ACGT - input check helper function
    - Exception_format - input check helper function
    - ParseSeqFile - main function
"""

def Exception_ACTG(string):
    """Checks if the nucleotides string only contains allowed characters
    input : string, sequence of nucleotides
    output: True if there's malformed input, False if everything is ok"""
    counter_errors = 0
    for letter in string:
        if letter not in ['A','a','c', 'C', 'T', 't','g', 'G', ' '] or not type(letter)==str:
            counter_errors += 1
            return True
        else:
            counter_errors = 0
    if counter_errors == 0:
        return False


def Exception_Format(string):
    """Checks if the string is of type string, starts with > followed by a label and data
     input : string : text file with sequences
     output: True if there's malformed input, False if everything is ok"""
    if type(string[1:])!= str:
        return True
    elif string[0] != '>':
        return True
    elif string[0] == '>' and string[1:] =="":
        return True
    elif string[0] == '>' and len(string[1:].split()) == 1:
        return True
    else:
        return False


def ParseSeqFile(string):
    """ Opens a file, reads it line by line and parses each line into a list of tuples
    Input: string, path for the text file with sequences
    Output: list(tuple(string,string))
    """
    strings = open(string, 'r')
    lines = strings.readlines()
    output = []
    for line in lines:
        if line.strip() != "":  # skip empty line
            if Exception_Format(line) == True:
                raise Exception("malformed input")
            else:
                strings_splitted = line[1:].split()
                Label = strings_splitted[0]
                sequence = ''
                for i in range(1, len(strings_splitted)):
                    sequence += strings_splitted[i]
                    if Exception_ACTG(sequence) == True:
                        raise Exception("malformed input")
                output.append((Label, sequence))
    return output















