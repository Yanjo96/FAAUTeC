#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import os
import datetime

from csv import DictReader
from Bio.Nexus import Nexus
from Bio import SeqIO

try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
__copyright__ = ''
__info__ = 'PlastomeGeneCongruenceTests'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

###########
# CLASSES #
###########

class Inp:
    ''' This class contains functions to conduct miscellaneous input
        operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def nexus2phylip(self, path_to_nex):
        ''' This function convert a NEXUS file to PHYLIP file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            matrix = aln.matrix
        except Nexus.NexusError as e:
            print(e)
            raise
        except Exception as e:
            msg = 'ERROR: %s:\n %s' % ('Parsing of '
            '.nex-file unsuccessful', e)
            warnings.warn(msg)
            raise Exception
        with open(path_to_nex + ".phy", 'w') as phylip:
            phylip.write(str(aln.ntax) + " " + str(aln.nchar) + "\n")
            print(str(aln.ntax) + " " + str(aln.nchar) + "\n")
            for m in matrix:
                phylip.write(str(m)[:10] + (9-len(str(m)[:9])) * " " + " ")
                size = 10
                parts = [str(matrix[m])[i:i+size] for i in range(0, len(str(matrix[m])), size)]
                phylip.write(' '.join(parts) + '\n')
        return path_to_nex + ".phy"

    def nexus2fasta(self, path_to_file):
        ''' This function convert a NEXUS file to FASTA file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_file)
            matrix = aln.matrix
            path_to_file = path_to_file.replace(".nex",".fasta")
        except Nexus.NexusError as e:
            print(e)
            raise
        except Exception as e:
            msg = 'ERROR: %s:\n %s' % ('Parsing of '
            '.nex-file unsuccessful', e)
            warnings.warn(msg)
            raise Exception
        with open(path_to_file, 'w') as fasta:
            for m in matrix:
                fasta.write('>' + str(m) + '\n')
                fasta.write(str(matrix[m]) + '\n')
        return path_to_file


class Outp:
    ''' This class contains two functions for various output operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def concatTrees(self, name):
        with open(name + "_multipleTrees.txt", "w") as multiTree:
            with open("RAxML_bestTree." + name + "_withoutConstraints", "r") as tree:
                for line in tree:
                    multiTree.write(line)
            with open("RAxML_bestTree." + name + "_withConstraints", "r") as tree:
                for line in tree:
                    multiTree.write(line)
