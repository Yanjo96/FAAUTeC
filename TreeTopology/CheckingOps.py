#!/usr/bin/env python
'''
Custom operations to check annotations
'''

#####################
# IMPORT OPERATIONS #
#####################
from Bio import SeqIO

###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
__copyright__ = ''

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

###########
# CLASSES #
###########

class InputFormatCheck():

    def __init__(self):
        pass

    def checkPHYLIB(self):
        pass

    def checkFASTA(self, alignment):
        with open(alignment, "r") as handle:
            
