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


    def checkFORMAT(self, alignment, format):
        try:
            return any(SeqIO.parse(alignment, format))
        except:
            return False
