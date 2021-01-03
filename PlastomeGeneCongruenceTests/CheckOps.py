#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import os
from Bio import SeqIO
from Bio.Nexus import Nexus

###############
# AUTHOR INFO #
###############

__author__ = 'Yannick Hartmaring'
__copyright__ = ''
__info__ = 'PlastomeGeneCongruenceTests'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

#############
# Functions #
#############


def checkAlignmentFile(alignment_path):
    nonEmptySeqs = []
    for seq_record in SeqIO.parse(alignment_path, "fasta"):
        if str(seq_record.seq).replace("-","") != "":
            nonEmptySeqs.append(seq_record)

    SeqIO.write(nonEmptySeqs, alignment_path, "fasta")

def checkPrerequisites():
    # Check RAxML

    # Check CONSEL

    # Check Biopython

    # Check dendropy

    # Check mmv
    return True
