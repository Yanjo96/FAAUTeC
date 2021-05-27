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
__info__ = 'FAAUTeC'

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

def checkPrerequisites(au_inference, path_iqtree2, path_consel, ml_inference):
    # Output folder
    if 'output' in os.listdir():
        print("there is already an output folder, please rename or remove it before running FAAUTeC")
        return(False)

    # Check AU Inference
    programs = au_inference.split(";")
    if(len(programs) == 0):
        print("Please specify at least one program for AU Test calculation")
        return(False)

    if("IQTree2" in programs and not path_iqtree2):
        print("Please specify the path of IQTree2 by the parameter --path_iqtree2")
        return(False)
    if("CONSEL" in programs and not path_consel):
        print("Please specify the path of CONSEL by the parameter --path_consel")
        return(False)

    try:
        programs.remove("CONSEL")
    except:
        pass

    try:
        programs.remove("IQTree")
    except:
        pass

    try:
        programs.remove("IQTree2")
    except:
        pass

    if(len(programs) > 0):
        print("'" + ' '.join(programs) + "' is not a supported Program for AU Test calculation, supported programs are: 'CONSEL', 'IQTree' and 'IQTree2'")
        return(False)

    # Check ML Inference
    if(ml_inference != "RAxML" and ml_inference != "IQTree"):
        print("'" + ml_inference + "' is not a supported program for ML inference, supported programs are: 'RAxML' and 'IQTree'")
        return(False)

    return(True)
