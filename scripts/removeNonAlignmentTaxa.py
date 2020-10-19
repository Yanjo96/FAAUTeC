#/bin/python3

__author__ = 'Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2020 Project Gruenstaeudl'
__info__ = 'Remove Taxa which are not part of the Contraints'
__version__ = '2020.10.18'

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Nexus import Nexus
from ete3 import Tree

## Parser
def parser():
    '''Simple argparse for an easier usage'''
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser._action_groups.pop()

    ### Arguments ###
    parser.add_argument('-a',
                        '--alignment',
                        help='absolute path to the folder containing the alignmentfiles; infile in PHYLIP or FASTA format; Example: /path_to_input/test.phy',
                        required=True)
    parser.add_argument('-c',
                        '--constraint',
                        help='absolute path to constraint file',
                        required=True)


    args = parser.parse_args()

    main(args.alignment,args.constraint)

def removableTaxa(alignment_path, allTaxa, format):
    Taxa = allTaxa.copy()
    if(format in ["fasta","nex","phy"]):
        for seq_record in SeqIO.parse(alignment_path, "nexus"):
            try:
                Taxa.remove(seq_record.id)
            except:
                pass
        for seq in Taxa:
            try:
                allTaxa.remove(seq)
            except:
                pass
    else:
        print("'" + alignment_path + "' is not in a supported format")

    return allTaxa

def removeTaxa(constraint_path, allTaxa, tree):
    tree = Tree(tree)
    tree.prune(allTaxa)
    return(tree.write(format=9))

## Main Function
def main(alignment_path, constraint_path):
    new_constraint_path = open(''.join(constraint_path.split(".")[:-1]) + "_new." + constraint_path.split(".")[-1], "w")
    with open(constraint_path, "r") as const:
        i = 0
        ## For each tree in tree file find the Taxa which are not part of the alignment
        for tree in const.readlines():
            print(i)
            allTaxa = tree.strip().replace("(","").replace(")","").replace(";","").split(",")
            ## Find the Taxa which are part of all alignmentfiles
            for ali in os.listdir(alignment_path):
                ali = alignment_path + "/" + ali
                allTaxa = removableTaxa(ali, allTaxa, ali.split(".")[-1])
            ## Write new constraint tree to file
            i = i + 1
            new_constraint_path.write(removeTaxa(constraint_path, allTaxa, tree) + "\n")
    new_constraint_path.close()


parser()