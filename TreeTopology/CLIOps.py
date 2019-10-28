#!/usr/bin/env python2.7
'''
Command-line execution of annonex2embl
'''

#####################
# IMPORT OPERATIONS #
#####################

import sys
import os

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'annonex2embl'))

# IMPORTANT: TFL must be after "sys.path.append"
import TreeTopologyMain
import argparse

###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2019 Project Gruenstaeudl'
__info__ = 'TreeTopology'
__version__ = '2019.10.25'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

############
# ARGPARSE #
############

class CLI():

    def __init__(self):
        self.client()

    def client(self):

        parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
        parser._action_groups.pop()
        required = parser.add_argument_group('required arguments')
        optional = parser.add_argument_group('optional arguments')


        ### REQUIRED ###
        required.add_argument('-a',
                            '--alignment',
                            help='absolute path to infile; infile in NEXUS, PHYLIB or FASTA format; Example: /path_to_input/test.nex',
                            required=True)

        required.add_argument('-c',
                            '--constraint',
                            help='absolute path to constraint file; infile in NEWICK format; Example: /path_to_input/tree.txt',
                            required=True)

        required.add_argument('-r',
                            '--raxml',
                            help='path to RAxML',
                            required=True)

        required.add_argument('-o',
                            '--output',
                            help='Output Tree in Newick format',
                            required=True)

        ### OPTIONAL ###

        optional.add_argument('--model',
                            help='Model for RAxML',
                            default="GTR+G",
                            required=False)

        optional.add_argument('--version',
                            help='print version number and exit',
                            action='version',
                            version='%(prog)s ' + __version__)

        args = parser.parse_args()

        TreeTopologyMain.treetopology(args.alignment,
                                      args.constraint,
                                      args.raxml,
                                      args.model,
                                      args.output)

########
# MAIN #
########

def start_treetopology():
    CLI()
