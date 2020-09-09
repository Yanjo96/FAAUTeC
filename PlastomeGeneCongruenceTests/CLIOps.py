#!/usr/bin/env python3
'''
Command-line execution
'''

#####################
# IMPORT OPERATIONS #
#####################

import sys
import os

# Add specific directory to sys.path in order to import its modules
# NOTE: THIS RELATIVE IMPORTING IS AMATEURISH.
# NOTE: COULD THE FOLLOWING IMPORT BE REPLACED WITH 'import annonex2embl'?
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plastomeGeneCongruenceTests'))

# IMPORTANT: TFL must be after "sys.path.append"
import PlastomeGeneCongruenceTestsMain
import argparse

###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2019 Project Gruenstaeudl'
__info__ = 'PlastomeGeneCongruenceTests'
__version__ = '2019.12.03'

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
                            help='absolute path to infile; infile in PHYLIP or FASTA format; Example: /path_to_input/test.phy',
                            required=True)

        required.add_argument('-c',
                            '--constraint',
                            help='absolute path to constraint file; infile in NEWICK format; Example: /path_to_input/tree.tre',
                            required=True)

        required.add_argument('--consel',
                              help='path to consel executables',
                              required=True)

        ### OPTIONAL ###

        optional.add_argument('--model',
                            help='Model for RAxML',
                            default="GTRGAMMAI",
                            required=False)

        optional.add_argument('--iqtree2',
                            help='absolute path to the iqtree2 executable',
                            default=False,
                            required=False)

        optional.add_argument('--latex',
                             help='',
                             default=False,
                             action='store_true',
                             required=False)

        optional.add_argument('--version',
                            help='print version number and exit',
                            action='version',
                            version='%(prog)s ' + __version__)

        args = parser.parse_args()

        PlastomeGeneCongruenceTestsMain.plastomeGeneCongruenceTests(args.alignment,
                                                                    args.constraint,
                                                                    args.consel,
                                                                    args.model,
                                                                    args.iqtree2,
                                                                    args.latex)

########
# MAIN #
########

def start_plastomeGeneCongruenceTests():
    CLI()
