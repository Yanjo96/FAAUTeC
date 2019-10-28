#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################

import IOOps as IOOps
import os

###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
__copyright__ = ''
__info__ = 'TreeTopology'
__version__ = '0.1'


#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()


#############
# FUNCTIONS #
#############


def treetopology(alignment,
                 constraint,
                 raxml,
                 model,
                 output):

    # 1. input handler
    #    Maybe check here if alignment file is in PHYLIB or FASTA format
    #    if alignment file is in NEXUS format it have to convert to phylib format
    #    because RAxML can not handle NEXUS files
    if (alignment.split(".")[-1] == 'nex'):
        alignment = IOOps.Inp().nexus2phylib(alignment)

    # 2. find optimal tree with RAxML or FastTree
    os.system(raxml + " --msa " + alignment + " --model " + model + " --threads 1")

    # 3. build trees from contraints and search most likely one

    # 4. is best tree calculated from contraints (3) significant worse than
    #    the optimal one (2); if yes return (2); otherwise return (3)
