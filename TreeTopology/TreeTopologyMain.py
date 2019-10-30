#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################

import IOOps as IOOps
import CheckingOps as CheckOps
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
    if CheckOps.InputFormatCheck().checkFORMAT(alignment, "phylip-relaxed"):
        pass
    elif CheckOps.InputFormatCheck().checkFORMAT(alignment, "fasta"):
        pass
    # if alignment file is in NEXUS format it have to convert to phylip format
    # because RAxML can not handle NEXUS files
    elif CheckOps.InputFormatCheck().checkFORMAT(alignment, "nexus"):
        alignment = IOOps.Inp().nexus2phylip(alignment)
    else:
        raise Exception("File " + alignment + " is not in supported file" +
                        " format (PHYLIP, FASTA, NEXUS)")


    # 2. find optimal tree with RAxML or FastTree (use pylogeny)
    os.system(raxml + " --msa " + alignment + " --model " + model + " --threads 1")

    # 3. build trees from contraints and search most likely one (use pylogeny)

    # 4. is best tree calculated from contraints (3) significant worse than
    #    the optimal one (2); if yes return (2); otherwise return (3)
    #    (use pylogeny)
