#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################

import IOOps as IOOps
import CheckingOps as CheckOps
import PylogenyOps as PylOps
import os
from Bio.Phylo.Applications import RaxmlCommandline
from pylogeny.executable import consel

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
                 consel,
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


    # 2. find optimal tree with RAxML or FastTree
    PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file=output + '_withoutConstraints',
                   parsimony_seed=10000).runFunction()

    # raxml_cline = RaxmlCommandline(sequences=alignment,
    #                                model=model,
    #                                name=output + "_withoutConstraints")
    # raxml_cline()
    #
    # 3. build trees from contraints and search most likely one
    PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file=output + "_withConstraints",
                   grouping_constraint=constraint,
                   parsimony_seed=10000).runFunction()
    # raxml_cline = RaxmlCommandline(sequences=alignment,
    #                                model=model,
    #                                grouping_constraint=constraint,
    #                                name=output + "_withConstraints")
    # raxml_cline()
    #
    IOOps.Outp().concatTrees(output)


    # # 4. is best tree calculated from contraints (3) significant worse than
    # #    the optimal one (2); if yes return (2); otherwise return (3)
    #
    # # 4.1 build site likelihoods
    PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   alg='g',
                   startingTree='RAxML_bestTree.' + output + '_withoutConstraints',
                   bipartition_filename=output + '_multipleTrees.txt',
                   out_file=output+ ".trees.sitelh",
                   parsimony_seed=10000).runFunction()
    # raxml_cline = RaxmlCommandline(sequences=alignment,
    #                                model=model,
    #                                algorithm='g',
    #                                starting_tree='RAxML_bestTree.' + output + '_withoutConstraints',
    #                                bipartition_filename=output + '_multipleTrees.txt',
    #                                name=output + ".sitelh")
    # raxml_cline()
    #
    # # 4.2 running several test using CONSEL
    consel('RAxML_perSiteLLs.' + output, "alignment", output + "_consel").getInterval()
    # os.system(consel + '/makermt --puzzle RAxML_perSiteLLs.' + output + '.sitelh')
    # os.system(consel + '/consel RAxML_perSiteLLs.' + output)
    # os.system(consel + '/catpv RAxML_perSiteLLs' + ' > consel_output.txt')
