#!/usr/bin/env python
'''
Custom operations input and output processes
'''

#####################
# IMPORT OPERATIONS #
#####################

import os
import datetime

from csv import DictReader
from Bio.Nexus import Nexus
from Bio import SeqIO

try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2019 Michael Gruenstaeudl'
__info__ = 'annonex2embl'
__version__ = '2019.10.11.1900'

#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()

###########
# CLASSES #
###########

class Inp:
    ''' This class contains functions to conduct miscellaneous input
        operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass

    def nexus2phylib(self, path_to_nex):
        ''' This function convert a NEXUS file to PHYLIB file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_nex)
            matrix = aln.matrix
        except Nexus.NexusError as e:
            print(e)
            raise
        except Exception as e:
            msg = 'ERROR: %s:\n %s' % ('Parsing of '
            '.nex-file unsuccessful', e)
            warnings.warn(msg)
            raise Exception
        with open(path_to_nex + ".phy", 'w') as phylib:
            phylib.write(str(aln.ntax) + " " + str(aln.nchar) + "\n")
            for m in matrix:
                phylib.write(str(m)[:10] + (9-len(str(m)[:9])) * " " + " ")
                size = 10
                parts = [str(matrix[m])[i:i+size] for i in range(0, len(str(matrix[m])), size)]
                phylib.write(' '.join(parts) + '\n')
        return path_to_nex + ".phy"


class Outp:
    ''' This class contains two functions for various output operations.
    Args:
        [specific to function]
    Returns:
        [specific to function]
    Raises:
        -
    '''

    def __init__(self):
        pass
