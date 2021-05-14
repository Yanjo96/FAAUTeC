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

__author__ = 'Yannick Hartmaring'
__ctb__ = 'Nils Jenke'
__copyright__ = ''
__info__ = 'FAAUTeC'

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

    def phylip2fasta(self, path_to_phy):
        ''' This function convert a PHYLIP file to FASTA file. '''
        phylipFile = open(path_to_phy, "r")
        alignment = phylipFile.readlines()
        phylipFile.close()
        path_to_fasta = path_to_phy.replace(".phy",".fasta")

        with open(path_to_fasta, "w") as fasta:
            for ali in alignment[1:]:
                ali = ali.split()
                fasta.write(">" + ali[0] + "\n")
                fasta.write(ali[1].replace("N","-") + "\n")
        return path_to_fasta

    def nexus2fasta(self, path_to_file):
        ''' This function convert a NEXUS file to FASTA file. '''
        try:
            aln = Nexus.Nexus()
            aln.read(path_to_file)
            matrix = aln.matrix
            path_to_file = path_to_file.replace(".nex",".fasta")
        except Nexus.NexusError as e:
            print(e)
            raise
        except Exception as e:
            msg = 'ERROR: %s:\n %s' % ('Parsing of '
            '.nex-file unsuccessful', e)
            warnings.warn(msg)
            raise Exception
        with open(path_to_file, 'w') as fasta:
            for m in matrix:
                fasta.write('>' + str(m) + '\n')
                fasta.write(str(matrix[m]) + '\n')
        return path_to_file


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

    def concatTrees(self, name):
        with open(name + "_multipleTrees.txt", "w") as multiTree:
            with open("RAxML_bestTree." + name + "_withoutConstraints", "r") as tree:
                for line in tree:
                    multiTree.write(line)
            with open("RAxML_bestTree." + name + "_withConstraints", "r") as tree:
                for line in tree:
                    multiTree.write(line)

    def edit_num(self, value):
        try:
            postSigns = ''.join(list(filter(lambda x: x not in "0123456789.-e", value)))
            postSigns = postSigns.replace("s","\\textsuperscript{s}").replace("*","\\textsuperscript{*}")
            preSigns = ""
            value = round(float(value.replace("s","").replace("*","")),3)
            if value < 0.001:
                value = 0.001
                preSigns = preSigns + "<"
            value = preSigns + str(value) + postSigns
            if "*" in str(value):
                return "\\textbf{" + str(value) + "}"
            return value
        except:
            return str(value)


    def createLatex(self, csv, output):
        '''
        This function converts the csv table output to a latex table
        It is optional but creates a more fancy look to present the data
        '''

        latexFile = open(output, "w")
        latexFile.write("\\documentclass[a4paper]{article}\n")
        latexFile.write("\\usepackage[T1]{fontenc}\n")
        latexFile.write("\\usepackage{geometry, longtable}\n")
        latexFile.write("\\usepackage{booktabs}\n")
        latexFile.write("\\usepackage[cmyk,table]{xcolor}\n")
        latexFile.write("\\geometry{paperheight=297mm, paperwidth=210mm, margin=2pt}\n")
        latexFile.write("\\pagenumbering{gobble}\n")
        latexFile.write("\\begin{document}\n")
        latexFile.write("%\\footnotesize\n")
        latexFile.write("%\\rowcolors{1}{white}{black!20}\n")
        with open(csv, "r") as csvFile:
            line = csvFile.readline()
            colNum = len(line.split(","))
            latexFile.write("\\begin{longtable}{|" + ('m{' + str(round(1/colNum,3)) + '\\textwidth}|') * (colNum) + "}\n")
            latexFile.write("\\caption[]{\\textsuperscript{s}tree with lowest distance to unconstraint tree; \\textsuperscript{*}p-value $\\leq$ 0.05\n}\\\\\n")
            latexFile.write("\\toprule\n")
            latexFile.write('&'.join([str(i) for i in line.strip().replace("_"," ").split(",")]) + "\\\\\n")
            latexFile.write("\\midrule\n")
            latexFile.write("\\endhead\n")
            latexFile.write("\\bottomrule\n")
            latexFile.write("\\endfoot\n")
            for line in csvFile.readlines():
                latexFile.write('&'.join([self.edit_num(str(i)) for i in line.strip().replace("_"," ").split(",")]) + "\\\\\n")
        latexFile.write("\\bottomrule\n")
        latexFile.write("\\end{longtable}\n")
        latexFile.write("\\end{document}\n")
        latexFile.close()

    def createLatex2(self, constNumber, programs):
        latexFile = open("output/SUMMARY/au_runtime_table.tex","w")
        latexFile.write("\\documentclass[a4paper]{article}\n")
        latexFile.write("\\usepackage{colortbl, geometry}\n")
        latexFile.write("\\usepackage[cmyk,table]{xcolor}\n")
        latexFile.write("\\geometry{paperheight=297mm, paperwidth=210mm, margin=2pt}\n")
        latexFile.write("\\pagenumbering{gobble}\n")
        latexFile.write("\\begin{document}\n")
        latexFile.write("\\footnotesize\n")
        latexFile.write("\\begin{tabular}{l|" + ''.join(["r" * len(programs) for i in range(constNumber+1)]) + "}\\\\\n")
        latexFile.write("gene " + ''.join([" & \\multicolumn{" + str(len(programs)) + "}{c}{Hypothesis " + str(i) + "}" for i in range(constNumber)]) + " & \\multicolumn{" + str(len(programs)) +"}{c}{Runtime in seconds}\\\\\n")
        latexFile.write("\\hline\\\\\n")
        latexFile.write(' & ' + '&'.join(['&'.join([program for program in programs]) for i in range(constNumber+1)]) + "\\\\\n")

        #llsFile = open("output/SUMMARY/likelihoods_table.tex","w")
        auFile = open("output/SUMMARY/au_runtime_table.csv", "r")
        au = auFile.readlines()
        auFile.close()
        for line in au[1:]:
            latexFile.write(' & '.join(line.replace("_","\_").split(",")) + "\\\\n")

        latexFile.write("\\end{tabular}\\\\\n")
        latexFile.write("\\textsuperscript{s}tree with lowest distance to unconstraint tree; \\textsuperscript{*}p-value $\\leq$ 0.05\n")
        latexFile.write("\\end{document}\n")
        #llsFile.write("\\documentclass[a4paper]{article}\n")
        #llsFile.write("\\usepackage{colortbl, geometry, longtable}\n")
        #llsFile.write("\\usepackage[cmyk,table]{xcolor}\n")
        #llsFile.write("\\geometry{paperheight=297mm, paperwidth=210mm, margin=2pt}\n")
        #llsFile.write("\\pagenumbering{gobble}\n")
        #llsFile.write("\\begin{document}\n")
        #llsFile.write("\\footnotesize\n")
        #llsFile.write("\\begin{longtable}{p{0.03\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}}\n")
        #llsFile.write("gene & \\multicolumn{4}{c}{Pylogeny} & \\multicolumn{4}{c}{IQTree}\\\\\n")
        #llsFile.write(" & Unconstraint & Hypothesis A & Hypothesis B & Hypothesis A & Unconstraint & Hypothesis A & Hypothesis B & Hypothesis C\\\\\n")
        #llsFile.write("\\endhead\n")

        latexFile.close()

        os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/au_runtime_table.tex ")
