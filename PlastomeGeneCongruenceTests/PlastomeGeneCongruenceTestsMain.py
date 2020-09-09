#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################


import os
import time
import dendropy

import IOOps as IOOps
import PylogenyOps as PylOps
import CheckOps as COps


###############
# AUTHOR INFO #
###############

__author__ = 'Yannick Hartmaring'
__ctb__ = 'Nils Jenke'
__copyright__ = ''
__info__ = 'PlastomeGeneCongruenceTests'
__version__ = '0.2'


#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()


#############
# FUNCTIONS #
#############


def plastomeGeneCongruenceTests(alignment,
                 constraint_path,
                 consel_path,
                 model,
                 latex):

    os.system("mkdir output/")
    os.system("mkdir output/SUMMARY")

    constraints = PylOps.readConstraints(constraint_path)
    place = ["hypo" + str(i) for i in range(len(constraints))]

    treeFile_raxml = open("output/SUMMARY/raxml_hypoTreeShortestDistUnconstTree.tre","w")
    treeFile_iqtree = open("output/SUMMARY/iqree_hypoTreeShortestDistUnconstTree.tre","w")
    auFile = open("output/SUMMARY/au_runtime_table.csv","w")

    auFile.write("gene," + ''.join(["consel_hypo" + str(i) + ",iqtree_hypo" + str(i) + "," for i in range(len(constraints))]) + "runtime_consel,runtime_iqtree\n")

    alis = os.listdir(alignment)
    colored = 0
    for ali in alis:
        gene = ali.split(".")[0]#.split("_")[1]
        ali = alignment + ali.strip()
        if ali.split(".")[-1] != "fasta" and ali.split(".")[-1] != "phy":
            continue

        log = ["#!/bin/bash",
               "# " + gene]

        ### Create a clear file system
        log.append(PylOps.commandline("mkdir output/" + gene))
        log.append(PylOps.commandline("mkdir output/" + gene + "/01_input"))
        log.append(PylOps.commandline("mkdir output/" + gene + "/02_output_RAxML"))
        log.append(PylOps.commandline("mkdir output/" + gene + "/03a_output_CONSEL"))
        log.append(PylOps.commandline("mkdir output/" + gene + "/03b_output_IQTree"))

        log.append(PylOps.commandline("cp " + ali + " output/" + gene + "/01_input/"))

        log = log + ["\n",
                     "# Calculate ML-Trees with RAxML"]

        ### Calculate the ML-Trees with RAxML
        log = log + PylOps.raxml(ali, constraints, model, gene)

        ### Find Tree which has the smallest euclidic distance to
        ### to the unconstraint tree
        trees = []
        with open('RAxML_' + gene + '_COMBINED.tre', "r") as multitree:
            for tree in multitree:
                trees.append(tree.strip())
        best_tree = PylOps.findBestTree(trees)
        trees = trees[1:]

        ## write the best tree to the best tree file
        treeFile_raxml.write(gene + place[best_tree] + " " + trees[best_tree])

        ### AU Test by CONSEL
        log = log + ["\n",
                     "# Calculate AU-Test with CONSEL"]

        start = time.time()
        log = log + PylOps.consel(ali, consel_path, model, gene)
        consel_runtime = round(time.time() - start,3)

        ## Save the AU Test values to a variable
        au_consel = [0] * (len(constraints) + 1)
        with open(gene + "_CONSEL.consel") as consel:
            for line in consel.readlines()[3:len(constraints)+4]:
                au_consel[int(line.split()[2]) - 1] = float(line.split()[4])
        au_consel = au_consel[1:]

        ## mark values below significance level
        for i in range(len(au_consel)):
            if au_consel[i] <= 0.05:
                 au_consel[i] = str(au_consel[i]) + "*"

        ## mark the tree with the smallest euclidic distance
        au_consel[best_tree] = str(au_consel[best_tree]) + "s"

        #llh_raxml = []
        #for llh_file in ["RAxML_log." + name + "_" + gene for name in  ["withoutConstraints"] + ["hypothesis" + str(i) for i in range(len(constraints))]]:
        #    with open(llh_file,"r") as llh:
        #        llh_raxml.append(float(llh.readlines()[-1].split()[-1]))

        log.append(PylOps.commandline("mv " + gene + "_CONSEL* output/" + gene + "/03a_output_CONSEL/"))

        # AU Test by IQTree
        log = log + ["\n",
                     "# Calculate AU-Test with IQTree"]
        start = time.time()
        log = log + PylOps.iqtree_autest(ali, gene)
        iqtree_runtime = round(time.time() - start,3)

        #trees_iqtree = []
        #with open(gene + '_COMBINED.tre', "r") as multitree:
        #    for tree in multitree:
        #        trees_iqtree.append(tree.strip())
        #best_iqtree = PylOps.findBestTree(trees_iqtree)
        #trees_iqtree = trees_iqtree[1:]

        au_iqtree = []
        with open(gene + "_IQTree.iqtree","r") as iqtree_out:
            for line in iqtree_out:
                if line.strip() == "-------------------------------------------------------------------------":
                    break
                else:
                    pass
            for line in iqtree_out:
                try:
                    au_iqtree.append(float(line.split()[11]))
                except:
                    break
        au_iqtree = au_iqtree[1:]

        #llh_iqtree = []
        #for llh_file in [gene + "_IQTree_" + name + ".iqtree" for name in ["unconst"] + ["hypo" + str(i) for i in range(len(constraints))]]:
        #    with open(llh_file,"r") as llh:
        #        for line in llh.readlines():
        #            if "Log-likelihood of the tree:" in line:
        #                llh_iqtree.append(float(line.split()[4]))
        #                break

        log.append(PylOps.commandline("mv " + gene + "_IQTree* output/" + gene + "/03b_output_IQTree/"))

        #treeFile_iqtree.write(gene + place[best_iqtree] + " " + trees_iqtree[best_iqtree] + "\n")

        ## mark values below significance level
        for i in range(len(au_iqtree)):
            if au_iqtree[i] <= 0.05:
                 au_iqtree[i] = str(au_iqtree[i]) + "*"

        ## mark the tree with the smallest euclidic distance
        au_iqtree[best_tree] = str(au_iqtree[best_tree]) + "s"

        auFile.write(gene + "," + ''.join([str(au_consel[i]) + "," + str(au_iqtree[i]) + "," for i in range(len(constraints))]) + str(consel_runtime) + "," + str(iqtree_runtime) + "\n")
        #if colored % 2 == 0:
        #    latexFile.write(" \\rowcolor{black!20} ")
        #    llsFile.write(" \\rowcolor{black!20} ")

        #latexFile.write("\\textit{" + gene + "}" + " & " + str(edit_num(au_pylogeny[0])) + " & " + str(edit_num(au_pylogeny[1])) + " & " + str(edit_num(au_pylogeny[2])) + " & " + str(edit_num(au_iqtree[0])) + " & " + str(edit_num(au_iqtree[1])) + " & " + str(edit_num(au_iqtree[2])) + " & " + str(round(pylogeny_runtime,3)) + " & " + str(round(iqtree_runtime,3)) + "\\\\" + "\n")
        #llsFile.write("\\textit{" + gene + "}" + " & " + str(llh_pylogeny[0]) + " & " + str(llh_pylogeny[1]) + " & " + str(llh_pylogeny[2]) + " & " + str(llh_pylogeny[3]) + " & " + str(llh_iqtree[0]) + " & " + str(llh_iqtree[1]) + " & " + str(llh_iqtree[2]) + " & " + str(llh_iqtree[3]) + "\\\\\n")

        #if colored % 2 == 0:
        #    llsFile.write(" \\rowcolor{black!20} ")
        #llsFile.write("$\\Delta$" + " & " + "" + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[1],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[2],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[3],4)) + " & " + "" + " & " + str(round(llh_iqtree[0] - llh_iqtree[1],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[2],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[3],4)) + "\\\\\n")

        log.append(PylOps.commandline("mv RAxML_* output/" + gene + "/02_output_RAxML/"))
        colored = colored + 1

        with open("output/" + gene + "/" + gene + "_log.sh","w") as logFile:
            for line in log:
                logFile.write(line + "\n")

        #os.remove(gene + "_COMBINED.tre")
        print(str(colored) + " / " + str(int(len(alis)/2)))


    #latexFile.write("\\end{tabular}\\\\\n")
    #latexFile.write("\\textsuperscript{s}tree with lowest distance to unconstraint tree; \\textsuperscript{*}p-value $\\leq$ 0.05\n")
    #latexFile.write("\\end{document}\n")

    #llsFile.write("\\end{longtable}\n")
    #llsFile.write("\\end{document}\n")

    #llsFile.close()
    #latexFile.close()
    treeFile_raxml.close()
    treeFile_iqtree.close()
    auFile.close()

    if(latex):
        PylOps.createLatex(len(constraints))

    #os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/au_runtime_table.tex ")
    #os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/likelihoods_table.tex")
