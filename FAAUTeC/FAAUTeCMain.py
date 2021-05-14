#!/usr/bin/env python3.7

''' Main operations in FAAUTeC '''

#####################
# IMPORT OPERATIONS #
#####################


import os
import time
import dendropy

import IOOps as IOOps
import FAAUTeCOps as FOps
import CheckOps as COps


###############
# AUTHOR INFO #
###############

__author__ = 'Yannick Hartmaring'
__ctb__ = 'Nils Jenke'
__copyright__ = ''
__info__ = 'FAAUTeC'
__version__ = '0.2'


#############
# DEBUGGING #
#############

#import ipdb
#ipdb.set_trace()


#############
# FUNCTIONS #
#############


def faautec(alignment,
            constraint_path,
            consel_path,
            model,
            mlcalc,
            threadNumber,
            iqtree2_path,
            iqtree_path,
            raxml_path,
            latex):

    os.system("mkdir output/")
    os.system("mkdir output/SUMMARY")

    constraints = FOps.readConstraints(constraint_path)
    place = ["hypo" + str(i) for i in range(len(constraints))]
    programs = ["CONSEL","IQTree"]
    if(iqtree2_path):
        programs.append("IQTree2")

    treeFile_raxml = open("output/SUMMARY/raxml_hypoTreeShortestDistUnconstTree.tre","w")
    treeFile_iqtree = open("output/SUMMARY/iqree_hypoTreeShortestDistUnconstTree.tre","w")
    auFile = open("output/SUMMARY/au_runtime_table.csv","w")

    auFile.write("gene," + ','.join([','.join([program + "_hypo" + str(i) for program in programs]) for i in range(len(constraints))]) + "," + ','.join(["runtime_" + program for program in programs]) + "\n")

    alis = os.listdir(alignment)
    colored = 0
    for ali in alis:
        print(ali)
        gene = ali.split(".")[0]
        ali = alignment + ali.strip()
        if ali.split(".")[-1] != "fasta":
            if ali.split(".")[-1] == "phy":
                ali = IOOps.Inp().phylip2fasta(ali)
            elif ali.split(".")[-1] == "nex":
                ali = IOOps.Inp().nexus2fasta(ali)
            else:
                print(ali + " was skipped because the file ending is not sopprted. Supported File endings: 'fasta', 'nex', 'phy'")
                continue

        COps.checkAlignmentFile(ali)

        log = ["#!/bin/bash",
               "# " + gene]

        ### Create a clear file system
        os.mkdir("output/" + gene)
        os.mkdir("output/" + gene + "/01_input")
        os.mkdir("output/" + gene + "/02_output_" + mlcalc)
        os.mkdir("output/" + gene + "/03a_output_CONSEL")
        os.mkdir("output/" + gene + "/03b_output_IQTree")

        if(iqtree2_path):
            log.append(FOps.commandline("mkdir output/" + gene + "/03c_output_IQTree2"))

        log.append(FOps.commandline("cp " + ali + " output/" + gene + "/01_input/"))

        if(mlcalc == "RAxML"):
            log = log + ["\n",
                        "# Calculate ML-Trees with RAxML"]

            ### Calculate the ML-Trees with RAxML
            log = log + FOps.raxml(ali, constraints, model, gene, threadNumber, raxml_path)

            ### Find Tree which has the smallest euclidic distance to
            ### to the unconstraint tree
            trees = []
            with open(gene + '_COMBINED.tre', "r") as multitree:
                for tree in multitree:
                    trees.append(tree.strip())
            best_tree = FOps.findBestTree(trees)
            trees = trees[1:]

            log.append(FOps.commandline("mv RAxML_* output/" + gene + "/02_output_RAxML/"))

        elif(mlcalc == "IQTree"):
            log = log + ["\n",
                        "# Calculate ML-Trees with IQTree"]

            ### Calculate the ML-Trees with IQTree
            log = log + FOps.iqtree_mltree(ali, constraints, gene, threadNumber,iqtree_path)

            ### Find Tree which has the smallest euclidic distance to
            ### to the unconstraint tree
            trees = []
            with open(gene + '_COMBINED.tre', "r") as multitree:
                for tree in multitree:
                    trees.append(tree.strip())
            best_tree = FOps.findBestTree(trees)
            trees = trees[1:]

            log.append(FOps.commandline("mv " + gene + "_IQTree* output/" + gene + "/02_output_IQTree/"))

        else:
            print("Error: The Program " + mlcalc + " is not supported to run ML Tree calculation use" +
                  " 'RAxML' or 'IQTree' instead.")
            break

        ## write the best tree to the best tree file
        treeFile_raxml.write(gene + place[best_tree] + " " + trees[best_tree])

        ### AU Test by CONSEL
        log = log + ["\n",
                     "# Calculate AU-Test with CONSEL"]

        start = time.time()
        log = log + FOps.consel(ali, consel_path, model, gene, mlcalc, threadNumber, raxml_path)
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

        log.append(FOps.commandline("mv " + gene + "_CONSEL* output/" + gene + "/03a_output_CONSEL/"))
        log.append(FOps.commandline("mv RAxML* output/" + gene + "/03a_output_CONSEL/"))

        # AU Test by IQTree
        print("AU Test by IQTree")
        log = log + ["\n",
                     "# Calculate AU-Test with IQTree"]
        start = time.time()
        log = log + FOps.iqtree_autest(ali, gene, mlcalc, threadNumber, iqtree_path)
        iqtree_runtime = round(time.time() - start,3)


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

        ## mark values below significance level
        for i in range(len(au_iqtree)):
            if au_iqtree[i] <= 0.05:
                 au_iqtree[i] = str(au_iqtree[i]) + "*"

        ## mark the tree with the smallest euclidic distance
        au_iqtree[best_tree] = str(au_iqtree[best_tree]) + "s"

        #llh_iqtree = []
        #for llh_file in [gene + "_IQTree_" + name + ".iqtree" for name in ["unconst"] + ["hypo" + str(i) for i in range(len(constraints))]]:
        #    with open(llh_file,"r") as llh:
        #        for line in llh.readlines():
        #            if "Log-likelihood of the tree:" in line:
        #                llh_iqtree.append(float(line.split()[4]))
        #                break

        log.append(FOps.commandline("mv " + gene + "_IQTree* output/" + gene + "/03b_output_IQTree/"))


        # AU Test by IQTree2
        if(iqtree2_path):
            print("AU Test by IQTree2")
            log = log + ["\n",
                         "# Calculate AU-Test with IQTree"]
            start = time.time()
            log = log + FOps.iqtree2_autest(ali, iqtree2_path, gene, mlcalc, threadNumber)
            iqtree2_runtime = round(time.time() - start,3)

            au_iqtree2 = []
            with open(gene + "_IQTree.iqtree","r") as iqtree2_out:
                for line in iqtree2_out:
                    if line.strip() == "-------------------------------------------------------------------------":
                        break
                    else:
                        pass
                for line in iqtree2_out:
                    try:
                        au_iqtree2.append(float(line.split()[11]))
                    except:
                        break
            au_iqtree2 = au_iqtree2[1:]

            ## mark values below significance level
            for i in range(len(au_iqtree2)):
                if au_iqtree2[i] <= 0.05:
                    au_iqtree2[i] = str(au_iqtree2[i]) + "*"

            ## mark the tree with the smallest euclidic distance
            au_iqtree2[best_tree] = str(au_iqtree2[best_tree]) + "s"

            log.append(FOps.commandline("mv " + gene + "_IQTree* output/" + gene + "/03c_output_IQTree2/"))


        #treeFile_iqtree.write(gene + place[best_iqtree] + " " + trees_iqtree[best_iqtree] + "\n")

        if(iqtree2_path):
            auFile.write(gene + "," + ''.join([str(au_consel[i]) + "," + str(au_iqtree[i]) + "," + str(au_iqtree2[i]) + "," for i in range(len(constraints))]) + str(consel_runtime) + "," + str(iqtree_runtime) + "," + str(iqtree2_runtime) + "\n")
        else:
            auFile.write(gene + "," + ''.join([str(au_consel[i]) + "," + str(au_iqtree[i]) + "," for i in range(len(constraints))]) + str(consel_runtime) + "," + str(iqtree_runtime) + "\n")
        #if colored % 2 == 0:
        #    latexFile.write(" \\rowcolor{black!20} ")
        #    llsFile.write(" \\rowcolor{black!20} ")

        #latexFile.write("\\textit{" + gene + "}" + " & " + str(edit_num(au_pylogeny[0])) + " & " + str(edit_num(au_pylogeny[1])) + " & " + str(edit_num(au_pylogeny[2])) + " & " + str(edit_num(au_iqtree[0])) + " & " + str(edit_num(au_iqtree[1])) + " & " + str(edit_num(au_iqtree[2])) + " & " + str(round(pylogeny_runtime,3)) + " & " + str(round(iqtree_runtime,3)) + "\\\\" + "\n")
        #llsFile.write("\\textit{" + gene + "}" + " & " + str(llh_pylogeny[0]) + " & " + str(llh_pylogeny[1]) + " & " + str(llh_pylogeny[2]) + " & " + str(llh_pylogeny[3]) + " & " + str(llh_iqtree[0]) + " & " + str(llh_iqtree[1]) + " & " + str(llh_iqtree[2]) + " & " + str(llh_iqtree[3]) + "\\\\\n")

        #if colored % 2 == 0:
        #    llsFile.write(" \\rowcolor{black!20} ")
        #llsFile.write("$\\Delta$" + " & " + "" + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[1],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[2],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[3],4)) + " & " + "" + " & " + str(round(llh_iqtree[0] - llh_iqtree[1],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[2],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[3],4)) + "\\\\\n")

        colored = colored + 1

        with open("output/" + gene + "/" + gene + "_log.sh","w") as logFile:
            for line in log:
                logFile.write(line + "\n")

        os.remove(gene + "_COMBINED.tre")
        os.remove("hypo.txt")
        print(str(colored) + " / " + str(int(len(alis))))


    treeFile_raxml.close()
    treeFile_iqtree.close()
    auFile.close()

    if(latex):
        IOOps.Outp().createLatex("output/SUMMARY/au_runtime_table.csv", latex)

    #os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/likelihoods_table.tex")