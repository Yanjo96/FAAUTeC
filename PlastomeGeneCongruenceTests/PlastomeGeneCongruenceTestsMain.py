#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################

import IOOps as IOOps
import os
import time
import dendropy


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
                 model):

    os.system("mkdir output/")
    os.system("mkdir output/SUMMARY")

    constraints = readConstraints(constraint_path)
    place = ["hypo" + str(i) for i in range(len(constraints))]

    treeFile_raxml = open("output/SUMMARY/raxml_hypoTreeShortestDistUnconstTree.tre","w")
    treeFile_iqtree = open("output/SUMMARY/iqree_hypoTreeShortestDistUnconstTree.tre","w")
    auFile = open("output/SUMMARY/au_runtime_table.csv","w")
    latexFile = open("output/SUMMARY/au_runtime_table.tex","w")
    llsFile = open("output/SUMMARY/likelihoods_table.tex","w")

    latexFile.write("\\documentclass[a4paper]{article}\n")
    latexFile.write("\\usepackage{colortbl, geometry}\n")
    latexFile.write("\\usepackage[cmyk,table]{xcolor}\n")
    latexFile.write("\\geometry{paperheight=297mm, paperwidth=210mm, margin=2pt}\n")
    latexFile.write("\\pagenumbering{gobble}\n")
    latexFile.write("\\begin{document}\n")
    latexFile.write("\\footnotesize\n")
    latexFile.write("\\begin{tabular}{p{0.03\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}}\n")
    latexFile.write("gene & \\multicolumn{3}{c}{RAxML}  & \\multicolumn{3}{c}{IQTree} & \\multicolumn{2}{c}{Runtime in seconds}\\\\\n")
    latexFile.write(" & Hypothesis A & Hypothesis B & Hypothesis C & Hypothesis A & Hypothesis B & Hypothesis C & RAxML & IQTree\\\\\n")

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

    auFile.write("gene," + ''.join(["hypo" + str(i) + "_raxml," for i in range(len(constraints))]) + ''.join(["hypo" + str(i) + "_iqtree," for i in range(len(constraints))]) + "runtime_pylogeny,runtime_iqtree\n")

    os.system("ls " + alignment + " > aligments.tmp")
    with open("aligments.tmp","r") as alignments:
        alis = alignments.readlines()
    colored = 0
    for ali in alis:
        gene = ali.split(".")[0]#.split("_")[1]
        ali = alignment + ali.strip()
        if ali.split(".")[-1] != "fasta" and ali.split(".")[-1] != "phy":
            continue

        log = ["#!/bin/bash",
               "# " + gene]

        log.append(commandline("mkdir output/" + gene))
        log.append(commandline("mkdir output/" + gene + "/01_input"))
        log.append(commandline("mkdir output/" + gene + "/02a_output_RAxML"))
        log.append(commandline("mkdir output/" + gene + "/02a_output_RAxML/CONSEL"))
        log.append(commandline("mkdir output/" + gene + "/02b_output_IQTree"))

        log.append(commandline("cp " + ali + " output/" + gene + "/01_input/"))

        log = log + ["\n",
                     "# RAxML"]

        start = time.time()
        log = log + raxml(ali, constraints, consel_path, model, gene)
        raxml_runtime = time.time() - start

        trees_raxml = []
        with open(gene + '_COMBINED.tre', "r") as multitree:
            for tree in multitree:
                trees_raxml.append(tree.strip())
        best_raxml = findBestTree(trees_raxml)
        trees_raxml = trees_raxml[1:]

        au_raxml = [0] * (len(constraints) + 1)
        with open(gene + "_CONSEL.consel") as consel:
            for line in consel.readlines()[3:len(constraints)+4]:
                au_raxml[int(line.split()[2]) - 1] = float(line.split()[4])
        au_raxml = au_raxml[1:]

        treeFile_raxml.write(gene + place[best_raxml] + " " + trees_raxml[best_raxml])

        for i in range(len(au_raxml)):
            if au_raxml[i] <= 0.05:
                 au_raxml[i] = str(au_raxml[i]) + "*"

        au_raxml[best_raxml] = str(au_raxml[best_raxml]) + "s"

        llh_raxml = []
        for llh_file in ["RAxML_log." + name + "_" + gene for name in  ["withoutConstraints"] + ["hypothesis" + str(i) for i in range(len(constraints))]]:
            with open(llh_file,"r") as llh:
                llh_raxml.append(float(llh.readlines()[-1].split()[-1]))

        log.append(commandline('mmv "RAxML*" ' + gene + '_RAxML#1'))
        log.append(commandline("mv " + gene + "_RAxML_* output/" + gene + "/02a_output_RAxML/"))
        log.append(commandline("mv " + gene + "_CONSEL* output/" + gene + "/02a_output_RAxML/CONSEL/"))

        log = log + ["\n\n", "# IQTree"]

        start = time.time()
        log = log + iqtree(ali, constraints, gene)
        iqtree_runtime = time.time() - start

        trees_iqtree = []
        with open(gene + '_COMBINED.tre', "r") as multitree:
            for tree in multitree:
                trees_iqtree.append(tree.strip())
        best_iqtree = findBestTree(trees_iqtree)
        trees_iqtree = trees_iqtree[1:]

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

        llh_iqtree = []
        for llh_file in [gene + "_IQTree_" + name + ".iqtree" for name in ["unconst"] + ["hypo" + str(i) for i in range(len(constraints))]]:
            with open(llh_file,"r") as llh:
                for line in llh.readlines():
                    if "Log-likelihood of the tree:" in line:
                        llh_iqtree.append(float(line.split()[4]))
                        break

        log.append(commandline("mv " + gene + "_IQTree* output/" + gene + "/02b_output_IQTree/"))

        treeFile_iqtree.write(gene + place[best_iqtree] + " " + trees_iqtree[best_iqtree] + "\n")

        for i in range(len(au_iqtree)):
            if au_iqtree[i] <= 0.05:
                 au_iqtree[i] = str(au_iqtree[i]) + "*"

        au_iqtree[best_iqtree] = str(au_iqtree[best_iqtree]) + "s"

        auFile.write(gene + "," + ''.join([str(au_raxml[i]) + "," for i in range(len(constraints))]) + ''.join([str(au_iqtree[i]) + "," for i in range(len(constraints))]) + str(raxml_runtime) + "," + str(iqtree_runtime) + "\n")
        #if colored % 2 == 0:
        #    latexFile.write(" \\rowcolor{black!20} ")
        #    llsFile.write(" \\rowcolor{black!20} ")

        #latexFile.write("\\textit{" + gene + "}" + " & " + str(edit_num(au_pylogeny[0])) + " & " + str(edit_num(au_pylogeny[1])) + " & " + str(edit_num(au_pylogeny[2])) + " & " + str(edit_num(au_iqtree[0])) + " & " + str(edit_num(au_iqtree[1])) + " & " + str(edit_num(au_iqtree[2])) + " & " + str(round(pylogeny_runtime,3)) + " & " + str(round(iqtree_runtime,3)) + "\\\\" + "\n")
        #llsFile.write("\\textit{" + gene + "}" + " & " + str(llh_pylogeny[0]) + " & " + str(llh_pylogeny[1]) + " & " + str(llh_pylogeny[2]) + " & " + str(llh_pylogeny[3]) + " & " + str(llh_iqtree[0]) + " & " + str(llh_iqtree[1]) + " & " + str(llh_iqtree[2]) + " & " + str(llh_iqtree[3]) + "\\\\\n")

        #if colored % 2 == 0:
        #    llsFile.write(" \\rowcolor{black!20} ")
        #llsFile.write("$\\Delta$" + " & " + "" + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[1],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[2],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[3],4)) + " & " + "" + " & " + str(round(llh_iqtree[0] - llh_iqtree[1],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[2],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[3],4)) + "\\\\\n")
        #colored = colored + 1

        with open("output/" + gene + "/" + gene + "_log.sh","w") as logFile:
            for line in log:
                logFile.write(line + "\n")

        os.system("rm " + gene + "_COMBINED.tre")
        print(str(colored) + " / " + str(len(alis)/2))


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

    #os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/au_runtime_table.tex ")
    #os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/likelihoods_table.tex")
    os.system("rm *.tmp")

def edit_num(value):
    value = str(value)
    value = value.replace("s","\\textsuperscript{s}").replace("*","\\textsuperscript{*}")
    if "*" in str(value):
        return "\\textbf{" + str(value) + "}"
    return value

def findBestTree(treeList):
    tns = dendropy.TaxonNamespace()
    rank = []
    unconst = dendropy.Tree.get_from_string(
            treeList[0],
            "newick",
            taxon_namespace=tns)
    unconst.encode_bipartitions()
    for tree in treeList[1:]:
        hypo = dendropy.Tree.get_from_string(
                tree,
                "newick",
                taxon_namespace=tns)
        hypo.encode_bipartitions()
        rank.append(dendropy.calculate.treecompare.euclidean_distance(unconst, hypo))
    return rank.index(min(rank))


def readConstraints(constraint_path):
    constraints_tmp = []

    ## put all constraint trees in one list
    with open(constraint_path) as constFile:
        for line in constFile:
            constraints_tmp.append(line.strip())
        constraints_tmp = ''.join(constraints_tmp).split(";")[:-1]

    ## deroot the trees
    constraints = []
    for const in constraints_tmp:
        tree = dendropy.Tree.get_from_string(const + ";", "newick")
        tree.deroot()
        constraints.append(tree)

    return constraints

def iqtree(alignment, constraints, gene_name):
    log = []
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -pre " + gene_name + "_IQTree_unconst -quiet"))

    for i in range(len(constraints)):
        with open("hypo.txt","w") as hypo:
            hypo.write(constraints[i].as_string(schema="newick"))
        log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -g hypo.txt -pre " + gene_name + "_IQTree_hypo" + str(i) + " -quiet"))
    log.append(commandline("cat " + gene_name + "_IQTree_unconst.treefile " + ''.join([gene_name + "_IQTree_hypo" + str(i) + ".treefile " for i in range(len(constraints))]) + "> " + gene_name + "_COMBINED.tre"))
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te " + gene_name + "_IQTree_unconst.treefile -zb 10000 -au -pre " + gene_name + "_IQTree -quiet"))
    return log

def raxml(alignment, constraints, consel_path, model, gene_name):
    log = []
    log.append(commandline("raxmlHPC -s " + alignment + " -n withoutConstraints_" + gene_name + " -m  " + model + " -p 10000 -w " + os.getcwd()) + " --silent")
    for i in range(len(constraints)):
        with open("hypo.txt","w") as hypo:
            hypo.write(constraints[i].as_string(schema="newick"))
        log.append(commandline("raxmlHPC -s " + alignment + " -n hypothesis" + str(i) + "_" + gene_name + " -m " + model + " -g hypo.txt -p 10000 -w " + os.getcwd()) + " --silent")

    log.append(commandline("cat RAxML_bestTree.withoutConstraints_" + gene_name + ''.join([" RAxML_bestTree.hypothesis" + str(i) + "_" + gene_name for i in range(len(constraints))]) + " > " + gene_name + "_COMBINED.tre"))
    log.append(commandline("raxmlHPC -s " + alignment + " -n " + gene_name + ".trees.sitelh -m "+ model + " -f g -t RAxML_bestTree.withoutConstraints_" + gene_name + " -z " + gene_name + "_COMBINED.tre -p 10000 -w " + os.getcwd()) + " --silent")
    log.append("\n## CONSEL")
    log.append(commandline("mv RAxML_perSiteLLs." + gene_name + ".trees.sitelh RAxML_perSiteLLs_" + gene_name + ".trees.sitelh"))
    log.append(commandline(consel_path + "/seqmt --puzzle RAxML_perSiteLLs_" + gene_name + ".trees.sitelh " + gene_name + "_CONSEL.mt"))
    log.append(commandline(consel_path + "/makermt " + gene_name + "_CONSEL.mt"))
    log.append(commandline(consel_path + "/consel " + gene_name + "_CONSEL.rmt"))
    log.append(commandline(consel_path + "/catpv " + gene_name + "_CONSEL.pv > " + gene_name + "_CONSEL.consel"))

    return log

def commandline(command):
    os.system(command)
    return command
