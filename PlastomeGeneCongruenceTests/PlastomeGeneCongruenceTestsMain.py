#!/usr/bin/env python3.7

''' Main operations in TreeTopology '''

#####################
# IMPORT OPERATIONS #
#####################

import IOOps as IOOps
import PylogenyOps as PylOps
import os
from pylogeny.tree import tree, treeSet
from pylogeny.alignment import phylipFriendlyAlignment as ali
import time
import dendropy
from dendropy.calculate import treecompare


###############
# AUTHOR INFO #
###############

__author__ = 'Nils Jenke; Yannick Hartmaring'
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
                 constraint,
                 consel_path,
                 model):

    os.system("mkdir output/")
    os.system("mkdir output/SUMMARY")

    place = ["hypoA", "hypoB", "hypoC"]
    treeFile_pylogeny = open("output/SUMMARY/pylogeny_hypoTreeShortestDistUnconstTree.tre","w")
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
    latexFile.write("gene & \\multicolumn{3}{c}{Pylogeny}  & \\multicolumn{3}{c}{IQTree} & \\multicolumn{2}{c}{Runtime in seconds}\\\\\n")
    latexFile.write(" & Hypothesis A & Hypothesis B & Hypothesis C & Hypothesis A & Hypothesis B & Hypothesis C & Pylogeny & IQTree\\\\\n")

    llsFile.write("\\documentclass[a4paper]{article}\n")
    llsFile.write("\\usepackage{colortbl, geometry, longtable}\n")
    llsFile.write("\\usepackage[cmyk,table]{xcolor}\n")
    llsFile.write("\\geometry{paperheight=297mm, paperwidth=210mm, margin=2pt}\n")
    llsFile.write("\\pagenumbering{gobble}\n")
    llsFile.write("\\begin{document}\n")
    llsFile.write("\\footnotesize\n")
    llsFile.write("\\begin{longtable}{p{0.03\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}|p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}p{0.095\\linewidth}}\n")
    llsFile.write("gene & \\multicolumn{4}{c}{Pylogeny} & \\multicolumn{4}{c}{IQTree}\\\\\n")
    llsFile.write(" & Unconstraint & Hypothesis A & Hypothesis B & Hypothesis A & Unconstraint & Hypothesis A & Hypothesis B & Hypothesis C\\\\\n")
    llsFile.write("\\endhead\n")

    auFile.write("gene,hypoA_pylogeny,hypoB_pylogeny,hypoC_pylogeny,hypoA_iqtree,hypoB_iqtree,hypoC_iqtree,runtime_pylogeny,runtime_iqtree\n")

    os.system("ls " + alignment + " > aligments.tmp")
    with open("aligments.tmp","r") as alignments:
        alis = alignments.readlines()
    colored = 0
    for ali in alis:
        gene = ali.split(".")[0].split("_")[1]
        ali = alignment + ali.strip()
        if ali.split(".")[-1] != "fasta" and ali.split(".")[-1] != "phy":
            continue

        log = ["#!/bin/bash",
               "# " + gene]

        log.append(commandline("mkdir output/" + gene))
        log.append(commandline("mkdir output/" + gene + "/01_input"))
        log.append(commandline("mkdir output/" + gene + "/02a_output_Pylogeny"))
        log.append(commandline("mkdir output/" + gene + "/02a_output_Pylogeny/RAxML"))
        log.append(commandline("mkdir output/" + gene + "/02a_output_Pylogeny/CONSEL"))
        log.append(commandline("mkdir output/" + gene + "/02b_output_IQTree"))

        log.append(commandline("cp " + ali + " output/" + gene + "/01_input/"))

        log = log + ["\n",
                     "# Pylogeny",
                     "## RAxML"]

        start = time.time()
        log = log + raxml(ali, constraint, consel_path, model, gene)
        pylogeny_runtime = time.time() - start


        trees_pylogeny = []
        with open(gene + '_COMBINED.tre', "r") as multitree:
            for tree in multitree:
                trees_pylogeny.append(tree.strip())
        best_pylogeny = findBestTree(trees_pylogeny)
        trees_pylogeny = trees_pylogeny[1:]

        au_pylogeny = [0,0,0,0]
        with open(gene + "_CONSEL.consel") as consel:
            for line in consel.readlines()[3:7]:
                au_pylogeny[int(line.split()[2]) - 1] = float(line.split()[4])
        au_pylogeny = au_pylogeny[1:]

        treeFile_pylogeny.write(gene + place[best_pylogeny] + " " + trees_pylogeny[best_pylogeny])

        for i in range(len(au_pylogeny)):
            if au_pylogeny[i] <= 0.05:
                 au_pylogeny[i] = str(au_pylogeny[i]) + "*"

        au_pylogeny[best_pylogeny] = str(au_pylogeny[best_pylogeny]) + "s"

        llh_pylogeny = []
        for llh_file in ["RAxML_log." + name + "_" + gene for name in ["withoutConstraints","hypothesisA","hypothesisB","hypothesisC"]]:
            print llh_file
            with open(llh_file,"r") as llh:
                llh_pylogeny.append(float(llh.readlines()[-1].split()[-1]))

        log.append(commandline('mmv "RAxML*" ' + gene + '_RAxML#1'))
        log.append(commandline("mv " + gene + "_RAxML_* output/" + gene + "/02a_output_Pylogeny/RAxML/"))
        log.append(commandline("mv " + gene + "_CONSEL* output/" + gene + "/02a_output_Pylogeny/CONSEL/"))

        log = log + ["\n\n", "# IQTree"]

        start = time.time()
        log = log + iqtree(ali, constraint, gene)
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
        for llh_file in [gene + "_IQTree_" + name + ".iqtree" for name in ["unconst","hypoA","hypoB","hypoC"]]:
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

        auFile.write(gene + "," + str(au_pylogeny[0]) + "," + str(au_pylogeny[1]) + "," + str(au_pylogeny[2]) + "," + str(au_iqtree[0]) + "," + str(au_iqtree[1]) + "," + str(au_iqtree[2]) + "," + str(pylogeny_runtime) + "," + str(iqtree_runtime) + "\n")
        if colored % 2 == 0:
            latexFile.write(" \\rowcolor{black!20} ")
            llsFile.write(" \\rowcolor{black!20} ")

        latexFile.write("\\textit{" + gene + "}" + " & " + str(edit_num(au_pylogeny[0])) + " & " + str(edit_num(au_pylogeny[1])) + " & " + str(edit_num(au_pylogeny[2])) + " & " + str(edit_num(au_iqtree[0])) + " & " + str(edit_num(au_iqtree[1])) + " & " + str(edit_num(au_iqtree[2])) + " & " + str(round(pylogeny_runtime,3)) + " & " + str(round(iqtree_runtime,3)) + "\\\\" + "\n")
        llsFile.write("\\textit{" + gene + "}" + " & " + str(llh_pylogeny[0]) + " & " + str(llh_pylogeny[1]) + " & " + str(llh_pylogeny[2]) + " & " + str(llh_pylogeny[3]) + " & " + str(llh_iqtree[0]) + " & " + str(llh_iqtree[1]) + " & " + str(llh_iqtree[2]) + " & " + str(llh_iqtree[3]) + "\\\\\n")

        if colored % 2 == 0:
            llsFile.write(" \\rowcolor{black!20} ")
        llsFile.write("$\\Delta$" + " & " + "" + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[1],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[2],4)) + " & " + str(round(llh_pylogeny[0] - llh_pylogeny[3],4)) + " & " + "" + " & " + str(round(llh_iqtree[0] - llh_iqtree[1],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[2],4)) + " & " + str(round(llh_iqtree[0] - llh_iqtree[3],4)) + "\\\\\n")
        colored = colored + 1

        with open("output/" + gene + "/" + gene + "_log.sh","w") as logFile:
            for line in log:
                logFile.write(line + "\n")

        os.system("rm " + gene + "_COMBINED.tre")
        print str(colored) + " / " + str(len(alis)/2)


    latexFile.write("\\end{tabular}\\\\\n")
    latexFile.write("\\textsuperscript{s}tree with lowest distance to unconstraint tree; \\textsuperscript{*}p-value $\\leq$ 0.05\n")
    latexFile.write("\\end{document}\n")

    llsFile.write("\\end{longtable}\n")
    llsFile.write("\\end{document}\n")

    llsFile.close()
    latexFile.close()
    treeFile_pylogeny.close()
    treeFile_iqtree.close()
    auFile.close()

    os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/au_runtime_table.tex ")
    os.system("xelatex -output-directory output/SUMMARY/ output/SUMMARY/likelihoods_table.tex")
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
        rank.append(treecompare.euclidean_distance(unconst, hypo))
    return rank.index(min(rank))


def iqtree(alignment, constraint, gene_name):
    log = []
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -pre " + gene_name + "_IQTree_unconst -quiet"))
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -g "+ constraint + "/constraints_hypothesisA.txt -pre " + gene_name + "_IQTree_hypoA -quiet"))
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -g "+ constraint + "/constraints_hypothesisB.txt -pre " + gene_name + "_IQTree_hypoB -quiet"))
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -g "+ constraint + "/constraints_hypothesisC.txt -pre " + gene_name + "_IQTree_hypoC -quiet"))
    log.append(commandline("cat " + gene_name + "_IQTree_unconst.treefile " + gene_name + "_IQTree_hypoA.treefile " + gene_name + "_IQTree_hypoB.treefile " + gene_name + "_IQTree_hypoC.treefile > " + gene_name + "_COMBINED.tre"))
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te " + gene_name + "_IQTree_unconst.treefile -zb 10000 -au -pre " + gene_name + "_IQTree -quiet"))
    return log

def raxml(alignment, constraint, consel_path, model, gene_name):
    log = []
    log.append(commandline("raxmlHPC -s " + alignment + " -n withoutConstraints_" + gene_name + " -m  " + model + " -p 10000 -w ."))
    log.append(commandline("raxmlHPC -s " + alignment + " -n hypothesisA_" + gene_name + " -m " + model + " -g " + constraint + "/_hypothesisA.txt -p 10000 -w ."))
    log.append(commandline("raxmlHPC -s " + alignment + " -n hypothesisB_" + gene_name + " -m " + model + " -g " + constraint + "/_hypothesisB.txt -p 10000 -w ."))
    log.append(commandline("raxmlHPC -s " + alignment + " -n hypothesisC_" + gene_name + " -m " + model + " -g " + constraint + "/_hypothesisC.txt -p 10000 -w ."))
    log.append(commandline("cat RAxML_bestTree.withoutConstraints_" + gene_name + " RAxML_bestTree.hypothesisA_" + gene_name + " RAxML_bestTree.hypothesisB_" + gene_name + " RAxML_bestTree.hypothesisC_" + gene_name + " > " + gene_name + "_COMBINED.tre"))
    log.append(commandline("raxmlHPC -s " + alignment + " -n " + gene_name + ".trees.sitelh -m "+ model + " -f g -t RAxML_bestTree.withoutConstraints_" + gene_name + " -z " + gene_name + "_COMBINED.tre -p 10000 -w ."))

    log.append("\n## CONSEL")

    a = ali(alignment)
    treeset = treeSet()
    with open(gene_name + '_COMBINED.tre', "r") as multitree:
        for tree in multitree:
            treeset.addTreeByNewick(tree.strip())

    log.append(commandline("mv RAxML_perSiteLLs." + gene_name + ".trees.sitelh RAxML_perSiteLLs_" + gene_name + ".trees.sitelh"))


    return log

def pylogeny(alignment, constraint, consel_path, model, gene_name):
    log = []
    # Pylogeny unconstraint
    log.append(PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file='withoutConstraints_' + gene_name,
                   parsimony_seed=10000).runFunction())

    # Pylogeny Hypothesis A
    log.append(PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file='hypothesisA_' + gene_name,
                   grouping_constraint=constraint + "/constraints_hypothesisA.txt",
                   parsimony_seed=10000).runFunction())

    # Pylogeny Hypothesis B
    log.append(PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file='hypothesisB_' + gene_name,
                   grouping_constraint=constraint + "/constraints_hypothesisB.txt",
                   parsimony_seed=10000).runFunction())

    # Pylogeny Hypothesis C
    log.append(PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   out_file='hypothesisC_' + gene_name,
                   grouping_constraint=constraint + "/constraints_hypothesisC.txt",
                   parsimony_seed=10000).runFunction())


    log.append(commandline("cat RAxML_bestTree.withoutConstraints_" + gene_name + " RAxML_bestTree.hypothesisA_" + gene_name + " RAxML_bestTree.hypothesisB_" + gene_name + " RAxML_bestTree.hypothesisC_" + gene_name + " > "+ gene_name + "_COMBINED.tre"))

    log.append(PylOps.raxml_2(inp_align=alignment,
                   model=model,
                   alg='g',
                   startingTree='RAxML_bestTree.withoutConstraints_' + gene_name,
                   bipartition_filename=gene_name + '_COMBINED.tre',
                   out_file=gene_name + ".trees.sitelh",
                   parsimony_seed=10000,
                   slow=True,
                   numboot=10000).runFunction())

    log.append("\n## CONSEL")

    a = ali(alignment)
    treeset = treeSet()
    with open(gene_name + '_COMBINED.tre', "r") as multitree:
        for tree in multitree:
            treeset.addTreeByNewick(tree.strip())

    log.append(commandline("mv RAxML_perSiteLLs." + gene_name + ".trees.sitelh RAxML_perSiteLLs_" + gene_name + ".trees.sitelh"))
    log = log + (PylOps.consel_2(treeset, a,  "RAxML_perSiteLLs_" + gene_name + ".trees.sitelh", gene_name + "_CONSEL", consel_path).getLog())
    return log

def commandline(command):
    os.system(command)
    return command
