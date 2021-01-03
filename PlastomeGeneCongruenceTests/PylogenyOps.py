import os
import dendropy
from Bio import SeqIO
from Bio.Nexus import Nexus
from ete3 import Tree

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

def removableTaxa(alignment_path, allTaxa, format):
    Taxa = allTaxa.copy()
    if format == "nex":
        format = "nexus"
    elif format == "phy":
        format == "phylip"

    if(format in ["fasta","nexus","phylip"]):
        for seq_record in SeqIO.parse(alignment_path, format):
            try:
                Taxa.remove(seq_record.id)
            except:
                pass
        for seq in Taxa:
            try:
                allTaxa.remove(seq)
            except:
                pass
    else:
        print("'" + alignment_path + "' is not in a supported format")
    return allTaxa

def removeTaxa(constraint_path, allTaxa, tree):
    tree = Tree(tree)
    tree.prune(allTaxa)
    return(tree.write(format=9))

## Main Function
def removeGenesFromConstTree(alignment_path, constraint_path, output_path):
    if not output_path:
        new_constraint_path = open(''.join(constraint_path.split(".")[:-1]) + "_new." + constraint_path.split(".")[-1], "w")
    else:
        new_constraint_path = open(output_path, "w")
    with open(constraint_path, "r") as const:
        ## For each tree in tree file find the Taxa which are not part of the alignment
        for tree in const.readlines():
            stayTaxa = tree.strip().replace("(","").replace(")","").replace(";","").split(",")
            allTaxa = stayTaxa.copy()
            ## Find the Taxa which are part of all alignmentfiles
            #try:
            #    for ali in os.listdir(alignment_path):
            #        ali = alignment_path + "/" + ali
            #        stayTaxa = removableTaxa(ali, stayTaxa, ali.split(".")[-1])
            #except:
            stayTaxa = removableTaxa(alignment_path, stayTaxa, alignment_path.split(".")[-1])
            ## Print all removed Taxa
            print("removed Taxa:")
            i = 0
            for taxa in allTaxa:
                if taxa not in stayTaxa:
                    i = i + 1
                    print(taxa)
            print("Number removed Taxa: " + str(i))

            ## Write new constraint tree to file
            new_constraint_path.write(removeTaxa(constraint_path, stayTaxa, tree) + "\n")
    new_constraint_path.close()


def iqtree_autest(alignment, gene_name, mlcalc, threadNumber):
    log = []
    if(mlcalc == "RAxML"):
        log.append(commandline("iqtree -s " + alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te output/" + gene_name + "/02_output_RAxML/RAxML_bestTree.withoutConstraints_" + gene_name + " -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber))
    else:
        log.append(commandline("iqtree -s " + alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te output/" + gene_name + "/02_output_IQTree/" + gene_name + "_IQTree_unconst.treefile -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber))
    return log

def iqtree2_autest(alignment, iqtree2_path, gene_name, threadNumber):
    log = []
    if(mlcalc == "RAxML"):
        log.append(commandline(iqtree2_path + " -s " + alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te output/" + gene_name + "/02_output_RAxML/RAxML_bestTree.withoutConstraints_" + gene_name + " -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber))
    else:
        log.append(commandline(iqtree2_path + " -s " + alignment +" -m GTR+I+G -z " + gene_name + "_COMBINED.tre -te output/" + gene_name + "/02_output_IQTree/" + gene_name + "_IQTree_unconst.treefile -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber))

    return log

def iqtree_mltree(alignment, constraints, gene_name, threadNumber):
    log = []
    log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -pre " + gene_name + "_IQTree_unconst -quiet -nt " + threadNumber))

    for i in range(len(constraints)):
        with open("hypo.txt","w") as hypo:
            hypo.write(constraints[i].as_string(schema="newick", suppress_rooting=True))
        removeGenesFromConstTree(alignment, "hypo.txt", "hypo_rem.txt")
        log.append(commandline("iqtree -s "+ alignment +" -m GTR+I+G -g hypo_rem.txt -pre " + gene_name + "_IQTree_hypo" + str(i) + " -quiet -nt " + threadNumber))
    log.append(commandline("cat " + gene_name + "_IQTree_unconst.treefile " + ''.join([gene_name + "_IQTree_hypo" + str(i) + ".treefile " for i in range(len(constraints))]) + "> " + gene_name + "_COMBINED.tre"))
    return log

def raxml(alignment, constraints, model, gene_name, threadNumber):
    log = []
    log.append(commandline("raxmlHPC -s " + alignment + " -n withoutConstraints_" + gene_name + " -m  " + model + " -p 10000 -w " + os.getcwd()) + " --silent -T " + threadNumber)
    for i in range(len(constraints)):
        with open("hypo.txt","w") as hypo:
            hypo.write(constraints[i].as_string(schema="newick", suppress_rooting=True))
        removeGenesFromConstTree(alignment, "hypo.txt", "hypo_rem.txt")
        log.append(commandline("raxmlHPC -s " + alignment + " -n hypothesis" + str(i) + "_" + gene_name + " -m " + model + " -g hypo_rem.txt -p 10000 -w " + os.getcwd()) + " --silent -T " + threadNumber)

    log.append(commandline("cat RAxML_bestTree.withoutConstraints_" + gene_name + ''.join([" RAxML_bestTree.hypothesis" + str(i) + "_" + gene_name for i in range(len(constraints))]) + " > " + gene_name + "_COMBINED.tre"))
    return log

def consel(alignment, consel_path, model, gene_name, mlcalc, threadNumber):
    log = []
    if(mlcalc == "RAxML"):
        log.append(commandline("raxmlHPC -s " + alignment + " -n " + gene_name + ".trees.sitelh -m "+ model + " -f g -t output/" + gene_name + "/02_output_RAxML/RAxML_bestTree.withoutConstraints_" + gene_name + " -z " + gene_name + "_COMBINED.tre -p 10000 -w " + os.getcwd()) + " --silent -T " + threadNumber)
    else:
        log.append(commandline("raxmlHPC -s " + alignment + " -n " + gene_name + ".trees.sitelh -m "+ model + " -f g -t output/" + gene_name + "/02_output_IQTree/" + gene_name + "_IQTree_unconst.treefile -z " + gene_name + "_COMBINED.tre -p 10000 -w " + os.getcwd()) + " --silent -T " + threadNumber)
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
