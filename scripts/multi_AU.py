######
#
# Quick and dirty python skript to run AU Test multiple times on RAxML and IQTree output
#
#####

import sys
import os
import time

def iqtree_autest(alignment, gene_name, threadNumber, iqtree_path, combinedtrees, starttree):
    os.system(iqtree_path + " -s " + alignment +" -m GTR+I+G -z " + combinedtrees + " -te " + starttree + " -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber)

def iqtree2_autest(alignment, gene_name, threadNumber, iqtree2_path, combinedtrees, starttree):
    os.system(iqtree2_path + " -s " + alignment +" -m GTR+I+G -z " + combinedtrees + " -te " + starttree + " -zb 10000 -au -pre " + gene_name + "_IQTree -quiet -nt " + threadNumber)

def consel_autest(consel_path, gene_name, sitelh):
    os.system(consel_path + "/seqmt --puzzle " + sitelh + " " + gene_name + "_CONSEL.mt")
    os.system(consel_path + "/makermt " + gene_name + "_CONSEL.mt")
    os.system(consel_path + "/consel " + gene_name + "_CONSEL.rmt")
    os.system(consel_path + "/catpv " + gene_name + "_CONSEL.pv > " + gene_name + "_CONSEL.consel")

def main(consel_path, iqtree_path, iqtree2_path, raxml_path, alignments, runs, constraints, mlcalc):
    genes = os.listdir(alignments)
    for i in range(runs):

        # for a better view create for each run a new file
        os.mkdir("run_" + str(i))
        auFile = open("run_" + str(i) + "/au_runtime_table.csv","w")
        auFile.write("gene," + ','.join([','.join([program + "_hypo" + str(i) for program in ["CONSEL","IQTree","IQTree2"]]) for i in range(constraints)]) + "," + ','.join(["runtime_" + program for program in ["CONSEL","IQTree","IQTree2"]]) + "\n")

        ## go throw all genes
        for gene in genes:
            # the SUMMARY file is not a gene
            if gene == "SUMMARY":
                continue

            # safe the input alignment for thhe current gene
            input = alignments + "/" + gene + "/01_input/" + gene + ".fasta"

            # for a better view create for each gene a new file
            os.mkdir("run_" + str(i) + "/" + gene)
            os.mkdir("run_" + str(i) + "/" + gene + "/CONSEL")
            os.mkdir("run_" + str(i) + "/" + gene + "/IQTREE")
            os.mkdir("run_" + str(i) + "/" + gene + "/IQTREE2")

            # create the combined tree file and set the starttree
            if(mlcalc == "IQTree"):
                starttree = alignments + "/" + gene + "/02_output_IQTree/" + gene + "_IQTree_unconst.treefile"
                os.system("cat " + starttree + " " + ''.join([alignments + "/" + gene + "/02_output_IQTree/" + gene + "_IQTree_hypo" + str(i) + ".treefile " for i in range(constraints)]) + "> " + gene + "_COMBINED.tre")
                combinedtrees = gene + "_COMBINED.tre"
            else:
                starttree = alignments + "/" + gene + "/02_output_RAxML/RAxML_bestTree.withoutConstraints_" + gene
                os.system("cat " + starttree + " " ''.join([" " + alignments + "/" + gene + "/02_output_RAxML/RAxML_bestTree.hypothesis" + str(i) + "_" + gene for i in range(constraints)]) + " > " + gene + "_COMBINED.tre")
                combinedtrees = gene + "_COMBINED.tre"

            # safe the file with sitelikelihoods to a variable, if the file does not exist: create it
            sitelh = alignments + "/" + gene + "/03a_output_CONSEL/RAxML_perSiteLLs_" + gene + ".trees.sitelh"
            if not os.path.exists(sitelh):
                os.system(raxml_path + " -s " + input + " -n " + gene + "_CONSEL.trees.sitelh -m GTRGAMMAI -f g -t " + starttree + " -z " + gene + "_COMBINED.tre -p 10000 -w " + os.getcwd())
                os.system("mv RAxML_perSiteLLs." + gene + "_CONSEL.trees.sitelh RAxML_perSiteLLs_" + gene + "_CONSEL.trees.sitelh")
                sitelh = "RAxML_perSiteLLs_" + gene + "_CONSEL.trees.sitelh"

            ## CONSEL
            start = time.time()
            consel_autest(consel_path, gene, sitelh)
            consel_runtime = round(time.time() - start,3)

            ## Save the AU Test values from CONSEL to a list
            au_consel = [0] * (constraints + 1)
            with open(gene + "_CONSEL.consel") as consel:
                for line in consel.readlines()[3:constraints+4]:
                    au_consel[int(line.split()[2]) - 1] = float(line.split()[4])
                au_consel = au_consel[1:]

            os.system("mv " + gene + "_CONSEL* run_" + str(i) + "/" + gene + "/CONSEL")

            ### IQTREE
            start = time.time()
            iqtree_autest(input, gene, "4", iqtree_path, combinedtrees, starttree)
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

            os.system("mv " + gene + "_IQTree* run_" + str(i) + "/" + gene + "/IQTREE")


            ## IQTREE2
            start = time.time()
            iqtree2_autest(input, gene, "4", iqtree2_path, combinedtrees, starttree)
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

            os.system("mv " + gene + "_IQTree* run_" + str(i) + "/" + gene + "/IQTREE2")
            os.remove(gene + "_COMBINED.tre")

            auFile.write(gene + "," + ''.join([str(au_consel[i]) + "," + str(au_iqtree[i]) + "," + str(au_iqtree2[i]) + "," for i in range(constraints)]) + str(consel_runtime) + "," + str(iqtree_runtime) + "," + str(iqtree2_runtime) + "\n")

        auFile.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]), int(sys.argv[7]), sys.argv[8])
