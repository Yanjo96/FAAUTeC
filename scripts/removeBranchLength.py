#/bin/python3

__author__ = 'Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2020 Project Gruenstaeudl'
__info__ = 'Remove Branch length from trees'
__version__ = '2020.10.18'

import argparse
from ete3 import Tree

## Parser
def parser():
    '''Simple argparse for an easier usage'''
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser._action_groups.pop()

    ### Arguments ###
    parser.add_argument('-t',
                        '--tree',
                        help='absolute path to the folder containing the alignmentfiles; infile in PHYLIP or FASTA format; Example: /path_to_input/test.phy',
                        required=True)

    parser.add_argument('-o',
                        '--out',
                        help='output file',
                        required=True)

    parser.add_argument('-r',
                        '--remove',
                        help='String with character wich should removed from the trees. For example: -r /&%$',
                        default='',
                        required=False)


    args = parser.parse_args()

    main(args.tree, args.remove, args.out)


## Main Function
def main(tree_path, remove, output_path):
    with open(output_path, "w") as out:
        with open(tree_path, "r") as infile:
            for tree in infile:
                tree = tree.strip()
                for c in remove:
                    tree = tree.replace(c, "")
            t = Tree(tree)
            s = t.write(format=8).replace("NoName","")
            out.write(s + '\n')


parser()
