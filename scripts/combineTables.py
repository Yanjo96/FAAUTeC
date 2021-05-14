__author__ = 'Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2020 Project Gruenstaeudl'
__info__ = 'combineTables -- Combine table with an equal column'
__version__ = '2020.10.30'

import sys
import argparse

## Parser
def parser():
    '''Simple argparse for an easier usage'''
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser._action_groups.pop()

    ### Arguments ###
    parser.add_argument('-1',
                        '--first',
                        help='absolute path to the folder containing the alignmentfiles; infile in PHYLIP or FASTA format; Example: /path_to_input/test.phy',
                        required=True)
    parser.add_argument('-2',
                        '--second',
                        help='absolute path to constraint file',
                        required=True)

    parser.add_argument('-o',
                        '--output',
                        help='absolute path to constraint file',
                        required=True)

    parser.add_argument('-d',
                        '--delimeter',
                        help='output file',
                        default=',',
                        required=False)


    args = parser.parse_args()

    main(args.first,args.second, args.output, args.delimeter)

def main(first_path, second_path, output_path, delimeter):
    try:
        first = open(first_path, "r")
        firstTable = [l.strip() for l in first.readlines()]
        first.close()
    except:
        print(first_path + " is not a correct path")
        sys.exit(0)
    try:
        second = open(second_path, "r")
        secondTable = [l.strip() for l in second.readlines()]
        second.close()
    except:
        print(second_path + " is no a correct path")
        sys.exit(0)

    with open(output_path, "w") as out:
        for line in firstTable:
            line = line.split(delimeter)
            for secondLine in secondTable:
                secondLine = secondLine.split(delimeter)
                if line[0] == secondLine[0]:
                    out.write(line[0] + delimeter + delimeter.join(line[1:]) + delimeter + delimeter.join(secondLine[1:]) + "\n")
                    break
parser()
