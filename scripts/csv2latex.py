#/bin/python3

__author__ = 'Yannick Hartmaring'
__copyright__ = 'Copyright (C) 2020 Project Gruenstaeudl'
__info__ = 'Remove Taxa which are not part of the Contraints'
__version__ = '2020.11.09'


## Parser
def parser():
    '''Simple argparse for an easier usage'''
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser._action_groups.pop()

    ### Arguments ###
    parser.add_argument('-c',
                        '--csv',
                        help='path to the csv table',
                        required=True)

    parser.add_argument('-o',
                        '--output',
                        help='output file',
                        required=True)

    args = parser.parse_args()

    IOOps.Outp().createLatex(args.csv, args.output)

if __name__ == '__main__':

    import sys, os
    import argparse
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'PlastomeGeneCongruenceTests'))
    #os.environ['LD_LIBRARY_PATH'] = "/home/yannick/Programs/tree-puzzle-5.3.rc16/src/"
    import IOOps

    parser()
