import sys
from ROOT import *
import argparse
import array

def index_tree(filename, treename):
    f = TFile(filename, "READ")
    tree = f.Get(treename)
    for branch in tree.GetListOfBranches():
        print branch
    f.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Applies the momentum, erec and evis cuts to the mtuple')
    parser.add_argument("-f", "--filename", dest="filename", type=str, help="ROOT filename")
    parser.add_argument("-t", "--treename", dest="treename", type=str, help="Tree name")
    args = parser.parse_args()

    index_tree(args.filename, args.treename)

