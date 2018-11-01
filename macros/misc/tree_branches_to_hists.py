import sys
from ROOT import *
import argparse
import array

def branches_to_hist(filename, treename, outname):
    fin  = TFile(filename, "READ")
    tree = fin.Get(treename)
    fout = TFile(outname , "RECREATE")
    fout.cd()
    for branch in tree.GetListOfBranches():
        bname = branch.GetName()
        draw_command = bname + ">>" + bname
        try:
            tree.Draw(draw_command)
            h = ROOT.gPad.GetPrimitive(bname)
            h.Write()
            print "Saved branch", bname
            del h
        except:
            print "Could not draw branch", bname, ". skipping this"
        del bname    
    fout.Close()
    fin.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Read all branches in tree into hists and saves to an output file')
    parser.add_argument("-f", "--filename", dest="filename", type=str, help="ROOT filename")
    parser.add_argument("-t", "--treename", dest="treename", type=str, help="Tree name")
    parser.add_argument("-o", "--outname" , dest="outname" , type=str, help="Name of output file")
    args = parser.parse_args()

    branches_to_hist(args.filename, args.treename, args.outname)

