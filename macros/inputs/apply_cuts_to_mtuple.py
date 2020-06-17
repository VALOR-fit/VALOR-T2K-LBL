import sys
from ROOT import *
import argparse

def elike_cut(evis, erec):
    return evis > 100 and erec < 1250

def mulike_cut(p):
    return p > 200

def apply_cuts(sk_tree):
    out_tree = sk_tree.CloneTree(0)
    for event in sk_tree:
        if event.iclass == 13: # numu fitqun selection
            if mulike_cut(event.fqmomm):
                out_tree.Fill()
        elif event.iclass == 11: # nue+1Rpi0 fitqun selection
            if elike_cut(event.fqevis, event.erec):
                out_tree.Fill()
        elif event.iclass == 14: # nue cc1pi fitqun selection
            if elike_cut(event.fqevis, event.erec):
                out_tree.Fill()
    return out_tree

def count_events(sk_tree):
    num_fhc1re = 0
    num_rhc1re = 0
    num_fhc1rmu = 0
    num_rhc1rmu = 0
    num_fhc1pi = 0
    for event in sk_tree:
        if event.iclass == 13:
            if event.horn_current > 0:
                num_fhc1rmu += 1
            else:
                num_rhc1rmu += 1
        elif event.iclass == 11:
            if event.horn_current > 0:
                num_fhc1re += 1
            else:
                num_rhc1re += 1
        elif event.iclass == 14:
            if event.horn_current > 0:
                num_fhc1pi += 1
    print "#FHC 1Re:", num_fhc1re
    print "#RHC 1Re:", num_rhc1re
    print "#FHC 1pi:", num_fhc1pi
    print "#FHC 1Rmu:", num_fhc1rmu
    print "#RHC 1Rmu:", num_rhc1rmu

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Applies the momentum, erec and evis cuts to the mtuple')
    parser.add_argument("-f", "--in-filename", dest="in_filename", type=str, help="Input mtuple filename")
    parser.add_argument("-o", "--out-filename", dest="out_filename", type=str, help="Output mtuple filename")
    args = parser.parse_args()

    in_file = TFile(args.in_filename, "OPEN")
    if not in_file:
        print "Error: could not open file ", args.in_filename
        sys.exit()
    sk_tree = in_file.Get("minituple")
    if not sk_tree:
        print "Error: could not find minituple object"
        sys.exit()

    out_file = TFile(args.out_filename, "RECREATE")
    out_tree = apply_cuts(sk_tree)
    count_events(out_tree)
    out_tree.Write("minituple")
    out_file.Close()
    in_file.Close()

