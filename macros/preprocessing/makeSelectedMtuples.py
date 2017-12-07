import os, sys
import numpy as np
from ROOT import *
import argparse

M_e = 0.510998
M_m = 105.65836
M_n = 939.565346
M_p = 938.272013
M_d = 1232.0
Eb = 27.0
neut_2p2h_mode = 2

def elike_cut(sk_tree):
    return sk_tree.fqevis > 100 and sk_tree.erec < 1250

def mulike_cut(sk_tree):
    return sk_tree.fqmomm > 200

def erec_ccqe(p_lep, costheta, M_lep, numode=True):
    if numode:
        M_1 = M_p
        M_2 = M_n
    else:
        M_1 = M_n
        M_2 = M_p
    E_lep = np.sqrt(p_lep**2 + M_lep**2)
    E_rec = (M_1**2 - (M_2 - Eb)**2 - M_lep**2 + 2 * (M_2 - Eb) * E_lep) / (2 * (M_2 - Eb - E_lep + p_lep * costheta))
    return E_rec

def erec_cc1pi(p_lep, costheta, M_lep):
    E_lep = np.sqrt(p_lep**2 + M_lep**2)
    E_rec = M_p * E_lep - 0.5 * M_lep**2
    E_rec += 0.5 * (M_d**2 - M_p**2)
    E_rec /= (M_p - E_lep + p_lep * costheta)
    return E_rec

def erec(sk_tree, numode=true):
    dirbeam = np.array([0.669764, -0.742179, 0.024223])
    dirlep = np.reshape(sk_tree.fq1rdir, (-1, 7, 3))
    if sk_tree.mode == neut_2p2h_mode:
        nlep = 3
    else:
        nlep = 2
    if sk_tree.iclass == 11:
        M_lep = M_e
        costheta = dirbeam[0] * dirlep[0][1][0] + dirbeam[1] * dirlep[0][1][1] + dirbeam[2] * dirlep[0][1][2]
        erec_new = erec_ccqe(sk_tree.fqmome, costheta, M_lep, numode)
        if not np.abs(sk_tree.erec - erec_new) < 1E-8:
            print "1Re erec_sktree:", sk_tree.erec, "erec_computed:", erec_new
    elif sk_tree.iclass == 14:
        M_lep = M_e
        costheta = dirbeam[0] * dirlep[0][1][0] + dirbeam[1] * dirlep[0][1][1] + dirbeam[2] * dirlep[0][1][2]
        erec_new = erec_cc1pi(sk_tree.fqmome, costheta, M_lep)
        if not np.abs(sk_tree.erec - erec_new) < 1E-8:
            print "CC1pi erec_sktree:", sk_tree.erec, "erec_computed:", erec_new
    elif sk_tree.iclass == 13:
        M_lep = M_m
        costheta = dirbeam[0] * dirlep[0][2][0] + dirbeam[1] * dirlep[0][2][1] + dirbeam[2] * dirlep[0][2][2]
        erec_new = erec_ccqe(sk_tree.fqmomm, costheta, M_lep, numode)
        if not np.abs(sk_tree.erec - erec_new) < 1E-8:
            print "Rmu erec_sktree:", sk_tree.erec, "erec_computed:", erec_new

if __name__ == "__main__":
    # Get arguments (ie. SK file name)
    parser = argparse.ArgumentParser(description='Splits an SK file into two files, containing the numu- and nue-selected events')
    parser.add_argument('filename', nargs=1)
    args = parser.parse_args()
    skfile = args.filename[0]

    # Get SK tree
    f_sk = TFile(skfile, "OPEN")
    if not f_sk:
        print "Error: could not open file ", skfile
        sys.exit()
    sk_tree = f_sk.Get("minituple")
    if not sk_tree:
        sk_tree = f_sk.Get("h1")

    if "fhc" in skfile:
        numode = True
    else:
        numode = False
    if not numode:
        print "Forcing RHC ereco calculation to use same equation as FHC"
        numode = True

    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print "Check: number of events in original SK tree = ", n_total

    # Get the SK filename (without '.root')
    # skfile = skfile[:-5]

    # Make new filenames for files containing only the numu/nue/cc1pi trees
    out_filenames = dict()
    out_filenames["numu"] = skfile[:-5]+"_numuselec.root"
    out_filenames["nue"] = skfile[:-5]+"_nueselec.root"
    out_filenames["cc1pi"] = skfile[:-5]+"_cc1piselec.root"

    print "Making three new files: "
    for key in out_filenames:
        print out_filenames[key]

    out_files = dict()
    out_trees = dict()
    for key in out_filenames:
        out_files[key] = TFile(out_filenames[key], "RECREATE")
        out_trees[key] = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.iclass == 13: # numu fitqun selection
            if mulike_cut(sk_tree):
                erec(sk_tree, numode)
                out_trees["numu"].Fill()
        elif sk_tree.iclass == 11: # nue+1Rpi0 fitqun selection
            if elike_cut(sk_tree):
                erec(sk_tree, numode)
                out_trees["nue"].Fill()
        elif sk_tree.iclass == 14: # nue cc1pi fitqun selection
            if elike_cut(sk_tree):
                erec(sk_tree, numode)
                out_trees["cc1pi"].Fill()

    for key in out_filenames:
        out_files[key].cd()
        out_trees[key].Write("mtuple")
        out_files[key].Close()

    # Check number of events in output trees
    n_events = dict()
    for key in out_filenames:
        file = TFile(out_filenames[key], "READ")
        n_events[key] = file.Get("mtuple").GetEntries()
        print "# selected", key, "events:", n_events[key]
        file.Close()

    print "Sum of selected events =", n_events["numu"] + n_events["nue"] + n_events["cc1pi"]
