import os, sys
import array
import numpy as np
from ROOT import *
import argparse

import csv

M_e = 0.510998
M_m = 105.65836
M_n = 939.565346
M_p = 938.272013
M_d = 1232.0
Eb = 27.0
neut_2p2h_mode = 2

erec_diff_m = []
erec_diff_e = []

def elike_cut(sk_tree):
    return sk_tree.fqevis > 100 and sk_tree.erec < 1250

def mulike_cut(fqmomm):
    return fqmomm > 200

# Current CCQE Erec calculation uses FHC equation in RHC mode as well
def erec_ccqe(p_lep, costheta, M_lep, numode=True):
    M_1 = M_p
    M_2 = M_n
#    if numode:
#        M_1 = M_p
#        M_2 = M_n
#    else:
#        M_1 = M_n
#        M_2 = M_p
    E_lep = np.sqrt(p_lep**2 + M_lep**2)
    E_rec = (M_1**2 - (M_2 - Eb)**2 - M_lep**2 + 2 * (M_2 - Eb) * E_lep) / (2 * (M_2 - Eb - E_lep + p_lep * costheta))
    return E_rec

def erec_cc1pi(p_lep, costheta, M_lep):
    E_lep = np.sqrt(p_lep**2 + M_lep**2)
    E_rec = M_p * E_lep - 0.5 * M_lep**2
    E_rec += 0.5 * (M_d**2 - M_p**2)
    E_rec /= (M_p - E_lep + p_lep * costheta)
    return E_rec

def coulomb_correction(sk_tree, numode=true):
    dirbeam = np.array([0.669764, -0.742179, 0.024223])
    dirlep = np.reshape(sk_tree.fq1rdir, (-1, 7, 3))
    # coulomb correction applies only to mu
    M_lep = M_m
    # Apply Coulomb correction to muon momentum (FHC: +4.3 MeV; RHC: -3.3 MeV - inverting the MC simulation)
    if numode:
        fqmomm = sk_tree.fqmomm + 4.3
    else:
        fqmomm = sk_tree.fqmomm - 3.3
    costheta = dirbeam[0] * dirlep[0][2][0] + dirbeam[1] * dirlep[0][2][1] + dirbeam[2] * dirlep[0][2][2]
    erec_new = erec_ccqe(fqmomm, costheta, M_lep, numode)
    return fqmomm, erec_new

def data_selection(sk_tree, numode):
    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print "Check: number of events in original SK tree = ", n_total

    # Make new filenames for files containing only the numu/nue/cc1pi trees
    out_filenames = dict()
    if numode:
        tag = "t2ksk.14a.neut5.3.2.13a_tuned_v1r1.250ka.fine."
    else:
        tag = "t2ksk.14a.neut5.3.2.13a_tuned_v1r1.m250ka.fine."

    out_filenames["numu"] = tag + skfile[:-9]+"_numuselec.root"
    out_filenames["nue"] = tag + skfile[:-9]+"_nueselec.root"
    out_filenames["cc1pi"] = tag + skfile[:-9]+"_cc1piselec.root"

    print "Making three new files: "
    for key in out_filenames:
        print out_filenames[key]

    out_files = dict()
    out_trees = dict()
    cc_mom = array.array('d', [ 0.0 ])
    cc_erec = array.array('d', [ 0.0 ])
    for key in out_filenames:
        out_files[key] = TFile(out_filenames[key], "RECREATE")
        out_trees[key] = sk_tree.CloneTree(0)
        out_trees[key].Branch("fqmomm_cc", cc_mom, "fqmomm_cc/D")
        out_trees[key].Branch("erec_cc", cc_erec, "erec_cc/D")

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.iclass == 13: # numu fitqun selection
            cc_mom[0], cc_erec[0] = coulomb_correction(sk_tree, numode)
            print sk_tree.fqmomm, sk_tree.erec, cc_mom[0], cc_erec[0]
            if mulike_cut(cc_mom[0]):
                out_trees["numu"].Fill()
                out_trees["numu"].fqmomm = cc_mom
                out_trees["numu"].erec = cc_erec
        elif sk_tree.iclass == 11: # nue+1Rpi0 fitqun selection
            cc_mom[0] = sk_tree.fqmomm
            cc_erec[0] = sk_tree.erec
            if elike_cut(sk_tree):
                out_trees["nue"].Fill()
        elif sk_tree.iclass == 14: # nue cc1pi fitqun selection
            cc_mom[0] = sk_tree.fqmomm
            cc_erec[0] = sk_tree.erec
            if elike_cut(sk_tree):
                out_trees["cc1pi"].Fill()

    for key in out_filenames:
        out_files[key].cd()
        out_trees[key].Write("mtuple")
        out_files[key].Close()

    # Check number of events in output trees
    n_events = dict()
    for key in out_filenames:
        file = TFile(out_filenames[key], "READ")
        tree = file.Get("mtuple")
        n_events[key] = tree.GetEntries()
        file.Close()

    print "Selected events: numu =", n_events["numu"], "nue =", n_events["nue"], "cc1pi =", n_events["cc1pi"]
    print "Sum of selected events =", n_events["numu"] + n_events["nue"] + n_events["cc1pi"]

def mc_selection(sk_tree, numode):
    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print "Check: number of events in original SK tree = ", n_total

    # Make new filenames for files containing only the numu/nue/cc1pi trees
    out_filenames = dict()
    if numode:
        tag = "t2ksk.14a.neut5.3.2.13a_tuned_v1r1.250ka.fine."
    else:
        tag = "t2ksk.14a.neut5.3.2.13a_tuned_v1r1.m250ka.fine."

    out_filenames["numu"] = tag + skfile[:-9]+"_numuselec.root"
    out_filenames["nue"] = tag + skfile[:-9]+"_nueselec.root"
    out_filenames["cc1pi"] = tag + skfile[:-9]+"_cc1piselec.root"

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
            if mulike_cut(sk_tree.fqmomm):
                out_trees["numu"].Fill()
        elif sk_tree.iclass == 11: # nue+1Rpi0 fitqun selection
            if elike_cut(sk_tree):
                out_trees["nue"].Fill()
        elif sk_tree.iclass == 14: # nue cc1pi fitqun selection
            if elike_cut(sk_tree):
                out_trees["cc1pi"].Fill()

    for key in out_filenames:
        out_files[key].cd()
        out_trees[key].Write("mtuple")
        out_files[key].Close()

    # Check number of events in output trees
    n_events = dict()
    for key in out_filenames:
        file = TFile(out_filenames[key], "READ")
        tree = file.Get("mtuple")
        n_events[key] = tree.GetEntries()
        file.Close()

    print "Selected events: numu =", n_events["numu"], "nue =", n_events["nue"], "cc1pi =", n_events["cc1pi"]
    print "Sum of selected events =", n_events["numu"] + n_events["nue"] + n_events["cc1pi"]


if __name__ == "__main__":
    # Get arguments (ie. SK file name)
    parser = argparse.ArgumentParser(description='Splits an SK file into two files, containing the numu- and nue-selected events')
    parser.add_argument("filename", nargs=1)
    parser.add_argument("type", nargs=1)
    args = parser.parse_args()
    skfile = args.filename[0]
    type = args.type[0]

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

    if type == "data":
        data_selection(sk_tree, numode)
    else:
        mc_selection(sk_tree, numode)


