import os, sys
import array
import numpy as np
from ROOT import *
import argparse

M_e = 0.510998
M_m = 105.65836
M_n = 939.565346
M_p = 938.272013
M_d = 1232.0
Eb = 27.0

def elike_cut(evis, erec):
    return evis > 100 and erec < 1250

def mulike_cut(p):
    return p > 200

# Current CCQE Erec calculation uses FHC equation in RHC mode as well
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

def lepton_angle(dirbeam, dirlep, is_electron):
    if is_electron:
        return dirbeam[0] * dirlep[0][1][0] + dirbeam[1] * dirlep[0][1][1] + dirbeam[2] * dirlep[0][1][2]
    else:
        return dirbeam[0] * dirlep[0][2][0] + dirbeam[1] * dirlep[0][2][1] + dirbeam[2] * dirlep[0][2][2]

def coulomb_correction(p, numode=true):
    if numode:
        return p + 4.3
    else:
        return p - 3.3

def select(sk_tree, out_filename):
    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print "Events in original SK tree:", n_total
    out_file = TFile(out_filename, "RECREATE")
    sk_tree_out = sk_tree.CloneTree(0)
    cc_mom = array.array('f', [ 0.0 ])
    cc_erec = array.array('d', [ 0.0 ])
    sk_tree_out.SetBranchAddress("fqmomm", cc_mom)
    sk_tree_out.SetBranchAddress("erec", cc_erec)
    for i in range(n_total):
        sk_tree.GetEntry(i)
        if sk_tree.iclass == 13: # numu fitqun selection
            dirbeam = np.array([0.669764, -0.742179, 0.024223])
            dirlep = np.reshape(sk_tree.fq1rdir, (-1, 7, 3))
            costheta = lepton_angle(dirbeam, dirlep, False)
            cc_mom[0] = coulomb_correction(sk_tree.fqmomm, sk_tree.horn_current > 0)
            cc_erec[0] = erec_ccqe(cc_mom[0], costheta, M_m, True)
            if mulike_cut(cc_mom[0]):
                sk_tree_out.Fill()
        elif sk_tree.iclass == 11: # nue+1Rpi0 fitqun selection
            if elike_cut(sk_tree.fqevis, sk_tree.erec):
                cc_mom[0] = sk_tree.fqmomm
                cc_erec[0] = sk_tree.erec
                sk_tree_out.Fill()
        elif sk_tree.iclass == 14: # nue cc1pi fitqun selection
            if elike_cut(sk_tree.fqevis, sk_tree.erec):
                cc_mom[0] = sk_tree.fqmomm
                cc_erec[0] = sk_tree.erec
                sk_tree_out.Fill()
    out_file.cd()
    sk_tree_out.Write("minituple")
    out_file.Close()

    # Check number of events in output trees
    n_events = dict()
    file = TFile(out_filename, "READ")
    tree = file.Get("minituple")
    n_events = tree.GetEntries()
    print "Events in coulomb corrected tree:", n_events
    file.Close()

if __name__ == "__main__":
    # Get arguments (ie. SK file name)
    parser = argparse.ArgumentParser(description='Splits an SK file into two files, containing the numu- and nue-selected events')
    parser.add_argument('-f', '--in-filename', dest='in_filename', help="Input SK minituple.")
    parser.add_argument('-o', '--out-filename', dest='out_filename', help="Output SK minituple.")
    args = parser.parse_args()

    # Get SK tree
    f_sk = TFile(args.in_filename, "OPEN")
    if not f_sk:
        print "Error: could not open file ", skfile
        sys.exit()
    sk_tree = f_sk.Get("minituple")
    if not sk_tree:
        sk_tree = f_sk.Get("h1")
    select(sk_tree, args.out_filename)
    f_sk.Close()

