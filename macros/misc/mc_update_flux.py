#!/usr/bin/env python

import ROOT
import array
import argparse

def check(mc_filename):
    mc_file = ROOT.TFile(mc_filename, "READ")
    tree = mc_file.Get("mtuple")
    N = tree.GetEntries()
    tree.GetEntry(0)
    print "Entries", N, "E", tree.pnu[0], "Flx weight", tree.wgtflx, "BeRPA", tree.wgtBeRPA_nom, "fqmome", tree.fqmome
    tree.GetEntry(N / 2)
    print "Entries", N, "E", tree.pnu[0], "Flx weight", tree.wgtflx, "BeRPA", tree.wgtBeRPA_nom, "fqmome", tree.fqmome
    tree.GetEntry(N - 1)
    print "Entries", N, "E", tree.pnu[0], "Flx weight", tree.wgtflx, "BeRPA", tree.wgtBeRPA_nom, "fqmome", tree.fqmome
    mc_file.Close()

def get_flux_weight_histogram(flux_filename, mc_filename):
    if "nueselec" in mc_filename or "cc1piselec" in mc_filename:
        selection = "nue"
    else:
        selection = "numu"
    if "m250ka" in mc_filename:
        selection = selection + "b"
    flux_file = ROOT.TFile(flux_filename, "READ")
    flux_obj = flux_file.Get("enu_sk_tuned13a_" + selection + "_ratio")
    flux_obj.SetDirectory(0)
    flux_file.Close()
    return flux_obj

def update_flux_weights(mc_filename, flux_filename, output_filename):
    flux_obj = get_flux_weight_histogram(flux_filename, mc_filename)
    flux_weight = array.array('d', [ 0.0 ])
    mc_file = ROOT.TFile(mc_filename, "READ")
    tree_in = mc_file.Get("mtuple")
    N = tree_in.GetEntries()
    output_file = ROOT.TFile(output_filename, "RECREATE")
    tree_out = tree_in.CloneTree(0)
    tree_out.SetBranchAddress("wgtflx", flux_weight)

    for e in range(0, N):
        tree_in.GetEntry(e)
        i = flux_obj.FindBin(tree_in.pnu[0])
        flux_weight[0] = flux_obj.GetBinContent(i)
        print "E", tree_in.pnu[0], "Bin", i, "Old", tree_in.wgtflx, "New", flux_weight[0]
        tree_out.Fill()

    output_file.cd()
    tree_out.Write("mtuple")
    output_file.Close()
    mc_file.Close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Update the MC flux weights to external tuning files")
    parser.add_argument("-f", "--filename", type=str, help="The MC mtuple file to update")
    parser.add_argument("--flux", type=str, help="The flux tuning file containing the new weights")
    parser.add_argument("-o", "--output", type=str, help="The output file")
    args = parser.parse_args()
    update_flux_weights(args.filename, args.flux, args.output)
    check(args.filename)
    check(args.output)
  
