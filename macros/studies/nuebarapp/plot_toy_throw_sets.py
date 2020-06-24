#!/usr/bin/env python

import ROOT
from copy import deepcopy
import collections
ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

def main(infile_name, ntoys, split):

    infile = ROOT.TFile(infile_name, "read")
    if ntoys % split != 0:
        print "The total of", ntoys, "toys does not divide evenly into", split, "parts. Exiting."
        util.Die(1)
    ntoys_split = ntoys / split

    # Initialise the hists that we'll fill
    hToyHisto_dcp   = []
    hToyHisto_Mix13 = []
    hToyHisto_Mix23 = []
    hToyHisto_Mix12 = []
    hToyHisto_Dm23  = []
    hToyHisto_Dm12  = []
    for i in range(0, split):
        h = infile.Get("hToyHisto_dcp_nh")
        hToyHisto_dcp  .append(ROOT.TH1D("dcp_nh_"+str(i)  , "dcp_nh_"+str(i)  , h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))
        h = infile.Get("hToyHisto_Mix13_nh")
        hToyHisto_Mix13.append(ROOT.TH1D("Mix13_nh_"+str(i), "Mix13_nh_"+str(i), h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))
        h = infile.Get("hToyHisto_Mix23_nh")
        hToyHisto_Mix23.append(ROOT.TH1D("Mix23_nh_"+str(i), "Mix23_nh_"+str(i), h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))
        h = infile.Get("hToyHisto_Mix12_nh")
        hToyHisto_Mix12.append(ROOT.TH1D("Mix12_nh_"+str(i), "Mix12_nh_"+str(i), h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))
        h = infile.Get("hToyHisto_Dm23_nh")
        hToyHisto_Dm23 .append(ROOT.TH1D("Dm23_nh_"+str(i) , "Dm23_nh_"+str(i) , h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))
        h = infile.Get("hToyHisto_Dm12_nh")
        hToyHisto_Dm12 .append(ROOT.TH1D("Dm12_nh_"+str(i) , "Dm12_nh_"+str(i) , h.GetNbinsX(), h.GetBinLowEdge(1), h.GetBinLowEdge(h.GetNbinsX()) + h.GetBinWidth(1)))

    # Read the PhysParamLists in the toy file to fill the hists
    for itoy in range(0, ntoys):
        isplit = -1
        for try_split in range(0, split):
            if itoy < (try_split + 1) * ntoys_split:
                isplit = try_split
                break
        params = infile.Get("iparm_" + str(itoy))
        hToyHisto_dcp  [isplit].Fill(params.Value(conv.kPhypCPV))
        hToyHisto_Mix13[isplit].Fill(params.Value(conv.kPhypMixing13))
        hToyHisto_Mix23[isplit].Fill(params.Value(conv.kPhypMixing23))
        hToyHisto_Mix12[isplit].Fill(params.Value(conv.kPhypMixing12))
        hToyHisto_Dm23 [isplit].Fill(params.Value(conv.kPhypDm23sq))
        hToyHisto_Dm12 [isplit].Fill(params.Value(conv.kPhypDm12sq))

    # Plot the hists
    draw_hists(hToyHisto_dcp  , "#delta_{CP}"          , split, True , False)
    draw_hists(hToyHisto_Mix13, "sin^{2} (#theta_{13})", split, True , False)
    draw_hists(hToyHisto_Mix23, "sin^{2} (#theta_{23})", split, False, False)
    draw_hists(hToyHisto_Mix12, "sin^{2} (#theta_{12})", split, False, False)
    draw_hists(hToyHisto_Dm23 , "#Delta m^{2}_{23}"    , split, False, False)
    draw_hists(hToyHisto_Dm12 , "#Delta m^{2}_{12}"    , split, False, True )

    print "Done!"


def draw_hists(hists, label, split, first, last):
    canv = ROOT.TCanvas(label, label, 800, 600)
    for i in range(0, split):
        draw_opt = "same"
        if i == 0:
            hists[i].GetXaxis().SetTitle(label)
            hists[i].GetYaxis().SetTitle("Entries")
            draw_opt = ""
        hists[i].SetLineColor(i+1)
        hists[i].Draw(draw_opt)
    canv.Update()
    outfile_name = "split_throws.pdf"
    if first and not last:
        outfile_name += "("
    elif last and not first:
        outfile_name += ")"
    canv.SaveAs(outfile_name)


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Combine many fake-data toy root files into one")
    parser.add_option('--infile', dest='infile', default="", help='Name of the input toy file')
    parser.add_option('--ntoys', dest='ntoys', default=20000, type=int, help='The number of toys in the input file')
    parser.add_option('--split', dest='split', default=4, type=int, help='Split the input toys into this number of equally sized groups')

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o, other = parser.parse_args()

    main(o.infile, o.ntoys, o.split)

    if other:
        print "Warning - ignored arguments: ",other


