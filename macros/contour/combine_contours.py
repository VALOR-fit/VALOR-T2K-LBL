#!/usr/bin/env python

import sys
import ROOT

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



#_______________________________________________________________________________________________________________________
def main(in_filenames, out_filename):

    if len(in_filenames) < 2:
        print "FATAL: Only", len(in_filenames), "input files. Need at least 2. Exiting."
        sys.exit()

    chisq_list  = []
    is_1D = True

    # Read the input files
    for filename in in_filenames:
        file = ROOT.TFile(filename, "READ")
        if file.GetListOfKeys().Contains("h_deltachisq"):
            # Is 1D
            is_1D = True
            hchisq = file.Get("h_chisq")
            hchisq.SetDirectory(0)
            chisq_list.append(hchisq)
            file.Close()
        elif file.GetListOfKeys().Contains("dchisq_map_2D"):
            # Is 2D
            is_1D = False
            hchisq = file.Get("chisq_map_2D")
            hchisq.SetDirectory(0)
            chisq_list.append(hchisq)
            file.Close()
        else:
            print "Error: No dchisq surface found"
            sys.exit()
        print "Read chisq hists from file", file, "; is 1D?", is_1D

    # Make the (empty) output hists
    hchisq_out  = None
    hdchisq_out = None
    if is_1D:
        hchisq_out  = ROOT.TH1D("h_chisq", "h_chisq",
                                chisq_list[0].GetNbinsX(),
                                chisq_list[0].GetXaxis().GetBinLowEdge(1),
                                chisq_list[0].GetXaxis().GetBinUpEdge(chisq_list[0].GetNbinsX()))
        hdchisq_out = ROOT.TH1D("h_deltachisq", "h_deltachisq",
                                chisq_list[0].GetNbinsX(),
                                chisq_list[0].GetXaxis().GetBinLowEdge(1),
                                chisq_list[0].GetXaxis().GetBinUpEdge(chisq_list[0].GetNbinsX()))
    else:
        hchisq_out  = ROOT.TH2D("chisq_map_2D", "chisq_map_2D",
                                chisq_list[0].GetNbinsX(),
                                chisq_list[0].GetXaxis().GetBinLowEdge(1),
                                chisq_list[0].GetXaxis().GetBinUpEdge(chisq_list[0].GetNbinsX()),
                                chisq_list[0].GetNbinsY(),
                                chisq_list[0].GetYaxis().GetBinLowEdge(1),
                                chisq_list[0].GetYaxis().GetBinUpEdge(chisq_list[0].GetNbinsY()))
        hdchisq_out = ROOT.TH2D("dchisq_map_2D", "dchisq_map_2D",
                                chisq_list[0].GetNbinsX(),
                                chisq_list[0].GetXaxis().GetBinLowEdge(1),
                                chisq_list[0].GetXaxis().GetBinUpEdge(chisq_list[0].GetNbinsX()),
                                chisq_list[0].GetNbinsY(),
                                chisq_list[0].GetYaxis().GetBinLowEdge(1),
                                chisq_list[0].GetYaxis().GetBinUpEdge(chisq_list[0].GetNbinsY()))

    # Fill the output chisq hist
    if is_1D:
        for ibinx in range(1, chisq_list[0].GetNbinsX() + 1):
            marg_like_bin = 0
            for ihist in range(0, len(chisq_list)):
                chisq = chisq_list[ihist].GetBinContent(ibinx)
                marg_like_bin += chisq_to_likelihood(chisq)
            marg_like_bin /= float(len(chisq_list))
            hchisq_out.SetBinContent(ibinx, likelihood_to_chisq(marg_like_bin))
    else:
        for ibinx in range(1, chisq_list[0].GetNbinsX() + 1):
            for ibiny in range(1, chisq_list[0].GetNbinsY() + 1):
                marg_like_bin = 0
                for ihist in range(0, len(chisq_list)):
                    chisq = chisq_list[ihist].GetBinContent(ibinx, ibiny)
                    marg_like_bin += chisq_to_likelihood(chisq)
                marg_like_bin /= float(len(chisq_list))
                hchisq_out.SetBinContent(ibinx, ibiny, likelihood_to_chisq(marg_like_bin))

    # Fill the output delta chisq hist
    chisq_min = hchisq_out.GetMinimum()
    if is_1D:
        for ibinx in range(1, chisq_list[0].GetNbinsX() + 1):
            hdchisq_out.SetBinContent(ibinx, hchisq_out.GetBinContent(ibinx) - chisq_min)
    else:
        for ibinx in range(1, chisq_list[0].GetNbinsX() + 1):
            for ibiny in range(1, chisq_list[0].GetNbinsY() + 1):
                hdchisq_out.SetBinContent(ibinx, ibiny, hchisq_out.GetBinContent(ibinx, ibiny) - chisq_min)

    # Save the combined contours to the output file
    outfile = ROOT.TFile(out_filename, "RECREATE")
    save_name_chisq  = None
    save_name_dchisq = None
    if is_1D:
        save_name_chisq  = "h_chisq"
        save_name_dchisq = "h_deltachisq"
    else:
        save_name_chisq  = "chisq_map_2D"
        save_name_dchisq = "dchisq_map_2D"
    hchisq_out .Write(save_name_chisq )
    hdchisq_out.Write(save_name_dchisq)
    outfile.Close()

    print "Done!"

#_______________________________________________________________________________________________________________________
def chisq_to_likelihood(chisq):
    return ROOT.TMath.Exp(-0.5*chisq)

#_______________________________________________________________________________________________________________________
def likelihood_to_chisq(likelihood):
    return -2.0*ROOT.TMath.Log(likelihood)

#_______________________________________________________________________________________________________________________
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filenames", dest="in_filenames", nargs="+", help="Input filenames for contours to combine.")
    parser.add_argument("-o", "--output"   , dest="out_filename"           , help="Output ROOT filename for combined contour.")
    args = parser.parse_args()

    main(args.in_filenames, args.out_filename)
