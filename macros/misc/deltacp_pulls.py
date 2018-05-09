import sys
from ROOT import *
import argparse
import array

gROOT.SetBatch(1)

def calculate_pulls(filename):
    f = TFile(filename, "READ")
    tree = f.Get("fitnt")
    hpulls = TH1D("deltacp_pulls", "Pulls", 100, -5, 5)
    hpulls.GetXaxis().SetTitle("Pull")
    hpulls.GetXaxis().SetTitle("Fraction of fits")
    for event in tree:
        bestfit_nh = event.bestfit0_nh
        trueval = event.truepar0
        hesse_err = event.bestfiterrSym0_nh
        pull = (bestfit_nh - trueval) / hesse_err
        hpulls.Fill(pull)
    hpulls.Scale(1.0 / hpulls.Integral())
    canvas = TCanvas("canvas", "canvas", 1024, 768)
    canvas.SetTopMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetBottomMargin(0.15)
    canvas.SetRightMargin(0.2)
    hpulls.Draw()
    canvas.SaveAs("pulls.pdf")
    f.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates pulls on deltacp')
    parser.add_argument("-f", "--filename", dest="filename", type=str, help="ROOT filename")
    args = parser.parse_args()

    calculate_pulls(args.filename)

