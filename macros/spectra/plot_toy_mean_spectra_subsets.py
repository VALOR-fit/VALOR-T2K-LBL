#!/usr/bin/env python

#
# Take in 4 different mean and RMS hists (generated from 4 independent sets of toys) and compare them all.
# Also take in a file containing all 4 sets combined, so that they can be compared to the mean
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# October 2019

import ROOT
import sys
import os

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

#_______________________________________________________________________

def main(subset_file_names, combined_file_name, out_filename):

  hist_suffix = "_mean_PDF"
  samples =  [
    "superk__JPARC_FHC__1rmu",
    "superk__JPARC_RHC__1rmu",
    "superk__JPARC_FHC__1re",
    "superk__JPARC_RHC__1re",
    "superk__JPARC_FHC__multire"
  ]

  combined_mean_hists     = dict()  # Mean
  combined_RMS_up_hists   = dict()  # Mean + RMS
  combined_RMS_down_hists = dict()  # Mean - RMS
  subset_mean_hists       = dict()  # Mean
  subset_RMS_up_hists     = dict()  # Mean + RMS
  subset_RMS_down_hists   = dict()  # Mean - RMS
  for sample in samples:
    subset_mean_hists    [sample] = []
    subset_RMS_up_hists  [sample] = []
    subset_RMS_down_hists[sample] = []

  # Open the input files and extract the relevant hists
  combined_file = ROOT.TFile(combined_file_name, "read")
  if not combined_file:
    print "ERROR: Could not open the input file:", combined_file_name
    sys.exit()
  for sample in samples:
    hname_mean     = sample + hist_suffix
    hname_RMS_up   = sample + hist_suffix + "_err_up"
    hname_RMS_down = sample + hist_suffix + "_err_down"
    h_mean     = combined_file.Get(hname_mean)
    h_RMS_up   = combined_file.Get(hname_RMS_up)
    h_RMS_down = combined_file.Get(hname_RMS_down)
    h_mean    .SetDirectory(0)
    h_RMS_up  .SetDirectory(0)
    h_RMS_down.SetDirectory(0)
    combined_mean_hists    [sample] = h_mean
    combined_RMS_up_hists  [sample] = h_RMS_up
    combined_RMS_down_hists[sample] = h_RMS_down
  combined_file.Close()

  for filename in subset_file_names:
    subset_file = ROOT.TFile(filename, "read")
    if not subset_file:
      print "ERROR: Could not open the input file:", filename
      sys.exit()
    print "Reading from subset file:", filename

    for sample in samples:
      hname_mean     = sample + hist_suffix
      hname_RMS_up   = sample + hist_suffix + "_err_up"
      hname_RMS_down = sample + hist_suffix + "_err_down"
      h_mean     = subset_file.Get(hname_mean)
      h_RMS_up   = subset_file.Get(hname_RMS_up)
      h_RMS_down = subset_file.Get(hname_RMS_down)
      h_mean    .SetDirectory(0)
      h_RMS_up  .SetDirectory(0)
      h_RMS_down.SetDirectory(0)
      subset_mean_hists    [sample].append(h_mean)
      subset_RMS_up_hists  [sample].append(h_RMS_up)
      subset_RMS_down_hists[sample].append(h_RMS_down)
    subset_file.Close()

  # Set styles and draw
  for sample in samples:
    canv = ROOT.TCanvas(sample, sample, 700, 700)
    leg = ROOT.TLegend(0.70, 0.80, 0.90, 0.95)

    # Get x and y max values
    max_x = 3.
    if not "1rmu" in sample:
      max_x = 1.25

    max_y = get_max_in_range_1D(combined_RMS_up_hists[sample], max_x)
    for isubset in range(len(subset_RMS_up_hists[sample])):
      max_y = max( max_y, get_max_in_range_1D(subset_RMS_up_hists  [sample][isubset], max_x) )
    max_y *= 1.2

    # Draw the frame
    this_pad = canv.cd(1)
    this_pad.SetLeftMargin(0.2)
    pdf_frame = this_pad.DrawFrame(0, 0, max_x, max_y)
    pdf_frame.SetDirectory(0)
    pdf_frame.GetXaxis().SetTitle("#nu Reconstructed Energy (GeV)")
    pdf_frame.GetYaxis().SetTitle("Number of Events")
    pdf_frame.Draw()

    # Set styles
    combined_mean_hists    [sample].SetLineWidth(3)
    combined_mean_hists    [sample].SetLineStyle(1)
    combined_mean_hists    [sample].SetLineColor(1)
    combined_RMS_up_hists  [sample].SetFillColor(ROOT.kRed)
    combined_RMS_up_hists  [sample].SetFillStyle(1001)
    combined_RMS_up_hists  [sample].SetLineColor(1)
    combined_RMS_up_hists  [sample].SetLineWidth(3)
    combined_RMS_down_hists[sample].SetFillColor(10)
    combined_RMS_down_hists[sample].SetFillStyle(1001)
    combined_RMS_down_hists[sample].SetLineColor(1)
    combined_RMS_down_hists[sample].SetLineWidth(3)
    subset_colours = [20, 30, 40, 46]
    #subset_styles  = [2 , 3 , 4 , 5 ]
    subset_styles  = [2 , 2 , 2 , 2 ]
    for isubset in range(len(subset_RMS_up_hists[sample])):
      subset_mean_hists    [sample][isubset].SetFillStyle(0)
      subset_RMS_up_hists  [sample][isubset].SetFillStyle(0)
      subset_RMS_down_hists[sample][isubset].SetFillStyle(0)
      subset_mean_hists    [sample][isubset].SetLineWidth(1)
      subset_RMS_up_hists  [sample][isubset].SetLineWidth(1)
      subset_RMS_down_hists[sample][isubset].SetLineWidth(1)
      subset_mean_hists    [sample][isubset].SetLineStyle(subset_styles [isubset])
      subset_RMS_up_hists  [sample][isubset].SetLineStyle(subset_styles [isubset])
      subset_RMS_down_hists[sample][isubset].SetLineStyle(subset_styles [isubset])
      subset_mean_hists    [sample][isubset].SetLineColor(subset_colours[isubset])
      subset_RMS_up_hists  [sample][isubset].SetLineColor(subset_colours[isubset])
      subset_RMS_down_hists[sample][isubset].SetLineColor(subset_colours[isubset])

    # Draw
    combined_RMS_up_hists  [sample].Draw("hist,same")
    combined_RMS_down_hists[sample].Draw("hist,same")
    combined_mean_hists    [sample].Draw("hist,same")
    leg.AddEntry(combined_mean_hists  [sample], "Mean, 40k toys"     , "l")
    leg.AddEntry(combined_RMS_up_hists[sample], "RMS error, 40k toys", "f")

    for isubset in range(len(subset_RMS_up_hists[sample])):
      subset_mean_hists    [sample][isubset].Draw("hist,same")
      subset_RMS_up_hists  [sample][isubset].Draw("hist,same")
      subset_RMS_down_hists[sample][isubset].Draw("hist,same")
      leg.AddEntry(subset_mean_hists[sample][isubset], "Set "+str(isubset+1)+", 10k toys", "l")
    leg.Draw()

    # Draw title
    title = ROOT.TPaveText(0, max_y, max_x, max_y * 1.07 )
    title.SetShadowColor(0)
    title.SetFillColor(0)
    title.SetFillStyle(0)
    title.SetLineColor(0)
    title.AddText(sample)
    title.SetTextSize(0.035)
    title.Draw("NB,same")

    # Save canvas to the output file
    save_name = out_filename
    if sample == samples[0]:
      save_name += "("
    elif sample == samples[-1]:
      save_name += ")"
    canv.SaveAs(save_name)


  print "Done!"


#_______________________________________________________________________
def get_max_in_range_1D(hist, xmax):
  max_n = 0.
  for xbin in range(1, hist.GetNbinsX()+1):
    if hist.GetBinCenter(xbin) < xmax:
      max_n = max(max_n, hist.GetBinContent(xbin))
  return max_n

#_______________________________________________________________________
if __name__=="__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Index splines")
  parser.add_argument("-f", "--subset-files" , dest="subset_files" , nargs="+", help="Input files to combine.")
  parser.add_argument("-c", "--combined-file", dest="combined_file"           , help="Output combined filename")
  parser.add_argument("-o", "--output"       , dest="out_filename"            , help="Output PDF image filename")
  args = parser.parse_args()

  if len(args.subset_files) != 4:
    print "ERROR: This macro is set up for 4 subset files, but", len(args.subset_files), "were input. Exiting."
    sys.exit()

  main(args.subset_files, args.combined_file, args.out_filename)
