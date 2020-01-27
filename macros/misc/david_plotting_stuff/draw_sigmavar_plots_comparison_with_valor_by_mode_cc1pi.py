#!/usr/bin/python
import sys
import ROOT
import os
import re

env_var = "T2KOSCANA"

# See end of file for arguments that need to be parsed
# List of parameters and/or modes are read in from files e.g. valor_mach3_param_list.arg valor_mach3_modes.arg
# For a comparison with no modes, supply the *no_modes.arg file ('','')

line_colors = dict()
line_colors["m1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["m3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def main(mach3_filename,valor_filename,tag,sample,parameter_filename,mode_filename):

  xmax = 1.25
  if "1Rmu" in sample:
    xmax = 3.0

  modes      = list_of_tuples(mode_filename)
  parameters = list_of_tuples(parameter_filename)

  # MaCh3 patterns
  mach3_pattern     = "Ratio_"+sample+"[MODE]_[PARAM]_[VAR]"
  ptheta_compare = False

  #  valor_nom_pattern = "Spectrum_"+sample+"_nominal"
#  valor_nom_pattern     = "Spectrum_"+sample+"[MODE]_[PARAM]_[VAR]" 
  valor_pattern     = "Ratio_"+sample+"[MODE]_[PARAM]_[VAR]" 
#  valor_pattern     = valor_nom_pattern
  
  # paired objects have valor, then ptheta
  # variations = [ ("p1","q1"),
  #                ("p3","q3"),
  #                ("m1","n1"),
  #                ("m3","n3") ]

  # paired objects have valor, then mach3
  variations = [ ("p1","p1"),
                 ("p3","p3"),
                 ("m1","m1"),
                 ("m3","m3") ]
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)

  mach3_file = ROOT.TFile(mach3_filename)
  valor_file = ROOT.TFile(valor_filename)

  assert(mach3_file)
  assert(valor_file)

  # valor_nom_histname = valor_nom_pattern
  
  # valor_nom_hist = valor_file.Get(valor_nom_histname)
  # if not valor_nom_hist:
  #   print "Failed to get VALOR nominal histogram ", valor_nom_histname

  # valor_nom_hist_binno = hist_to_bin_number(valor_nom_hist)

  valor_hists             = dict() # only valor names will be used as keys
  mach3_ratio_hists       = dict() # only valor names will be used as keys
  valor_hists_binno       = dict() # only valor names will be used as keys
  mach3_ratio_hists_binno = dict() # only valor names will be used as keys
  canvases                = dict()

  pdf_printing = [""]*(len(parameters)*len(modes))
  pdf_printing[ 0] = "("
  pdf_printing[-1] = ")"
  pdf_printing_iter_ereco = iter(pdf_printing)
  pdf_printing_iter_binno = iter(pdf_printing)

  for imode,(valor_mode,mach3_mode) in enumerate(modes):
    valor_hists             [valor_mode] = dict()
    mach3_ratio_hists       [valor_mode] = dict()
    valor_hists_binno       [valor_mode] = dict()
    mach3_ratio_hists_binno [valor_mode] = dict()

    for iparam,(valor_param,mach3_param) in enumerate(parameters):
      
      valor_hists             [valor_mode][valor_param] = dict()
      mach3_ratio_hists       [valor_mode][valor_param] = dict()
      valor_hists_binno       [valor_mode][valor_param] = dict()
      mach3_ratio_hists_binno [valor_mode][valor_param] = dict()

      #if util.IsXSec(util.PhysParamIdFromString(valor_param))

      hist_max_var = -999
      hist_min_var = 999
      for valor_var,mach3_var in variations:

        valor_histname = valor_pattern
        for old,new in [("[PARAM]", valor_param),
                        ("[MODE]" , valor_mode ),
                        ("[VAR]"  , valor_var  )] :
          valor_histname = valor_histname.replace(old,new)

        valor_hists[valor_mode][valor_param][valor_var] = valor_file.Get(valor_histname)
        mach3_ratio_hists[valor_mode][valor_param][valor_var] = mach3_file.Get(valor_histname)

        if not valor_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get VALOR histogram ",valor_histname
        if not mach3_ratio_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get MaCh3 histogram ",valor_histname

        valor_hists[valor_mode][valor_param][valor_var].SetLineWidth(1)
        valor_hists[valor_mode][valor_param][valor_var].SetLineColor(line_colors[valor_var])

        mach3_ratio_hists[valor_mode][valor_param][valor_var].SetName(valor_hists[valor_mode][valor_param][valor_var].GetName()+"_mach3")
        mach3_ratio_hists[valor_mode][valor_param][valor_var].SetLineStyle(2)
        mach3_ratio_hists[valor_mode][valor_param][valor_var].SetLineWidth(2)
        mach3_ratio_hists[valor_mode][valor_param][valor_var].SetLineColor(valor_hists[valor_mode][valor_param][valor_var].GetLineColor()+2)

        valor_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_mode][valor_param][valor_var])
        mach3_ratio_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(mach3_ratio_hists[valor_mode][valor_param][valor_var])

        set_zero_to_one(mach3_ratio_hists      [valor_mode][valor_param][valor_var])
        set_zero_to_one(mach3_ratio_hists_binno[valor_mode][valor_param][valor_var])

        hist_max_var = max(hist_max_var, get_max_var_in_range(valor_hists[valor_mode][valor_param][valor_var], xmax))
        hist_max_var = max(hist_max_var, get_max_var_in_range(mach3_ratio_hists[valor_mode][valor_param][valor_var], xmax))
        hist_min_var = min(hist_min_var, get_min_var_in_range(valor_hists[valor_mode][valor_param][valor_var], xmax))
        hist_min_var = min(hist_min_var, get_min_var_in_range(mach3_ratio_hists[valor_mode][valor_param][valor_var], xmax))

      hist_max_var = max(hist_max_var, ROOT.TMath.Abs(hist_min_var))
      hist_var = 1.1*hist_max_var
      hist_max = hist_var + 1
      hist_min = 1-hist_var

      canv_name = valor_mode+"_"+valor_param

      canvases[canv_name] = ROOT.TCanvas(canv_name,canv_name,1280,720)
      option=""
      for valor_var,mach3_var in variations:
        valor_hists[valor_mode][valor_param][valor_var].SetTitle(mach3_mode+" Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
        valor_hists[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
        valor_hists[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists[valor_mode][valor_param][valor_var].GetXaxis().SetRangeUser(valor_hists[valor_mode][valor_param][valor_var].GetXaxis().GetBinLowEdge(1), xmax)
        valor_hists[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        mach3_ratio_hists[valor_mode][valor_param][valor_var].Draw(option)

      title_box = make_title_box(canv_name)
      title_box.Draw("SAME")

      canvases[canv_name].SetLeftMargin(0.125)   
      canvases[canv_name].SetBottomMargin(0.125)          
      canvases[canv_name].Update()
      canvases[canv_name].Print("syst_var_mode_ereco_"+tag+"_"+sample+".pdf"+next(pdf_printing_iter_ereco))

      canv_name_binno = canv_name+"_binno"

      canvases[canv_name_binno] = ROOT.TCanvas(canv_name_binno,canv_name_binno,1280,720)
      option=""
      for valor_var,mach3_var in variations:
        valor_hists_binno[valor_mode][valor_param][valor_var].SetTitle(mach3_mode+" Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
        valor_hists_binno[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("Bin number")
        valor_hists_binno[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists_binno[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        mach3_ratio_hists_binno[valor_mode][valor_param][valor_var].Draw(option)

      canvases[canv_name].SetLeftMargin(0.125)   
      canvases[canv_name].SetBottomMargin(0.125)  
      canvases[canv_name_binno].Update()
      canvases[canv_name_binno].Print("syst_var_mode_binno_"+tag+"_"+sample+".pdf"+next(pdf_printing_iter_binno))

  print "Done!"
  return

#_______________________________________________________________________

def hist_to_bin_number(h):
  nbins = h.GetXaxis().GetNbins()
  hnew = ROOT.TH1D(h.GetName()+"_binno",h.GetTitle()+" by bins",nbins,0,nbins)

  hnew.SetLineColor(h.GetLineColor())
  hnew.SetLineStyle(h.GetLineStyle())
  hnew.SetLineWidth(h.GetLineWidth())

  for ibin in xrange(1,nbins+1):
    hnew.Fill(ibin-0.5, h.GetBinContent(ibin))

  hnew.GetXaxis().SetTitle("Bin number")

  return hnew

def ratio_hist(hvar,hnom):
  nbins = hvar.GetXaxis().GetNbins()
  hratio = ROOT.TH1D(hvar)
  hratio.SetName (hvar.GetName ()+"_ratio")
  hratio.SetTitle(hvar.GetTitle()+" Ratio")
  hratio.GetYaxis().SetTitle("Ratio")
  hratio.Divide(hnom)
  return hratio

def set_zero_to_one(h):
  nbins = h.GetXaxis().GetNbins()
  for ibin in xrange(1,nbins+1):
    if abs(h.GetBinContent(ibin) < 1E-10):
      h.SetBinContent(ibin,1)
  return

def hist_to_list(h):
  nbins = h.GetXaxis().GetNbins()
  return [h.GetBinContent(ibin) for ibin in xrange(1,nbins+1)]

def list_of_tuples(filename):
  tuplelist = []
  with open(filename) as file:
    for line in file:
      filtered = re.sub('[\"\'\(]|(\)\s*\,)|[\)]|(#.*$)|^\s*$|\s', '', line)
      if re.match(r'^\s*$', filtered):
        continue
      tuplelist.append(tuple(filtered.split(',')))

  return tuplelist

def get_max_var_in_range(h, x_max):
  max_var = -999
  for i in range(1, h.GetNbinsX()):
    if h.GetXaxis().GetBinUpEdge(i) <= x_max:
      max_var = max(max_var, ROOT.TMath.Abs(h.GetBinContent(i)-1))  # -1 to make the plot symmetric around ratio = 1
  return max_var

def get_min_var_in_range(h, x_max):
  max_var = 999
  for i in range(1, h.GetNbinsX()):
    if h.GetXaxis().GetBinUpEdge(i) <= x_max:
      max_var = min(max_var, ROOT.TMath.Abs(h.GetBinContent(i)-1))  # -1 to make the plot symmetric around ratio = 1
  return max_var

def make_title_box(title):
  title_box = ROOT.TPaveText(0.3, 1.0-0.1+0.005, 0.7, 1.0, "NDC")
  title_box.AddText(title)
  title_box.SetTextSize(0.05)
  title_box.SetFillColor(0)
  title_box.SetFillStyle(0)
  title_box.SetLineStyle(0)
  title_box.SetLineColor(0)
  ROOT.SetOwnership(title_box, False)       # This stops the object being deleted from memory when it goes out of scope
  return title_box

if __name__=="__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Get input files from the command line")
  parser.add_argument('-m', '--mach3-file', dest='mach3filename', help="MaCh3 file containing ratio of variations")
  parser.add_argument('-v', '--valor-file', dest='valorfilename', help="VALOR file containing ratio of variations")

  parser.add_argument('-t', '--tag', dest='tag', help='oscillated or unoscillated')
  parser.add_argument('-r', '--valor-sample', dest='sample', help='valor sample')
  parser.add_argument('-p', '--parameter-file', dest='parameter_filename', help="file containing parameter mappings")
  parser.add_argument('-o', '--mode-file', dest='mode_filename', help="file containing mode mappings")
  
  args = parser.parse_args()

  main(args.mach3filename,args.valorfilename,args.tag,args.sample,args.parameter_filename,args.mode_filename)


