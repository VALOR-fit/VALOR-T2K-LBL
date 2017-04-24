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

def main(mach3_filename,valor_filename,mach3_selection,tag,sample,parameter_filename,mode_filename):

  modes      = list_of_tuples(mode_filename)
  parameters = list_of_tuples(parameter_filename)

  # MaCh3 patterns
  mach3_nom_pattern = "[CURRENT][SELECTION][MODE]_nom"
  mach3_pattern     = "[CURRENT][SELECTION]_[PARAM][MODE]_sig_[VAR]"

  # ptheta patterns
  # mach3_nom_pattern = "[SELECTION]_[PARAM]_nom_osc"
  # mach3_pattern     = "[SELECTION]_[PARAM]_[VAR]sig_osc"
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
                 ("m1","n1"),
                 ("m3","n3") ]

  
  mach3_current = ""  
  if "RHC" in sample:
    mach3_current = "anu_" # only looking at the FHC current samples, and are not labelled in the current round
  else:
    mach3_current = ""
  mach3_selection_1 = mach3_selection # nue CC 1pi sample "sknue1pi" ( "sknumu" "sknue" )
  
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
  
  mach3_nom_hists         = dict() # only valor names will be used as keys
  mach3_nom_hists_binno   = dict() # only valor names will be used as keys

  valor_hists             = dict() # only valor names will be used as keys
  mach3_hists             = dict() # only valor names will be used as keys
  mach3_ratio_hists       = dict() # only valor names will be used as keys
  valor_hists_binno       = dict() # only valor names will be used as keys
  mach3_hists_binno       = dict() # only valor names will be used as keys
  mach3_ratio_hists_binno = dict() # only valor names will be used as keys
  canvases                = dict()

  pdf_printing = [""]*(len(parameters)*len(modes))
  pdf_printing[ 0] = "("
  pdf_printing[-1] = ")"
  pdf_printing_iter_ereco = iter(pdf_printing)
  pdf_printing_iter_binno = iter(pdf_printing)

  for imode,(valor_mode,mach3_mode) in enumerate(modes):
    valor_hists             [valor_mode] = dict()
    mach3_hists             [valor_mode] = dict()
    mach3_ratio_hists       [valor_mode] = dict()
    valor_hists_binno       [valor_mode] = dict()
    mach3_hists_binno       [valor_mode] = dict()
    mach3_ratio_hists_binno [valor_mode] = dict()

    if mach3_mode == "" and not(ptheta_compare) :
      mach3_selection = mach3_selection_1 #"h" + mach3_selection_1
    else :
      mach3_selection = mach3_selection_1

    mach3_nom_histname = mach3_nom_pattern

    for old,new in [("[CURRENT]"  , mach3_current  ),
                    ("[MODE]"     , mach3_mode     ),
                    ("[SELECTION]", mach3_selection) ]:
        mach3_nom_histname = mach3_nom_histname.replace(old,new)

    if not(ptheta_compare):
      mach3_nom_hists[valor_mode] = mach3_file.Get(mach3_nom_histname)
      if not mach3_nom_hists[valor_mode]:
        print "Failed to get MaCh3 nominal histogram ",mach3_nom_histname

      print mach3_nom_hists[valor_mode]

      mach3_nom_hists_binno[valor_mode] = hist_to_bin_number(mach3_nom_hists[valor_mode])


    for iparam,(valor_param,mach3_param) in enumerate(parameters):
      
      valor_hists             [valor_mode][valor_param] = dict()
      mach3_hists             [valor_mode][valor_param] = dict()
      mach3_ratio_hists       [valor_mode][valor_param] = dict()
      valor_hists_binno       [valor_mode][valor_param] = dict()
      mach3_hists_binno       [valor_mode][valor_param] = dict()
      mach3_ratio_hists_binno [valor_mode][valor_param] = dict()

      #if util.IsXSec(util.PhysParamIdFromString(valor_param))

      for valor_var,mach3_var in variations:

        valor_histname = valor_pattern
        for old,new in [("[PARAM]", valor_param),
                        ("[MODE]" , valor_mode ),
                        ("[VAR]"  , valor_var  )] :
          valor_histname = valor_histname.replace(old,new)

        mach3_histname = mach3_pattern
        for old,new in [("[PARAM]"    , mach3_param    ),
                        ("[CURRENT]"  , mach3_current  ),
                        ("[MODE]"     , mach3_mode     ),
                        ("[SELECTION]", mach3_selection_1),
                        ("[VAR]"      , mach3_var      ) ]:
          mach3_histname = mach3_histname.replace(old,new)

        valor_hists[valor_mode][valor_param][valor_var] = valor_file.Get(valor_histname)

        if not valor_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get VALOR histogram ",valor_histname

        mach3_hists[valor_mode][valor_param][valor_var] = mach3_file.Get(mach3_histname)
        if not mach3_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get MaCh3 histogram ",mach3_histname

        valor_hists[valor_mode][valor_param][valor_var].SetLineWidth(1)
        valor_hists[valor_mode][valor_param][valor_var].SetLineColor(line_colors[valor_var])

        mach3_hists[valor_mode][valor_param][valor_var].SetName(valor_hists[valor_mode][valor_param][valor_var].GetName()+"_mach3")
        mach3_hists[valor_mode][valor_param][valor_var].SetLineStyle(2)
        mach3_hists[valor_mode][valor_param][valor_var].SetLineWidth(2)
        mach3_hists[valor_mode][valor_param][valor_var].SetLineColor(valor_hists[valor_mode][valor_param][valor_var].GetLineColor()+2)

        valor_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_mode][valor_param][valor_var])
        mach3_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(mach3_hists[valor_mode][valor_param][valor_var])

        mach3_ratio_hists      [valor_mode][valor_param][valor_var] = ratio_hist(mach3_hists      [valor_mode][valor_param][valor_var],mach3_nom_hists[valor_mode])
        mach3_ratio_hists_binno[valor_mode][valor_param][valor_var] = ratio_hist(mach3_hists_binno[valor_mode][valor_param][valor_var],mach3_nom_hists_binno[valor_mode])

        set_zero_to_one(mach3_ratio_hists      [valor_mode][valor_param][valor_var])
        set_zero_to_one(mach3_ratio_hists_binno[valor_mode][valor_param][valor_var])

      hist_limits = []
      for valor_var,mach3_var in variations:
        for bin_content in (hist_to_list(valor_hists[valor_mode][valor_param][valor_var])+hist_to_list(mach3_ratio_hists[valor_mode][valor_param][valor_var])):
          if abs(bin_content) > 1E-5:
            hist_limits.append(abs(bin_content-1))

      hist_var = 1.1*max(hist_limits)
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
        valor_hists[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        mach3_ratio_hists[valor_mode][valor_param][valor_var].Draw(option)
#        mach3_hists[valor_mode][valor_param][valor_var].Draw(option)

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

if __name__=="__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Get input files from the command line")
  parser.add_argument('-m', '--mach3-file', dest='mach3filename', help="MaCh3 file containing ratio of variations")
  parser.add_argument('-v', '--valor-file', dest='valorfilename', help="VALOR file containing ratio of variations")
  parser.add_argument('-s', '--mach3-selection', dest='mach3selection', help="MaCh3 selection")

  parser.add_argument('-t', '--tag', dest='tag', help='oscillated or unoscillated')
  parser.add_argument('-r', '--valor-sample', dest='sample', help='valor sample')
  parser.add_argument('-p', '--parameter-file', dest='parameter_filename', help="file containing parameter mappings")
  parser.add_argument('-o', '--mode-file', dest='mode_filename', help="file containing mode mappings")
  
  args = parser.parse_args()

  main(args.mach3filename,args.valorfilename,args.mach3selection,args.tag,args.sample,args.parameter_filename,args.mode_filename)


