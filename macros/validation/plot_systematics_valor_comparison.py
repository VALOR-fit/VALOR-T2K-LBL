#!/usr/bin/python
import sys
import ROOT
import os
import re
import argparse

env_var = "T2KOSCANA"

# See end of file for arguments that need to be parsed
# List of parameters and/or modes are read in from files e.g. valor_other_param_list.arg valor_other_modes.arg
# For a comparison with no modes, supply the *no_modes.arg file ('','')

line_colors = dict()
line_colors["m1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["m3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def configure(compare, selection):
  if compare.lower() == "mach3":
    nom_pattern = "[CURRENT]" + selection + "[MODE]_nom"
    pattern = "[CURRENT]" + selection + "_[PARAM][MODE]_sig_[VAR]"
    # paired objects have valor, then mach3
    variations = [ ("p1","p1"),
                   ("p3","p3"),
                   ("m1","n1"),
                   ("m3","n3") ]
    fhc_current = ""
    rhc_current = "anu_"
  elif compare.lower() == "ptheta":
    nom_pattern = selection + "_[PARAM]_nom_osc" 
    pattern = selection + "_[PARAM]_[VAR]sig_osc"
    # paired objects have valor, then ptheta
    variations = [ ("p1","q1"),
                   ("p3","q3"),
                   ("m1","n1"),
                   ("m3","n3") ]
    fhc_current = ""
    rhc_current = ""
  elif compare.lower() == "valor":
    nom_pattern = "Spectrum_" + selection + "[MODE]_[PARAM]_[VAR]" 
    pattern = "Ratio_" + selection + "[MODE]_[PARAM]_[VAR]"
    # paired objects have valor, then valor
    variations = [ ("p1","p1"),
                   ("p3","p3"),
                   ("m1","m1"),
                   ("m3","m3") ]
    fhc_current = ""
    rhc_current = ""
  else:
    print "Error, unknown comparison fitter '" + compare + "'"
    sys.exit()

  return nom_pattern, pattern, variations, fhc_current, rhc_current

def main(compare, valor_filename, other_filename, valor_selection, other_selection, tag, parameter_filename, mode_filename):
  modes      = list_of_tuples(mode_filename)
  parameters = list_of_tuples(parameter_filename)

  nom_pattern, pattern, variations, fhc_current, rhc_current = configure(compare, other_selection)
  _, valor_pattern, _, _, _ = configure("valor", valor_selection)

  if "RHC" in valor_selection:
    current = rhc_current
  else:
    current = fhc_current
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)

  other_file = ROOT.TFile(other_filename)
  valor_file = ROOT.TFile(valor_filename)

  assert(other_file)
  assert(valor_file)

  if compare.lower() == "mach3":
    mach3_nom_hists         = dict()
    mach3_nom_hists_binno   = dict()

  valor_hists             = dict()
  valor_hists_binno       = dict()

  other_hists             = dict()
  other_ratio_hists       = dict()
  other_hists_binno       = dict()
  other_ratio_hists_binno = dict()

  canvases                = dict()

  pdf_printing = [""]*(len(parameters)*len(modes))
  pdf_printing[ 0] = "("
  pdf_printing[-1] = ")"
  pdf_printing_iter_ereco = iter(pdf_printing)
  pdf_printing_iter_binno = iter(pdf_printing)

  for imode,(valor_mode,other_mode) in enumerate(modes):
    valor_hists             [valor_mode] = dict()
    valor_hists_binno       [valor_mode] = dict()

    other_hists             [valor_mode] = dict()
    other_ratio_hists       [valor_mode] = dict()
    other_hists_binno       [valor_mode] = dict()
    other_ratio_hists_binno [valor_mode] = dict()

    if compare.lower() == "mach3":
      mach3_nom_histname = nom_pattern
      for old,new in [("[CURRENT]"  , current),
                      ("[MODE]"     , other_mode),
                      ("[SELECTION]", other_selection) ]:
        mach3_nom_histname = mach3_nom_histname.replace(old,new)

    if compare.lower() == "mach3":
      mach3_nom_hists[valor_mode] = other_file.Get(mach3_nom_histname)
      if not mach3_nom_hists[valor_mode]:
        print "Failed to get other nominal histogram ",mach3_nom_histname
      print mach3_nom_hists[valor_mode]
      mach3_nom_hists_binno[valor_mode] = hist_to_bin_number(mach3_nom_hists[valor_mode])

    for iparam,(valor_param,other_param) in enumerate(parameters):      
      valor_hists             [valor_mode][valor_param] = dict()
      valor_hists_binno       [valor_mode][valor_param] = dict()
      other_hists             [valor_mode][valor_param] = dict()
      other_ratio_hists       [valor_mode][valor_param] = dict()
      other_hists_binno       [valor_mode][valor_param] = dict()
      other_ratio_hists_binno [valor_mode][valor_param] = dict()

      for valor_var,other_var in variations:
        valor_histname = valor_pattern
        for old,new in [("[PARAM]", valor_param),
                        ("[MODE]" , valor_mode ),
                        ("[VAR]"  , valor_var  )] :
          valor_histname = valor_histname.replace(old,new)
        other_histname = pattern
        if compare.lower() == "valor":
            other_histname = valor_histname
        elif compare.lower() == "mach3":
            for old,new in [("[PARAM]"    , other_param    ),
                            ("[CURRENT]"  , current  ),
                            ("[MODE]"     , other_mode     ),
                            ("[VAR]"      , other_var      ) ]:
              other_histname = other_histname.replace(old,new)
        valor_hists[valor_mode][valor_param][valor_var] = valor_file.Get(valor_histname)

        if not valor_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get VALOR histogram ",valor_histname

        other_hists[valor_mode][valor_param][valor_var] = other_file.Get(other_histname)
        if not other_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get other histogram ",other_histname

        valor_hists[valor_mode][valor_param][valor_var].SetLineWidth(1)
        valor_hists[valor_mode][valor_param][valor_var].SetLineColor(line_colors[valor_var])

        other_hists[valor_mode][valor_param][valor_var].SetName(
            valor_hists[valor_mode][valor_param][valor_var].GetName() + "_mach3")
        other_hists[valor_mode][valor_param][valor_var].SetLineStyle(2)
        other_hists[valor_mode][valor_param][valor_var].SetLineWidth(2)
        other_hists[valor_mode][valor_param][valor_var].SetLineColor(
            valor_hists[valor_mode][valor_param][valor_var].GetLineColor() + 2)

        valor_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_mode][valor_param][valor_var])
        other_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(other_hists[valor_mode][valor_param][valor_var])

        if compare.lower() == "valor":
          other_ratio_hists[valor_mode][valor_param][valor_var] = other_hists[valor_mode][valor_param][valor_var]
          other_ratio_hists_binno[valor_mode][valor_param][valor_var] = other_hists_binno[valor_mode][valor_param][valor_var]
        elif compare.lower() == "mach3":
          other_ratio_hists[valor_mode][valor_param][valor_var] = ratio_hist(other_hists[valor_mode][valor_param][valor_var],mach3_nom_hists[valor_mode])
          other_ratio_hists_binno[valor_mode][valor_param][valor_var] = ratio_hist(other_hists_binno[valor_mode][valor_param][valor_var],mach3_nom_hists_binno[valor_mode])
          set_zero_to_one(other_ratio_hists[valor_mode][valor_param][valor_var])
          set_zero_to_one(other_ratio_hists_binno[valor_mode][valor_param][valor_var])

      hist_limits = []
      for valor_var, other_var in variations:
        for bin_content in (hist_to_list(valor_hists[valor_mode][valor_param][valor_var])+hist_to_list(other_ratio_hists[valor_mode][valor_param][valor_var])):
          if abs(bin_content) > 1E-5:
            hist_limits.append(abs(bin_content-1))

      hist_var = 1.1*max(hist_limits)
      hist_max = hist_var + 1
      hist_min = 1-hist_var

      canv_name = valor_mode+"_"+valor_param

      canvases[canv_name] = ROOT.TCanvas(canv_name,canv_name,1280,720)
      option=""
      for valor_var,other_var in variations:
        valor_hists[valor_mode][valor_param][valor_var].SetTitle(other_mode+" Systematic tweak ratio: "+valor_param+" vs "+other_param)
        valor_hists[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
        valor_hists[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        other_ratio_hists[valor_mode][valor_param][valor_var].Draw(option)
#        other_hists[valor_mode][valor_param][valor_var].Draw(option)

      canvases[canv_name].SetLeftMargin(0.125)   
      canvases[canv_name].SetBottomMargin(0.125)          
      canvases[canv_name].Update()
      canvases[canv_name].Print("syst_var_mode_ereco_"+tag+"_"+valor_selection+".pdf"+next(pdf_printing_iter_ereco))

      canv_name_binno = canv_name+"_binno"

      canvases[canv_name_binno] = ROOT.TCanvas(canv_name_binno,canv_name_binno,1280,720)
      option=""
      for valor_var,other_var in variations:
        valor_hists_binno[valor_mode][valor_param][valor_var].SetTitle(other_mode+" Systematic tweak ratio: "+valor_param+" vs "+other_param)
        valor_hists_binno[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("Bin number")
        valor_hists_binno[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists_binno[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        other_ratio_hists_binno[valor_mode][valor_param][valor_var].Draw(option)

      canvases[canv_name].SetLeftMargin(0.125)   
      canvases[canv_name].SetBottomMargin(0.125)  
      canvases[canv_name_binno].Update()
      canvases[canv_name_binno].Print("syst_var_mode_binno_"+tag+"_"+valor_selection+".pdf"+next(pdf_printing_iter_binno))

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
  parser = argparse.ArgumentParser(description="Get input files from the command line")
  parser.add_argument('-c', '--compare', help="The fitter to compare against. mach3 | ptheta | valor")
  parser.add_argument('-o', '--other-file', dest='other_filename', help="file containing ratio of variations for comparison fitter")
  parser.add_argument('-v', '--valor-file', dest='valor_filename', help="VALOR file containing ratio of variations")
  parser.add_argument('-s', '--other-selection', dest='other_selection', help="selection for comparison fitter. e.g. sknumu | sknue | sknue1pi | sknumubar | sknuebar")
  parser.add_argument('-t', '--tag', dest='tag', help='oscillated or unoscillated')
  parser.add_argument('-r', '--valor-selection', dest='valor_selection', help='valor selection. i.e. 1Rmu | 1Re | MultiRe | RHC1Rmu | RHC1Re')
  parser.add_argument('-p', '--parameter-file', dest='parameter_filename', help="file containing parameter mappings")
  parser.add_argument('-m', '--mode-file', dest='mode_filename', help="file containing mode mappings")
  
  args = parser.parse_args()

  main(args.compare, args.valor_filename, args.other_filename, args.valor_selection, args.other_selection, args.tag,
      args.parameter_filename, args.mode_filename)


