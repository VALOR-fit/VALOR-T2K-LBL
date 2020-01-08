#!/usr/bin/env python

#
# Plot +/- 1 and 3 sigma variations of selected systematic parameters.
# The plots show (varied spectra) / (nominal spectra).
# Produce the varied spectra that are input into this macro using valor_dune-tweak_pdf_sigma_var.
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# July 2019

import ROOT
import sys
import os
import re

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util  = valor.util
conv  = valor.conv
opt   = valor.Opt.Instance()
env   = valor.Env.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

_name_prefix = "env_sb"
_det     = None
_beam    = None
_selec   = None
_by_mode = False
_x_title = "E_{rec.} (GeV)"

_line_col = dict()
_line_col["p3"] = ROOT.kBlue
_line_col["p1"] = ROOT.kCyan
_line_col["m1"] = ROOT.kViolet
_line_col["m3"] = ROOT.kRed

margin_l = 0.1
margin_r = 0.1
margin_t = 0.1
margin_b = 0.1
title_size = 0.05
label_size = 0.04
title_offset = 0.9

_ratio_diff_warn_threshold_percent = 0.1

# The legacy T2K analysis does not track CC1piC and CC1pi0 modes separately, like we do in the SDK. So when comparing
# systematic variations to the legacy analysis, sum these modes.
_do_sum_cc1pi_modes = False

#_______________________________________________________________________

def main(output_filename, in_file, param_set_str, config, det, beam, selec, comp_file, comp_format, x_max, by_mode, x_title, sum_t2k_cc1pi, modelist_file):
  global _det
  global _beam
  global _selec
  global _by_mode
  global _x_title
  global _do_sum_cc1pi_modes
  _det     = det
  _beam    = beam
  _selec   = selec
  _by_mode = by_mode
  _x_title = x_title

  # TODO: plot twice: 1) with the input kv on x-axis; 2) with the bin number on the x-axis

  elements = valor.CompoundFitElements()
  elements.Read(util.CompleteConfigPath(config), "nomc+nosplines")
  compound_model = elements.FitModel()

  # Get the modes
  if _by_mode:
    if not modelist_file:
      print "ERROR: could not open the input mode list file!"
      sys.exit()
    ml = valor.ModeList(util.CompleteConfigPath(modelist_file))
    read_modes = ml.Modes()

  print "Reading from input VALOR file", in_file
  fin = ROOT.TFile(in_file, "read")
  if not fin.IsOpen():
    print "ERROR: Could not open input file", in_file
    sys.exit()

  draw_comparison = False
  fin_comp = None
  if comp_file != "":
    draw_comparison = True
    fin_comp = ROOT.TFile(comp_file, "READ")

  param_set_temp = list(util.GetParameterSet(param_set_str, True))   # Only get exact matches to any input param names (otherwise problems can be caused)
  param_set = [ipar for ipar in param_set_temp if util.IsSystematic(ipar)]
  if len(param_set) == 0:
    print "ERROR: Got no prams to plot"
    sys.exit()
  print "Plotting systematic variations for the following params:\n", \
      [util.NameRegistry.Instance().ParamAsString(ipar) for ipar in param_set]

  for i in range(compound_model.GetNSubModels()):
    # Check if this is the right sample, if not, then skip it
    selection_model = compound_model.GetSubModel(i)
    selection_name = selection_model.Name()
    if not(selection_model.GetBeam()      == util.BeamConfFromString(_beam) and \
           selection_model.GetDetector()  == util.DetectorFromString(_det)  and \
           selection_model.GetSample()    == util.SampleFromString  (_selec)):
      continue

    if _by_mode and sum_t2k_cc1pi and selection_model.GetDetector() == conv.kdSuperK \
            and comp_format == "t2k_legacy":
      _do_sum_cc1pi_modes = True
      print "Will sum T2K CC1piC and CC1pi0 modes"
    else:
      _do_sum_cc1pi_modes = False

    hnom      = get_nominal(fin, "valor", selection_name)
    hnom_comp = None
    if draw_comparison and not _by_mode:
      hnom_comp = get_nominal(fin_comp, comp_format, selection_name)

    for ip in range(0, len(param_set)):
      p = param_set[ip]
      param_name = util.NameRegistry.Instance().ParamAsStringVerySimple(p)

      modes = ["all"]
      if _by_mode:
        modes = read_modes
        if i == 0 and ip == 0:
          print "Plotting systematic variations for the following modes:\n", \
              [util.ModeAsString(imode) for imode in modes]

      for im in range(0, len(modes)):
        m = modes[im]
        if m == conv.kmUndefined or m == python_undefined or m == python_undefined_short:
          continue
        if _do_sum_cc1pi_modes and util.ReacModeToType(m) == conv.kmtCC1pi0:
          # If comparing to the legacy T2K analysis, skip the CC1pi0 modes, because we've already summed them with CC1piC
          continue

        canv_name = "sigma_var_ratio_" + selection_name + "_" + param_name
        if _by_mode:
          canv_name += "_" + util.ModeAsString(m)
        canv = ROOT.TCanvas(canv_name, canv_name, 800, 600)
        canv.SetMargin(margin_l, margin_r, margin_t, margin_b)
        print "Making plot named", canv_name

        leg = ROOT.TLegend(1.0-margin_r, 0.70, 0.995, 1.0-margin_t)
        hvar_list      = []
        hvar_comp_list = []
        hist_max_var = -999

        # Get the +/- 1 and 3 sigma variations and take the ratio to the nominal
        for sigma_var in ["p3", "p1", "m1", "m3"]:
          # Draw the VALOR hist
          h = get_ratio_nom(fin, sigma_var, param_name, selection_name, m, hnom, False, "")  # Get the VALOR sigma var ratio
          h_comp = None
          draw_opt = ""
          if sigma_var != "p3":
            draw_opt = "SAME"
          else:
            format_drawn_hist(h, x_max)
          h.Draw(draw_opt)
          leg.AddEntry(h, get_label(sigma_var), "l")
          hvar_list.append(h)
          hist_max_var = max(hist_max_var, get_max_var_in_range(h, x_max))

          # Draw the comparison hist
          if draw_comparison:
            h_comp = get_ratio_nom(fin_comp, sigma_var, param_name, selection_name, m, hnom_comp, True, comp_format)  # Get the comparison sigma var ratio
            if h_comp:
              h_comp.Draw("SAME")
              hvar_comp_list.append(h_comp)
              hist_max_var = max(hist_max_var, get_max_var_in_range(h_comp, x_max))
              if sigma_var == "m3":
                leg.AddEntry(h_comp, "comparison", "l")

            # Throw a warning if the ratios have a difference greater than some threshold
            for ibin in range(1, h.GetNbinsX()):
              diff = 100 * ROOT.TMath.Abs(1.0 - (h_comp.GetBinContent(ibin) / h.GetBinContent(ibin)))
              if diff > _ratio_diff_warn_threshold_percent:
                if _by_mode:
                  print "WARNING: For param =", param_name, "; sample =", selection_name, "; mode =", util.ModeAsString(m), "; kinevar val =", \
                    h.GetBinCenter(ibin), "; variation =", sigma_var, "; difference of", diff, \
                    "% in the ratios of input and compaison hist"
                else:
                  print "WARNING: For param =", param_name, "; sample =", selection_name, "; kinevar val =", \
                    h.GetBinCenter(ibin), "; variation =", sigma_var, "; difference of", diff, \
                    "% in the ratios of input and compaison hist"

        hvar_list[0].SetMaximum(1 + 1.1 * hist_max_var)
        hvar_list[0].SetMinimum(1 - 1.1 * hist_max_var)

        title_text = selection_name + "  " + param_name
        if _by_mode:
          if not (_do_sum_cc1pi_modes and util.ReacModeToType(m) == conv.kmtCC1piC):
            title_text = selection_name + "  " + util.ModeAsString(m) + "  " + param_name
          else:
            title_text = selection_name + "  " + util.ModeAsString(m).replace("CC1piC", "CC1pi(C+0)") + "  " + param_name
        title_box = make_title_box(title_text)
        title_box.Draw("SAME")
        leg.Draw("SAME")

        # Save the canvas to the output file
        canv.Update()
        save_name = output_filename + "_" + selection_name + "_" + param_set_str + ".pdf"
        if _by_mode:
          save_name = output_filename + "_mode_" + "_" + selection_name + "_" + param_set_str + ".pdf"
        start_pdf = False
        end_pdf   = False
        if _by_mode:
          if len(modes) > 1:
            if im == 0:
              start_pdf = True
            elif im == len(modes) -1:
              end_pdf   = True
        else:
          if len(param_set) > 1:
            if ip == 0:
              start_pdf = True
            elif ip == len(param_set) -1:
              end_pdf   = True

        if start_pdf:
          save_name += "("
        elif end_pdf:
          save_name += ")"
        canv.SaveAs(save_name, "RECREATE")
        del canv

        # Free up memory
        for h in hvar_list:
          del h
        for h in hvar_comp_list:
          del h
        del title_box
    del hnom


  print "Done!"

#_______________________________________________________________________
def get_ratio_nom(fin, sigma_var, param_name, selection_name, m, hnom, is_comparison, comparison_type):
  hname = ""
  h = None
  if not is_comparison:
    hname = _name_prefix + "_" + selection_name + "_" + param_name + "_" + sigma_var
    if _by_mode:
      hname = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(m) + "_" + param_name + "_" + sigma_var
    h = fin.Get(hname)

    if _do_sum_cc1pi_modes and util.ReacModeToType(m) == conv.kmtCC1piC:
      # If comparing to the legacy T2K analysis, sum the PDFs of the CC1pi0 and CC1piC modes
      hname = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(m).replace("CC1piC", "CC1pi0") + "_" + param_name + "_" + sigma_var
      hCC1pi0 = fin.Get(hname)
      if not hCC1pi0:
        print "WARNING: Could not get hist named [", hname, "] from the input file"
      h.Add(hCC1pi0)

    # If we're considering the contribution from each mode, get the nominal hist here instead of using the input nominal hist
    if _by_mode:
      hnom = get_nominal(fin, "valor", selection_name, m)
    h.Divide(hnom)
  else:
    if comparison_type == "valor":
      hname = _name_prefix + "_" + selection_name + "_" + param_name + "_" + sigma_var
      if _by_mode:
        hname = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(m) + "_" + param_name + "_" + sigma_var
      h = fin.Get(hname)

      if _do_sum_cc1pi_modes and util.ReacModeToType(m) == conv.kmtCC1piC:
        # If comparing to the legacy T2K analysis, sum the PDFs of the CC1pi0 and CC1piC modes
        hname = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(m).replace("CC1piC", "CC1pi0") + "_" + param_name + "_" + sigma_var
        hCC1pi0 = fin.Get(hname)
        if not hCC1pi0:
          print "WARNING: Could not get hist named [", hname, "] from the input file"
        h.Add(hCC1pi0)

      # If we're considering the contribution from each mode, get the nominal hist here instead of using the input nominal hist
      if _by_mode:
        hnom = get_nominal(fin, "valor", selection_name, m)
      h.Divide(hnom)
    elif comparison_type == "t2k_legacy":
      hname = get_t2k_legacy_hist_name(param_name, sigma_var, m)
      h = fin.Get(hname)  # This is already a ratio, so don't divide it.
    else:
      print "ERROR: Unsupported comparison type [", comparison_type, "]"
      sys.exit()

  if not h:
    print "WARNING: Could not get hist named [", hname, "] from the input file"
    return None
  ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope

  # Set empty bins to 1
  for ibin in range(1, h.GetNbinsX()+1):
    if h.GetBinContent(ibin) < conv.kASmallNum:
      h.SetBinContent(ibin, 1)

  # Formatting
  if not is_comparison:
    h.SetLineColor(_line_col[sigma_var])
    h.SetLineWidth(1)
  else:
    h.SetLineColor(_line_col[sigma_var] +2)
    h.SetLineStyle(2)
    h.SetLineWidth(2)

  return h

#_______________________________________________________________________
def get_nominal(fin, format, selection_name, mode = conv.kmUndefined):
  if format == "valor":
    nom_name = _name_prefix + "_" + selection_name + "_nom"
    if _by_mode:
      nom_name = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(mode) + "_nom"
    hnom = fin.Get(nom_name)

    if _do_sum_cc1pi_modes and util.ReacModeToType(mode) == conv.kmtCC1piC:
      # If comparing to the legacy T2K analysis, sum the PDFs of the CC1pi0 and CC1piC modes
      hname = _name_prefix + "_" + selection_name + "_" + util.ModeAsString(mode).replace("CC1piC", "CC1pi0") + "_nom"
      hCC1pi0 = fin.Get(hname)
      if not hCC1pi0:
        print "WARNING: Could not get hist named [", hname, "] from the input file"
      hnom.Add(hCC1pi0)
  elif format == "t2k_legacy":
    # The legacy T2K syst var hists have already had the ratio taken, so do nothing here
    return None
  else:
    print "ERROR: Unsupported format [", format, "]"
    sys.exit()

  ROOT.SetOwnership(hnom, False)       # This stops the object being deleted from memory when it goes out of scope
  return hnom

#_______________________________________________________________________
def get_label(sigma_var):
  if sigma_var   == "p3":
    return "+3 #sigma"
  elif sigma_var == "p1":
    return "+1 #sigma"
  elif sigma_var == "m1":
    return "-1 #sigma"
  elif sigma_var == "m3":
    return "-3 #sigma"
  else:
    print "ERROR: could not convert variation [", sigma_var, "] to label"
    sys.exit()

#_______________________________________________________________________
def get_t2k_legacy_hist_name(param_name, sigma_var, m):
  # Open up a file containing a parameter name map
  parameters = list_of_tuples(util.CompleteConfigPath("maps/t2xk_sdk_param_name_map.arg"))
  comp_param_name = ""
  for iparam,(sdk_param,legacy_param) in enumerate(parameters):
    #print "considering SDK name", sdk_param, "; t2xk name", legacy_param
    if sdk_param == param_name:
      #print "   matches!"
      comp_param_name = legacy_param
      break

  # Get the T2K legacy sample
  comp_sample = ""
  is_RHC = False
  if "rhc" in util.ToLowerCase(_beam):
    is_RHC = True
  if "1re" in util.ToLowerCase(_selec):
    if is_RHC:
      comp_sample = "RHC1Re"
    else:
      comp_sample = "1Re"
  elif "1rmu" in util.ToLowerCase(_selec):
    if is_RHC:
      comp_sample = "RHC1Rmu"
    else:
      comp_sample = "1Rmu"
  elif "multire" in util.ToLowerCase(_selec):
    if is_RHC:
      print "ERROR: no RHC MultiRe samples!"
      sys.exit()
    else:
      comp_sample = "MultiRe"

  comp_mode_name = ""
  if _by_mode:
    sdk_mode_str = util.ToLowerCase(util.ModeAsString(m))
    comp_mode_name = sdk_mode_str  .replace("ccmec", "mec")
    comp_mode_name = comp_mode_name.replace("cc1pic", "cc1pi")
    comp_mode_name = comp_mode_name.replace("cc1pi0", "cc1pi")
    comp_mode_name = comp_mode_name.replace("nc1pic", "nc1piC")
    comp_mode_name = comp_mode_name.replace("numutonue", "oscnue")
    comp_mode_name = comp_mode_name.replace("numubartonuebar", "oscnuebar")

  hname = "Ratio_" + comp_sample + "_" + comp_param_name + "_" + sigma_var
  if _by_mode:
    hname = "Ratio_" + comp_sample + "_" + comp_mode_name +"_" + comp_param_name + "_" + sigma_var
  return hname

#_______________________________________________________________________
def list_of_tuples(filename):
  tuplelist = []
  with open(filename) as file:
    for line in file:
      filtered = re.sub('[\"\'\(]|(\)\s*\,)|[\)]|(#.*$)|^\s*$|\s', '', line)
      if re.match(r'^\s*$', filtered):
        continue
      tuplelist.append(tuple(filtered.split(',')))

  return tuplelist

#_______________________________________________________________________
def format_drawn_hist(h, x_max):
  h.GetXaxis().SetTitle(_x_title)
  h.GetYaxis().SetTitle("Ratio to Nominal")
  h.GetYaxis().SetLabelSize(label_size)
  h.GetXaxis().SetLabelSize(label_size)
  h.GetYaxis().SetTitleSize(title_size)
  h.GetXaxis().SetTitleSize(title_size)
  h.GetYaxis().SetTitleOffset(title_offset)
  h.GetXaxis().SetTitleOffset(title_offset)
  h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(1), x_max)

#_______________________________________________________________________
def get_max_var_in_range(h, x_max):
  max_var = -999
  for i in range(1, h.GetNbinsX()):
    if h.GetXaxis().GetBinUpEdge(i) <= x_max:
      max_var = max(max_var, ROOT.TMath.Abs(h.GetBinContent(i)-1))  # -1 to make the plot symmetric around ratio = 1
  return max_var

#_______________________________________________________________________
def make_title_box(title):
  title_box = ROOT.TPaveText(0.3, 1.0-margin_t+0.005, 0.7, 1.0, "NDC")
  title_box.AddText(title)
  title_box.SetTextSize(0.05)
  title_box.SetFillColor(0)
  title_box.SetFillStyle(0)
  title_box.SetLineStyle(0)
  title_box.SetLineColor(0)
  ROOT.SetOwnership(title_box, False)       # This stops the object being deleted from memory when it goes out of scope
  return title_box

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Plot single-sample +/- 1 and 3 sigma variations produced by tweak_pdf_sigma_var.cxx")
  
  parser.add_option('-o', '--out-file'     , dest='output_filename', default=""     , help="Name of the output image file")
  parser.add_option('-i', '--in-file'      , dest='in_file'        , default=""     , help="Input file name produced by tweak_pdf_sigma_var.cxx")
  parser.add_option('-p', '--params'       , dest='params'         , default=""     , help="Name of the tweak parameter set")
  parser.add_option('-c', '--config'       , dest='config'         , default=""     , help="Name of the analysis config file")
  parser.add_option('-d', '--det'          , dest='det'            , default=""     , help="The detector for this sample")
  parser.add_option('-b', '--beam'         , dest='beam'           , default=""     , help="The beam mode for this sample")
  parser.add_option('-s', '--selec'        , dest='selec'          , default=""     , help="The selection for this sample")
  parser.add_option('-m', '--by-mode'      , dest='by_mode'        , default=False  , action="store_true", help="Display the variations broken down by mode")
  parser.add_option(      '--x-title'      , dest='x_title'        , default="E_{rec.} (GeV)", help="The x-axis title of the 1D PDFs")
  parser.add_option('-x', '--x-max'        , dest='x_max'          , type=float     , default=3.0, help="Maximum x-axis value")
  parser.add_option(      '--comp-file'    , dest='comp_file'      , default=""     , help="Name of file to compare to")
  parser.add_option(      '--comp-format'  , dest='comp_format'    , default=""     , help="Choose from: \'valor\', \'t2k_legacy\'")
  parser.add_option(      '--sum-t2k-cc1pi', dest='sum_t2k_cc1pi'  , default=False  , action="store_true", help="Sum T2k CC1piC and CC1pi0 modes when comparing to T2K legacy")
  parser.add_option(      '--mode-list'    , dest='modelist_file'  , default=""     , help="File containing the list of used reaction modes")

  o,other = parser.parse_args()

  main(o.output_filename, o.in_file, o.params, o.config, o.det, o.beam, o.selec, o.comp_file, o.comp_format, o.x_max, o.by_mode, o.x_title, o.sum_t2k_cc1pi, o.modelist_file)

  if other:
    print "Warning - ignored arguments: ", other
