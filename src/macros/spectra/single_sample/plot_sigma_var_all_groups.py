#!/usr/bin/env python

# Plot +/- 1 and 3 sigma variations of selected systematic parameters.
# The plots show (varied spectra) / (nominal spectra).
#
# Input one ROOT file per group, each containing TH1D's of either the varied and
# nominal spectra or the ratio of the two (depending on the options set below).
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# Feburary 2020

import ROOT
import sys
import re
import copy
from copy import deepcopy

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

# Names used internally within this macro
_valid_formats    = ["VALOR", "MaCh3", "p-theta", "VALOR_legacy"]
_valid_selections = ["FHC1Rmu", "FHC1Re", "RHC1Rmu", "RHC1Re", "FHCMultiRe"]
_variations       = ["p3", "p1", "m1", "m3"]

#===================================================================================
# STUFF THAT CHANGES DEPENDING ON THE HIST NAMES AND FORMATS FOR EACH GROUP (START)
#===================================================================================
# In the following dictionaries, the keys are the names used internally within this
# macro (and should stay as they are), and the values are the group-specific names
# used to retrieve the stored hists from the input files (and can be altered as needed).

# The naming pattern used for the nominal PDFs
_hist_nom_name_patterns = dict()
_hist_nom_name_patterns["VALOR"]        = "env_sb_[SAMPLE][MODE]_nom"
_hist_nom_name_patterns["MaCh3"]        = "[SAMPLE][MODE]_nom"
_hist_nom_name_patterns["p-theta"]      = "[SAMPLE]_[PARAM][MODE]_nom_[OSC]"
_hist_nom_name_patterns["VALOR_legacy"] = None

# The naming pattern used for the PDFs produced by varying one of the parameters and
# keeping the others at their nominal values
_hist_tweaked_name_patterns = dict()
_hist_tweaked_name_patterns["VALOR"]        = "env_sb_[SAMPLE][MODE]_[PARAM]_[VAR]"
_hist_tweaked_name_patterns["MaCh3"]        = "[SAMPLE]_[PARAM][MODE]_sig_[VAR]"
_hist_tweaked_name_patterns["p-theta"]      = "[SAMPLE]_[PARAM][MODE]_[VAR]sig_[OSC]"
_hist_tweaked_name_patterns["VALOR_legacy"] = "Ratio_[SAMPLE][MODE]_[PARAM]_[VAR]"

# If the hists are stored as the ratio of the nominal PDF to the tweaked PDF, set this to True
_ratio_already_taken = dict()
_ratio_already_taken["VALOR"]        = False
_ratio_already_taken["MaCh3"]        = False
_ratio_already_taken["p-theta"]      = False
_ratio_already_taken["VALOR_legacy"] = True

# Pattern used in the hist name to indicate whether the hist corresponds to the -3, -1, +1, or +3 sigma tweak.
# In the dict below "m3" is the -3 sigma tweak and "p1" is the +1 sigma tweak etc.
_var_names = dict()
_var_names["VALOR"]         = { "m3":"m3", "m1":"m1", "p1":"p1", "p3":"p3" }
_var_names["MaCh3"]         = { "m3":"n3", "m1":"n1", "p1":"p1", "p3":"p3" }
_var_names["p-theta"]       = { "m3":"m3", "m1":"m1", "p1":"p1", "p3":"p3" }
_var_names["VALOR_legacy"]  = { "m3":"m3", "m1":"m1", "p1":"p1", "p3":"p3" }

# Beam mode + event topology
_selecs_all = dict()
_selecs_all["VALOR"       ] = { "FHC1Rmu":"superk__JPARC_FHC__1rmu", "FHC1Re":"superk__JPARC_FHC__1re", "RHC1Rmu":"superk__JPARC_RHC__1rmu", "RHC1Re":"superk__JPARC_RHC__1re", "FHCMultiRe":"superk__JPARC_FHC__multire" }
_selecs_all["MaCh3"       ] = { "FHC1Rmu":"sknumu"                 , "FHC1Re":"sknue"                 , "RHC1Rmu":"anu_sknumubar"          , "RHC1Re":"anu_sknuebar"          , "FHCMultiRe":"sknue1pi"                   }
_selecs_all["p-theta"     ] = { "FHC1Rmu":"numu1R"                 , "FHC1Re":"nue1R"                 , "RHC1Rmu":"numubar1R"              , "RHC1Re":"nuebar1R"              , "FHCMultiRe":"nue1RD"                     }
_selecs_all["VALOR_legacy"] = { "FHC1Rmu":"1Rmu"                   , "FHC1Re":"1Re"                   , "RHC1Rmu":"RHC1Rmu"                , "RHC1Re":"RHC1Re"                , "FHCMultiRe":"MultiRe"                    }

#==================================================================================
# STUFF THAT CHANGES DEPENDING ON THE HIST NAMES AND FORMATS FOR EACH GROUP (END)
#==================================================================================

_formats     = None
_input_files = dict()
_selecs      = dict()
_by_mode     = False
_mode_map    = None
_n_modes     = None
_param_map   = None

_line_col = dict()
_line_col["p3"] = ROOT.kBlue
_line_col["p1"] = ROOT.kCyan
_line_col["m1"] = ROOT.kViolet
_line_col["m3"] = ROOT.kRed

_margin_l     = 0.1
_margin_r     = 0.1
_margin_t     = 0.1
_margin_b     = 0.1
_title_size   = 0.05
_label_size   = 0.04
_title_offset = 0.9
_a_small_num  = 1E-10

# Optionally throw a warning when the ratio hists from each group differ by more than this percentage
_ratio_diff_warn_threshold_percent = 0.1

# The legacy VALOR T2K analysis does not track CC1piC and CC1pi0 modes separately, like we do in the VALOR SDK.
# So when comparing systematic variations to the legacy analysis, sum these modes.
_do_sum_cc1pi_modes = False

# The legacy VALOR T2K analysis does not track CCmisc, CCMpi and CCDIS modes separately, like we do in the VALOR SDK.
# So when comparing systematic variations to the legacy analysis, sum these modes into CCoth.
_do_sum_ccoth_modes = False

#_______________________________________________________________________
def main(output_filename, in_files, formats, selec, params_file, modes_file, x_max, suppress_warn, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth):
  setup(in_files, formats, selec, params_file, modes_file, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth)

  for ip in range(0, len(_param_map)):
    param_name_dict = get_param_name_dict(ip)
    for im in range(0, _n_modes):
      mode_name_dict = None
      if _by_mode:
        mode_name_dict = get_mode_name_dict(im)
        if(should_skip_mode(mode_name_dict)):
          continue

      canv_name = "sigma_var_ratio_" + _selecs[_formats[0]] + "_" + param_name_dict[_formats[0]]
      if _by_mode:
        canv_name += "_" + mode_name_dict[_formats[0]]
      canv_ereco = ROOT.TCanvas(canv_name + "_Ereco", canv_name + "_Ereco", 800, 600)
      canv_binno = ROOT.TCanvas(canv_name + "_binno", canv_name + "_binno", 800, 600)
      canv_ereco.SetMargin(_margin_l, _margin_r, _margin_t, _margin_b)
      canv_binno.SetMargin(_margin_l, _margin_r, _margin_t, _margin_b)
      print "Making plot named", canv_name

      leg = ROOT.TLegend(1.0-_margin_r, 0.65, 0.995, 1.0-_margin_t)
      leg.AddEntry(None, "Variations", "")
      h_list       = []
      h_list_binno = []
      hist_max_var = -999

      # Get the +/- 1 and 3 sigma variations and take the ratio to the nominal
      for ivar, sigma_var in enumerate(_variations):
        comparison_dict = dict()
        for iformat, f in enumerate(_formats):
          h = get_ratio_hist(sigma_var, param_name_dict, mode_name_dict, iformat)
          if not h:
            # Skip any hists that we couldn't read from the input file
            continue
          h.SetName(canv_name + "_" + sigma_var + "_" + f)
          h_binno = hist_to_bin_number(h)   # Also plot with the bin number on the x-axis

          draw_opt = "HIST"
          if ivar != 0 or iformat != 0:
            draw_opt = "HIST,SAME"
          else:
            format_drawn_hist(h, x_max, False)
            format_drawn_hist(h_binno, h.FindBin(x_max) -1, True)
          canv_ereco.cd()
          h.Draw(draw_opt)
          canv_binno.cd()
          h_binno.Draw(draw_opt)

          if iformat == 0:
            # Label for each variation
            leg.AddEntry(h, get_label(sigma_var), "l")
          if ivar == len(_variations) -1:
            # Label for each group
            if iformat == 0:
              leg.AddEntry(None, "", "")  # Blank line
              leg.AddEntry(None, "Groups", "")
            h_copy = ROOT.TH1D(h)
            h_copy.SetLineColor(1)  # Keep the group labels in black, so the focus is just on the line style
            leg.AddEntry(h_copy, f, "l")
            h_list.append(h_copy)

          h_list.append(h)
          h_list_binno.append(h_binno)
          comparison_dict[f] = h
          hist_max_var = max(hist_max_var, get_max_var_in_range(h, x_max))

        # Throw a warning if any of the ratio hists have a difference greater than some threshold
        if not suppress_warn:
          compare_ratio_hists(param_name_dict, mode_name_dict, sigma_var, comparison_dict)
        del comparison_dict

      # Set the y-axis range to be symmetric around ratio = 1
      h_list      [0].SetMaximum(1 + 1.1 * hist_max_var)
      h_list      [0].SetMinimum(1 - 1.1 * hist_max_var)
      h_list_binno[0].SetMaximum(1 + 1.1 * hist_max_var)
      h_list_binno[0].SetMinimum(1 - 1.1 * hist_max_var)

      title_box = make_title_box(get_title_box_text(param_name_dict, mode_name_dict))
      canv_ereco.cd()
      title_box.Draw("SAME")
      leg.Draw("SAME")
      canv_binno.cd()
      title_box.Draw("SAME")
      leg.Draw("SAME")

      canv_ereco.cd()
      canv_ereco.Update()
      canv_ereco.SaveAs(get_output_name(output_filename, ip, im, False), "RECREATE")
      canv_binno.cd()
      canv_binno.Update()
      canv_binno.SaveAs(get_output_name(output_filename, ip, im, True), "RECREATE")

      # Free up memory
      del canv_ereco
      del canv_binno
      del h_list[:]
      del h_list
      del h_list_binno[:]
      del h_list_binno
      del title_box
      del leg

  print "Done!"

#_______________________________________________________________________
def setup(in_files, formats, selec, params_file, modes_file, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth):
  global _input_files
  global _formats
  global _selecs
  global _by_mode
  global _mode_map
  global _n_modes
  global _param_map
  global _is_unosc
  global _do_sum_cc1pi_modes
  global _do_sum_ccoth_modes
  _by_mode  = modes_file != ""
  _formats  = formats
  _is_unosc = is_unosc

  # Check that the inputs are valid
  if(len(in_files) != len(formats)):
    print "ERROR: The number of input files and formats should be the same, but instead we have", \
      len(in_files), "input files,", len(formats), "formats."
    sys.exit()

  for format_name in formats:
    if format_name not in _valid_formats:
      print "ERROR: input format \'{}\'".format(format_name), "is invalid. Choose from:", _valid_formats
      sys.exit()

  if selec not in _valid_selections:
    print "ERROR: input selection \'{}\'".format(selec), "is invalid. Choose from:", _valid_selections
    sys.exit()

  for i,f in enumerate(_formats):
    _selecs[f] = _selecs_all[f][selec]
    print "\nFor", f, ", using the selection", _selecs[f]

    _input_files[f] = ROOT.TFile(in_files[i], "read")
    if not _input_files[f].IsOpen():
      print "ERROR: Could not open the input file", in_files[i]
      sys.exit()

  # Get the list of parameters
  if not params_file:
    print "ERROR: Could not open the input parameter map file!"
    sys.exit()
  _param_map = list_of_tuples(params_file)
  _param_map = sort_map(_param_map, _formats)
  if len(_param_map) == 0:
    print "ERROR: Got no prams to plot"
    sys.exit()

  print "\nPlotting systematic variations for the following params:\n", \
    [param_tuple[0] for param_tuple in _param_map]

  # Get the list of modes, if required
  _n_modes = 1
  if _by_mode:
    if not modes_file:
      print "ERROR: Could not open the input mode map file!"
      sys.exit()
    _mode_map = list_of_tuples(modes_file)
    _mode_map = sort_map(_mode_map, _formats)
    if len(_mode_map) == 0:
      print "ERROR: Got no modes to plot"
      sys.exit()

    _n_modes = len(_mode_map)

    print "\nPlotting systematic variations for the following modes:\n", \
      [mode_tuple[0] for mode_tuple in _mode_map]

  # In VALOR, the CC1piC and CC1pi0 treated separately, so their PDFs must be summed when comparing to
  # groups who keep track of the sum of these modes
  if _by_mode and sum_t2k_cc1pi and "VALOR" in _formats:
    _do_sum_cc1pi_modes = True
    print "\nWill sum T2K CC1piC and CC1pi0 modes"

  if _by_mode and sum_t2k_ccoth and "VALOR" in _formats and "VALOR_legacy" in _formats:
    _do_sum_ccoth_modes = True
    print "\nWill sum T2K CCmisc, CCMpi and CCDIS modes"

#_______________________________________________________________________
def get_param_name_dict(ip):
  param_name = dict()
  for iformat, f in enumerate(_formats):
    param_name[f] = _param_map[ip][iformat]
  if len(param_name) != len(_formats):
    print "ERROR: The number of entries in each row of the parameter map file should equal the length of the --formats argument. Exiting."
    sys.exit()
  return param_name

#_______________________________________________________________________
def get_mode_name_dict(im):
  mode_name = dict()
  for iformat, f in enumerate(_formats):
    mode_name[f] = _mode_map[im][iformat]
  if len(mode_name) != len(_formats):
    print "ERROR: The number of entries in each row of the mode map file should equal the length of the --formats argument. Exiting."
    sys.exit()
  return mode_name

#_______________________________________________________________________
def should_skip_mode(mode_name_dict):
  # If comparing to the legacy T2K analysis, skip the CC1pi0 modes, because we've already summed them with CC1piC
  if _do_sum_cc1pi_modes and _formats[0] == "VALOR" and "CC1pi0" in mode_name_dict["VALOR"]:
    return True

  # If comparing to the legacy T2K analysis, skip the kmtCCDIS and kmtCCmisc modes, because we've already summed them with CCmisc
  if _do_sum_ccoth_modes and _formats[0] == "VALOR" and "CCMpi" in mode_name_dict["VALOR"]:
    return True

  return False

#_______________________________________________________________________
def get_title_box_text(param_name_dict, mode_name_dict):
  print_selec_name = _selecs[_formats[0]]
  print_param_name = param_name_dict[_formats[0]]
  title_text = print_selec_name + "  " + print_param_name
  if _by_mode:
    print_mode_name  = mode_name_dict[_formats[0]]
    if not (_do_sum_cc1pi_modes and "CC1piC" in print_mode_name) or \
            (_do_sum_ccoth_modes and "CCmisc" in print_mode_name):
      title_text = print_selec_name + "  " + print_mode_name + "  " + print_param_name

    elif (_do_sum_ccoth_modes and "CCmisc" in print_mode_name):
      title_text = print_selec_name + "  " + print_mode_name.replace("CCmisc", "CCoth (CCmisc+CCDIS+CCMpi)") + "  " + print_param_name

    elif (_do_sum_cc1pi_modes and "CC1piC" in print_mode_name):
      title_text = print_selec_name + "  " + print_mode_name.replace("CC1piC", "CC1pi(C+0)") + "  " + print_param_name

    else:
      title_text = print_selec_name + "  " + print_param_name
  return title_text

#_______________________________________________________________________
def get_ratio_hist(sigma_var, param_name_dict, mode_name_dict, iformat):
  format = _formats[iformat]
  hratio = None
  if not _ratio_already_taken[format]:
    # Make the ratio hist: tweaked / nominal
    hnom   = get_nominal(mode_name_dict, format)
    htweak = get_tweaked(sigma_var, param_name_dict, mode_name_dict, format)
    if not hnom or not htweak:
      return None
    hratio = ROOT.TH1D(htweak)
    hratio.Divide(hnom)
  else:
    # The ratio hist is already in the input file, so just retrieve it
    hratio = get_tweaked(sigma_var, param_name_dict, mode_name_dict, format)
    if not hratio:
      return None

  set_empty_bins_to_one(hratio)

  # Set the line style and colour differently for each group
  col_offset = (2 * iformat)
  # The colour offset values in TColor range from -10 to +4 (inclusive), so keep within range
  if col_offset > 4:
    col_offset *= -1
  if col_offset < -10:
    print "WARNING: Setting the TColor offset to", col_offset, ", which may not be a valid colour."
  hratio.SetLineColor(_line_col[sigma_var] + col_offset)
  hratio.SetLineStyle(1 + iformat)
  hratio.SetLineWidth(1 + iformat)

  ROOT.SetOwnership(hratio, False)       # This stops the object being deleted from memory when it goes out of scope
  return hratio

#_______________________________________________________________________
def set_empty_bins_to_one(h):
  for ibin in range(1, h.GetNbinsX()+1):
    if h.GetBinContent(ibin) < _a_small_num:
      h.SetBinContent(ibin, 1)

#_______________________________________________________________________
def get_nominal(mode_name_dict, format):
  return get_hist(False, None, None, mode_name_dict, format)

#_______________________________________________________________________
def get_tweaked(sigma_var, param_name_dict, mode_name_dict, format):
  return get_hist(True, sigma_var, param_name_dict, mode_name_dict, format)

#_______________________________________________________________________
def get_hist(get_tweaked_hist, sigma_var, param_name_dict, mode_name_dict, format):
  hname = None

  if get_tweaked_hist:
    hname = _hist_tweaked_name_patterns[format]
  else:
    hname = _hist_nom_name_patterns[format]

  hname = hname.replace("[SAMPLE]", _selecs[format])

  if get_tweaked_hist:
    hname = hname.replace("[PARAM]" , param_name_dict[format])
    hname = hname.replace("[VAR]"   , _var_names[format][sigma_var])

  if mode_name_dict:
    hname = hname.replace("[MODE]", "_" + mode_name_dict[format])
  else:
    hname = hname.replace("[MODE]", "")

  if _is_unosc:
    hname = hname.replace("[OSC]", "unosc")
  else:
    hname = hname.replace("[OSC]", "osc")

  h = None
  if format == "VALOR":
    h = get_valor_hist(hname)
  elif format == "MaCh3":
    h = get_mach3_hist(hname)
  elif format == "p-theta":
    h = get_ptheta_hist(hname)
  elif format == "VALOR_legacy":
    h = get_t2k_legacy_hist(hname)
  else:
    print "ERROR: Unrecognised format \'{}\'".format(format)
    sys.exit()

  if h:
    ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope
  return h

#_______________________________________________________________________
def get_valor_hist(hname):
  # Put in any special cases related to getting VALOR hists here
  fin = _input_files["VALOR"]
  h = fin.Get(hname)
  if not h:
    print "WARNING: Could not get hist named [", hname, "] from the input file,", fin.GetName()

  if _do_sum_cc1pi_modes and "CC1piC" in hname:
    # If comparing to the legacy T2K analysis, sum the PDFs of the CC1pi0 and CC1piC modes
    hname_add = hname.replace("CC1piC", "CC1pi0")
    hCC1pi0 = fin.Get(hname_add)
    if not hCC1pi0:
      print "WARNING: Could not get hist named [", hname_add, "] from the input file,", fin.GetName()
    h.Add(hCC1pi0)

  elif _do_sum_ccoth_modes and "CCmisc" in hname:
    # If comparing to the legacy T2K analysis, sum the PDFs of the CCMpi, CCMisc and CCDIS modes
    hname_add = hname.replace("CCmisc", "CCDIS") + "_nom"
    hCCDIS = fin.Get(hname_add)
    if not hCCDIS:
      print "WARNING: Could not get hist named [", hname_add, "] from the input file,", fin.GetName()
    h.Add(hCCDIS)

    hname_add = hname.replace("CCmisc", "CCMpi") + "_nom"
    hCCMpi = fin.Get(hname_add)
    if not hCCMpi:
      print "WARNING: Could not get hist named [", hname_add, "] from the input file,", fin.GetName()
    h.Add(hCCMpi)

  if h:
    ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope
  return h

#_______________________________________________________________________
def get_mach3_hist(hname):
  # Put in any special cases related to getting MaCh3 hists here
  fin = _input_files["MaCh3"]
  h = fin.Get(hname)
  if not h:
    print "WARNING: Could not get hist named [", hname, "] from the input file,", fin.GetName()
  else:
    ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope
  return h

#_______________________________________________________________________
def get_ptheta_hist(hname):
  # Put in any special cases related to getting p-theta hists here
  fin = _input_files["p-theta"]
  h = fin.Get(hname)
  if not h:
    print "WARNING: Could not get hist named [", hname, "] from the input file,", fin.GetName()
  else:
    ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope
  return h

#_______________________________________________________________________
def get_t2k_legacy_hist(hname):
  # Put in any special cases related to getting VALOR T2K legacy analysis hists here
  fin = _input_files["VALOR_legacy"]
  h = fin.Get(hname)
  if not h:
    print "WARNING: Could not get hist named [", hname, "] from the input file,", fin.GetName()
  else:
    ROOT.SetOwnership(h, False)       # This stops the object being deleted from memory when it goes out of scope
  return h

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
def sort_map(old_list_in, new_order):
  old_list = copy.deepcopy(old_list_in)
  new_list = []
  old_order = old_list[0]
  if len(list(set(new_order))) != len(new_order):
      print "WARNING: Tried to sort a naming map, but found duplicate names in the ordering, so will not attempt to sort"
      del old_list[0]
      return old_list
  if(len(new_order) > len(old_order)):
    print "ERROR: The length of the new order list was longer than the length old order list. \n New order =", \
      new_order, "\n old order =", old_order
    sys.exit()
  if(len(new_order) == 0):
    print "ERROR: The length of the new order list was 0. New order =", new_order
    sys.exit()
  order_map = []
  for o in new_order:
      if o not in old_order:
          print "ERROR: New order \'{}\'".format(o), "is not in the old order list:", old_order
          sys.exit()
      order_map.append(old_order.index(o))
  for i in range(1, len(old_list)):
      # Skip index 0, as that corresponds to the order list
      new_row = []
      for j in range(0, len(new_order)):
          new_row.append(old_list[i][order_map[j]])
      new_list.append(new_row)
  return new_list

#_______________________________________________________________________
def format_drawn_hist(h, x_max, is_binno):
  if is_binno:
    h.GetXaxis().SetTitle("Bin Number")
  else:
    h.GetXaxis().SetTitle("E_{rec.} (GeV)")
  h.GetYaxis().SetTitle("Ratio to Nominal")
  h.GetYaxis().SetLabelSize(_label_size)
  h.GetXaxis().SetLabelSize(_label_size)
  h.GetYaxis().SetTitleSize(_title_size)
  h.GetXaxis().SetTitleSize(_title_size)
  h.GetYaxis().SetTitleOffset(_title_offset)
  h.GetXaxis().SetTitleOffset(_title_offset)
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
  title_box = ROOT.TPaveText(0.3, 1.0-_margin_t+0.005, 0.7, 1.0, "NDC")
  title_box.AddText(title)
  title_box.SetTextSize(0.05)
  title_box.SetFillColor(0)
  title_box.SetFillStyle(0)
  title_box.SetLineStyle(0)
  title_box.SetLineColor(0)
  ROOT.SetOwnership(title_box, False)       # This stops the object being deleted from memory when it goes out of scope
  return title_box

#_______________________________________________________________________
def get_output_name(output_filename, ip, im, is_binno):
  save_name = output_filename
  if is_binno:
    save_name += "_binno"
  else:
    save_name += "_ereco"
  if _by_mode:
    save_name += "_mode"
  save_name += ".pdf"
  start_pdf = (ip == 0) and (im == 0)
  end_pdf   = (ip == len(_param_map) - 1) and (im == _n_modes - 1)
  if start_pdf:
    save_name += "("
  elif end_pdf:
    save_name += ")"
  return save_name

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

#_______________________________________________________________________
def compare_ratio_hists(param_name_dict, mode_name_dict, sigma_var, comparison_dict):
  already_compared = []
  for iformat, fi in enumerate(_formats):
    for jformat, fj in enumerate(_formats):
      if iformat == jformat:
        continue
      if [jformat, iformat] in already_compared:
        # If we've already compared A to B there's no need to compare B to A
        continue
      for ibin in range(1, comparison_dict[fi].GetNbinsX()):
        diff = 100 * ROOT.TMath.Abs(1.0 - (comparison_dict[fj].GetBinContent(ibin) / comparison_dict[fi].GetBinContent(ibin)))
        if diff > _ratio_diff_warn_threshold_percent:
          if _by_mode:
            print "WARNING: For param =", param_name_dict[fi], "; sample =", _selecs[fi], "; mode =", \
              mode_name_dict[fi], "; Ereco =", comparison_dict[fi].GetBinCenter(ibin), "; variation =", \
              sigma_var, "; difference of", diff, "% in the ratios of input and comparison hist"
          else:
            print "WARNING: For param =", param_name_dict[fi], "; sample =", _selecs[fi], "; Ereco =", \
              comparison_dict[fi].GetBinCenter(ibin), "; variation =", sigma_var, "; difference of", diff, \
              "% in the ratios of input and comparison hist"
      already_compared.append([iformat, jformat])

#_______________________________________________________________________
if __name__=="__main__":
  
  import argparse
  parser = argparse.ArgumentParser(description="Plot single-sample +/- 1 and 3 sigma variations and compare between groups")
  
  parser.add_argument('-o', '--out-file'         , dest='output_filename', default=""     , help="Name of the output image file")
  parser.add_argument('-i', '--in-files'         , dest='in_files'       , nargs="+"      , help="Input file names containing spectra to plot. One file per group")
  parser.add_argument('-f', '--formats'          , dest='formats'        , nargs="+"      , help="Name used to identify the format used for hist naming. Choose from: \'VALOR\', \'MaCh3\', \'p-theta\', \'VALOR_legacy\'")
  parser.add_argument('-s', '--selection'        , dest='selec'          , default=""     , help="The selection name. Choose from: \'FHC1Rmu\', \'FHC1Re\', \'RHC1Rmu\', \RHC1Re\', \'FHCMultiRe\'")
  parser.add_argument('-p', '--params'           , dest='params'         , default=""     , help="File containing a map between param names used by each group. Also used to decide which plots to produce.")
  parser.add_argument('-m', '--modes'            , dest='modes'          , default=""     , help="File containing a map between mode names used by each group. Also used to decide which plots to produce.")
  parser.add_argument('-x', '--x-max'            , dest='x_max'          , type=float     , default=3.0, help="Maximum Ereco value")
  parser.add_argument('-w', '--suppress-warnings', dest='suppress_warn'  , default=False  , action="store_true", help="Usually warnings are printed if any of the ratio hist bin contents are sufficiently different from that of other groups. Use this option to disable these warnings.")
  parser.add_argument(      '--is-unosc'         , dest='is_unosc'       , default=False  , action="store_true", help="The p-theta hists names contain the string \'osc\' or \'unosc\', this flag switches between them. ")
  parser.add_argument(      '--sum-t2k-cc1pi'    , dest='sum_t2k_cc1pi'  , default=False  , action="store_true", help="Sum VALOR T2K CC1piC and CC1pi0 modes when comparing to VALOR T2K legacy")
  parser.add_argument(      '--sum-t2k-ccoth'    , dest='sum_t2k_ccoth'  , default=False  , action="store_true", help="Sum VALOR T2K CCmisc, CCMpi and CCDIS modes when comparing to VALOR T2K legacy")

  o = parser.parse_args()

  main(o.output_filename, o.in_files, o.formats, o.selec, o.params, o.modes, o.x_max, o.suppress_warn, o.is_unosc, o.sum_t2k_cc1pi, o.sum_t2k_ccoth)
