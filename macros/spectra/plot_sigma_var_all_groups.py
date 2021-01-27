#!/usr/bin/env python

# Plot +/- 1 and 3 sigma variations of selected systematic parameters.
# The plots show (varied spectra) / (nominal spectra).
#
# Input one ROOT file per group, each containing TH1D's of either the varied and
# nominal spectra or the ratio of the two (depending on the options set below).
#
# Example usage (see the bottom of this file for explanations of each argument):
#   Compare VALOR, MaCh3 and p-theta (non-mode; FHC single-ring mu-like sample):
#     python plot_sigma_var_all_groups.py \
#     --out-file sigma_var_VALOR_vs_MaCh3_vs_p-theta \
#     --in-files \
#       valor_file.root \
#       mach3_file.root \
#       p-theta_file.root \
#     --formats VALOR MaCh3 p-theta \
#     --selection FHC1Rmu \
#     --x-max 3 \
#     --params maps/valor_mach3_ptheta_xsec_param_name_map.arg
#
#  Compare p-theta and MaCh3 (by mode; FHC nueCC1pi-like sample):
#    python plot_sigma_var_all_groups.py \
#    --out-file sigma_var_p-theta_vs_MaCh3 \
#    --in-files \
#      p-theta_file.root \
#      mach3_file.root \
#    --formats p-theta MaCh3 \
#    --selection FHCMultiRe \
#    --x-max 1.25 \
#    --params maps/valor_mach3_ptheta_xsec_param_name_map.arg \
#    --modes $VALOR_EXPT/config/maps/valor_mach3_ptheta_mode_name_map.arg
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# February 2020

import ROOT
import sys
import re
import copy
from copy import deepcopy
from array import *

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptTitle(0)

# Names used internally within this macro
_valid_formats    = ["VALOR", "MaCh3", "p-theta", "VALOR_legacy"]
_valid_selections = ["FHC1Rmu", "FHC1Re", "RHC1Rmu", "RHC1Re", "FHCMultiRe", "FHCMultiRmu1de", "FHCMultiRmu2de"]
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
_hist_nom_name_patterns["p-theta"]      = "[SAMPLE]_nom[MODE]_[OSC]"
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
_selecs_all["VALOR"       ] = { "FHC1Rmu":"superk__JPARC_FHC__1rmu", "FHC1Re":"superk__JPARC_FHC__1re", "RHC1Rmu":"superk__JPARC_RHC__1rmu", "RHC1Re":"superk__JPARC_RHC__1re", "FHCMultiRe":"superk__JPARC_FHC__multire", "FHCMultiRmu1de":"superk__JPARC_FHC__multirmu1de", "FHCMultiRmu2de":"superk__JPARC_FHC__multirmu2de"}
_selecs_all["MaCh3"       ] = { "FHC1Rmu":"sknumu"                 , "FHC1Re":"sknue"                 , "RHC1Rmu":"anu_sknumubar"          , "RHC1Re":"anu_sknuebar"          , "FHCMultiRe":"sknue1pi"                  , "FHCMultiRmu1de":"sknumucc1pi1de"                , "FHCMultiRmu2de":"sknumucc1pi2de"                }
_selecs_all["p-theta"     ] = { "FHC1Rmu":"numu1R"                 , "FHC1Re":"nue1R"                 , "RHC1Rmu":"numubar1R"              , "RHC1Re":"nuebar1R"              , "FHCMultiRe":"nue1RD"                     }
_selecs_all["VALOR_legacy"] = { "FHC1Rmu":"1Rmu"                   , "FHC1Re":"1Re"                   , "RHC1Rmu":"RHC1Rmu"                , "RHC1Re":"RHC1Re"                , "FHCMultiRe":"MultiRe"                    }

# Whether or not to compare the tweaked parameter values between groups
_compare_param_values = dict()
_compare_param_values["VALOR"         ] = True
_compare_param_values["MaCh3"         ] = False
_compare_param_values["p-theta"       ] = True
_compare_param_values["p-VALOR_legacy"] = False

# The naming format used for the comparison mentioned above
_param_value_array_patterns = dict()
_param_value_array_patterns["VALOR"       ] = "nom_sig_values_vector_[PARAM]"
_param_value_array_patterns["MaCh3"       ] = None
_param_value_array_patterns["p-theta"     ] = "sig_var_pts_[PARAM]"
_param_value_array_patterns["VALOR_legacy"] = None

# copy all the p-theta settings for p-theta internally calculated using either PTheta or ETheta instead of Erec (the actual output files are the same structure)
# the Erec is obviously the default, but adding here for clarity in plots
_ptheta_internal_format_suffixes = ["_Erec", "_ETheta", "_PTheta", "_old", "_new"]
for suffix in _ptheta_internal_format_suffixes:
  _valid_formats += ["p-theta"+suffix]
  _hist_nom_name_patterns    ["p-theta"+suffix] = _hist_nom_name_patterns    ["p-theta"]
  _hist_tweaked_name_patterns["p-theta"+suffix] = _hist_tweaked_name_patterns["p-theta"]
  _ratio_already_taken       ["p-theta"+suffix] = _ratio_already_taken       ["p-theta"]
  _var_names                 ["p-theta"+suffix] = _var_names                 ["p-theta"]
  _selecs_all                ["p-theta"+suffix] = _selecs_all                ["p-theta"]


#==================================================================================
# STUFF THAT CHANGES DEPENDING ON THE HIST NAMES AND FORMATS FOR EACH GROUP (END)
#==================================================================================

_formats          = None
_input_files      = dict()
_selecs           = dict()
_by_mode          = False
_mode_map         = None
_n_modes          = None
_param_map        = None
_rebin_to_format  = ""
_binning_ref_axis = None

_line_col = dict()
_line_col["p3"] = ROOT.kBlue
_line_col["p1"] = ROOT.kCyan
_line_col["m1"] = ROOT.kViolet
_line_col["m3"] = ROOT.kRed

_default_nom_ratio_col = ROOT.kGray

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
def main(output_filename, in_files, formats, selec, params_file, modes_file, x_max, suppress_warn, is_etrue, compare_sigma, plot_nom_ratio, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth, rebin_to_format):
  setup(in_files, formats, selec, params_file, modes_file, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth, rebin_to_format)

  # Should we check that each group has evaluated each parameter at the same values?
  compare_counter = 0
  for format in _formats:
    if _compare_param_values[format] and _param_value_array_patterns[format]:
      compare_counter += 1
  do_compare_param_values = (compare_counter > 1 and compare_sigma)

  if do_compare_param_values:
    compare_sigma_between_groups()

  # Should we superimpose ratio plots of the nominal spectra onto the sigma variation plots? (e.g. plot valor nominal / mach3 nominal etc.)
  nom_ratio_format_pairs = []
  if plot_nom_ratio:
    ROOT.gStyle.SetPadTickY(0)  # Don't draw y-axis tick marks on the right hand side, as we'll draw a separate y-axis there for the plots of the ratio of nominal spectra
    for iformat, fi in enumerate(_formats):
      for jformat, fj in enumerate(_formats):
        if iformat == jformat:
          continue
        if [fj, fi] in nom_ratio_format_pairs:
          # If we've already compared A to B there's no need to compare B to A
          continue
        nom_ratio_format_pairs.append([fi, fj])

  for ip in range(-1, len(_param_map)):
    plot_nominal = (ip == -1)
    param_name_dict = get_param_name_dict(ip if not plot_nominal else 0)
    for im in range(0, _n_modes):
      mode_name_dict = None
      if _by_mode:
        mode_name_dict = get_mode_name_dict(im)
        if(should_skip_mode(mode_name_dict)):
          continue

      canv_name = "sigma_var_ratio_" + _selecs[_formats[0]] + "_" + (param_name_dict[_formats[0]] if not plot_nominal else "nominal")
      if _by_mode:
        canv_name += "_" + mode_name_dict[_formats[0]]
      if is_etrue:
        canv_ereco = ROOT.TCanvas(canv_name + "_Etrue", canv_name + "_Etrue", 800, 600)
      else:
        canv_ereco = ROOT.TCanvas(canv_name + "_Ereco", canv_name + "_Ereco", 800, 600)
      canv_binno = ROOT.TCanvas(canv_name + "_binno", canv_name + "_binno", 800, 600)
      canv_ereco.SetMargin(_margin_l, _margin_r, _margin_t, _margin_b)
      canv_binno.SetMargin(_margin_l, _margin_r, _margin_t, _margin_b)
      print "Making plot named", canv_name

      leg = ROOT.TLegend(1.0-_margin_r, 0.65 - (0.25 * int(plot_nom_ratio and not plot_nominal)), 0.995, 1.0-_margin_t)
      if not plot_nominal:
        leg.AddEntry(None, "Variations", "")
      h_list       = []
      h_list_binno = []
      hist_max_var = -999

      # Get the +/- 1 and 3 sigma variations and take the ratio to the nominal
      first_drawn_var = True
      for ivar, sigma_var in enumerate(_variations):
        comparison_dict = dict()
        first_drawn_format = True
        for iformat, f in enumerate(_formats):
          if not plot_nominal:
            h = get_ratio_hist(sigma_var, param_name_dict, mode_name_dict, iformat)
          else:
            h = get_nominal_for_plot(mode_name_dict, iformat)

          if not h:
            # Skip any hists that we couldn't read from the input file
            continue
          h.SetName(canv_name + ("_" + sigma_var if not plot_nominal else "") + "_" + f)
          h_binno = hist_to_bin_number(h)   # Also plot with the bin number on the x-axis

          draw_opt = "HIST"
          if not first_drawn_var or not first_drawn_format:
            draw_opt = "HIST,SAME"
          else:
            format_drawn_hist(h, x_max, False, is_etrue, plot_nominal)
            format_drawn_hist(h_binno, h.FindBin(x_max) -1, True, is_etrue, plot_nominal)
          canv_ereco.cd()
          h.Draw(draw_opt)
          canv_binno.cd()
          h_binno.Draw(draw_opt)

          if first_drawn_format and not plot_nominal:
            # Label for each variation
            leg.AddEntry(h, get_label(sigma_var), "l")
          if ivar == len(_variations) -1:
            # Label for each group
            if first_drawn_format:
              if not plot_nominal:
                leg.AddEntry(None, "", "")  # Blank line
              leg.AddEntry(None, "Groups", "")
            h_copy = ROOT.TH1D(h)
            h_copy.SetLineColor(1)  # Keep the group labels in black, so the focus is just on the line style
            leg.AddEntry(h_copy, f, "l")
            h_list.append(h_copy)

          h_list.append(h)
          h_list_binno.append(h_binno)
          comparison_dict[f] = h
          hist_max_var = max(hist_max_var, get_max_var_in_range(h, x_max, not plot_nominal))
          first_drawn_format = False
          first_drawn_var = False

        # Throw a warning if any of the ratio hists have a difference greater than some threshold
        if not suppress_warn and not plot_nominal:
          compare_ratio_hists(param_name_dict, mode_name_dict, sigma_var, comparison_dict)
        del comparison_dict

      if not plot_nominal:
        # Set the y-axis range to be symmetric around ratio = 1
        h_list      [0].SetMaximum(1 + 1.1 * hist_max_var)
        h_list      [0].SetMinimum(1 - 1.1 * hist_max_var)
        h_list_binno[0].SetMaximum(1 + 1.1 * hist_max_var)
        h_list_binno[0].SetMinimum(1 - 1.1 * hist_max_var)
      else:
        h_list      [0].SetMaximum(1.1 * hist_max_var)
        h_list      [0].SetMinimum(0.)
        h_list_binno[0].SetMaximum(1.1 * hist_max_var)
        h_list_binno[0].SetMinimum(0.)

      # Superimpose plots of (group 1 nominal) / (group 2 nominal) etc.
      h_list_nom_ratio = []
      if plot_nom_ratio:
        # Iterate over all pairs of groups twice: First to get the max and min values so we can scale the hists
        # to fit in the same plot as the sigma variation ratios, then second to actually draw the hists.
        max_nom_ratio = -999
        for i, format_pair in enumerate(nom_ratio_format_pairs):
          h_ratio_of_nominals = get_ratio_of_nominal_hists(format_pair, mode_name_dict, i, plot_nominal)
          h_ratio_of_nominals.SetName(canv_name + "_nom_ratio_" + format_pair[0] + "_" + format_pair[1])
          max_nom_ratio = max(max_nom_ratio, get_max_var_in_range(h_ratio_of_nominals, x_max, not plot_nominal))
          h_list_nom_ratio.append(h_ratio_of_nominals)
        axis_value_ymax = 1 + 1.1 * max_nom_ratio
        axis_value_ymin = 1 - 1.1 * max_nom_ratio
        if plot_nominal:
          axis_value_ymax = 1.1 * max_nom_ratio
          axis_value_ymin = 0
        if axis_value_ymax == 1 and axis_value_ymin == 1:
          axis_value_ymin = 0.9
          axis_value_ymax = 1.1

        for i, format_pair in enumerate(nom_ratio_format_pairs):
          if i == 0:
            leg.AddEntry(None, "", "")
            leg.AddEntry(None, "#splitline{Ratio of}{nominal spectra}", "")

          h_ratio_of_nominals = h_list_nom_ratio[i]
          scale_factor = 1.0
          if max_nom_ratio != 0:
            if hist_max_var == 0:
              scale_factor = (ROOT.gPad.GetUymax() - ROOT.gPad.GetUymin()) / max_nom_ratio
            else:
              scale_factor = hist_max_var / max_nom_ratio

            if plot_nominal:
              h_ratio_of_nominals.Scale(scale_factor)
            else:
              rescale_around_y_eq_one(h_ratio_of_nominals, scale_factor) # Scale so it will fit on the sigma variations plot

          canv_ereco.cd()
          h_ratio_of_nominals.Draw("HIST,SAME")

          canv_binno.cd()
          h_ratio_of_nominals_binno = hist_to_bin_number(h_ratio_of_nominals)
          h_ratio_of_nominals_binno.Draw("HIST,SAME")

          leg.AddEntry(h_ratio_of_nominals, "#frac{" + format_pair[0] + "}{" + format_pair[1] + "}", "l")
          h_list_nom_ratio.append(h_ratio_of_nominals)

        axis_ymin = h_list[0].GetMinimum()
        axis_ymax = h_list[0].GetMaximum()
        ndiv      = h_list[0].GetYaxis().GetNdivisions()
        if axis_ymin == axis_ymax:
          axis_ymin = 0
          axis_ymax = 2
          ndiv = 510
        axis_ratio_of_nominals = ROOT.TGaxis(x_max, axis_ymin, x_max, axis_ymax,    # The position of the axis, in the coordinate system used by the sigma var ratio hists
                                             axis_value_ymin, axis_value_ymax, ndiv, "+C")
        axis_ratio_of_nominals.SetLineColor  (_default_nom_ratio_col)
        axis_ratio_of_nominals.SetLabelColor (_default_nom_ratio_col)
        axis_ratio_of_nominals.SetLabelSize  (h_list[0].GetYaxis().GetLabelSize())
        axis_ratio_of_nominals.SetLabelFont  (h_list[0].GetYaxis().GetLabelFont())
        axis_ratio_of_nominals.SetLabelOffset(-0.05)

        xmax_binno = h_list_nom_ratio[0].FindBin(x_max) -1
        axis_ratio_of_nominals_binno = ROOT.TGaxis(xmax_binno, axis_ymin, xmax_binno, axis_ymax,    # The position of the axis, in the coordinate system used by the sigma var ratio hists
                                                   axis_value_ymin, axis_value_ymax, ndiv, "+C")
        axis_ratio_of_nominals_binno.SetLineColor  (_default_nom_ratio_col)
        axis_ratio_of_nominals_binno.SetLabelColor (_default_nom_ratio_col)
        axis_ratio_of_nominals_binno.SetLabelSize  (h_list[0].GetYaxis().GetLabelSize())
        axis_ratio_of_nominals_binno.SetLabelFont  (h_list[0].GetYaxis().GetLabelFont())
        axis_ratio_of_nominals_binno.SetLabelOffset(-0.05)

        canv_ereco.cd()
        axis_ratio_of_nominals.Draw()

        canv_binno.cd()
        axis_ratio_of_nominals_binno.Draw()

        title_box_ratio_of_nominals = ROOT.TPaveText(0.92, 0.10, 1.0, 0.37, "NDC")
        text = title_box_ratio_of_nominals.AddText("#splitline{Ratio of}{Nominal Spectra}")
        text.SetTextAngle(270)
        text.SetTextAlign(22)
        title_box_ratio_of_nominals.SetTextSize(h_list[0].GetYaxis().GetTitleSize())
        title_box_ratio_of_nominals.SetTextColor(_default_nom_ratio_col)
        title_box_ratio_of_nominals.SetFillColor(0)
        title_box_ratio_of_nominals.SetFillColorAlpha(0, 0)
        title_box_ratio_of_nominals.SetFillStyle(0)
        title_box_ratio_of_nominals.SetLineStyle(0)
        title_box_ratio_of_nominals.SetLineColor(0)
        ROOT.SetOwnership(title_box_ratio_of_nominals, False)

        canv_ereco.cd()
        title_box_ratio_of_nominals.Draw("SAME")

        canv_binno.cd()
        title_box_ratio_of_nominals.Draw("SAME")

      title_box = make_title_box(get_title_box_text(param_name_dict if not plot_nominal else None, mode_name_dict))
      canv_ereco.cd()
      title_box.Draw("SAME")
      leg.Draw("SAME")
      canv_binno.cd()
      title_box.Draw("SAME")
      leg.Draw("SAME")

      canv_ereco.cd()
      canv_ereco.Update()
      print_opt = ""
      if plot_nominal:
        if _by_mode:
          print_opt = "Title:Nominal_" + mode_name_dict[_formats[0]]
        else:
          print_opt = "Title:Nominal"
      else:
        if _by_mode:
          print_opt = "Title:" + param_name_dict[_formats[0]] + "_" + mode_name_dict[_formats[0]]
        else:
          print_opt = "Title:" + param_name_dict[_formats[0]]
      canv_ereco.Print(get_output_name(output_filename, selec, ip, im, False, is_etrue), print_opt)
      canv_binno.cd()
      canv_binno.Update()
      canv_binno.Print(get_output_name(output_filename, selec, ip, im, True, is_etrue), print_opt)

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
def setup(in_files, formats, selec, params_file, modes_file, is_unosc, sum_t2k_cc1pi, sum_t2k_ccoth, rebin_to_format):
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
  global _rebin_to_format
  _by_mode  = modes_file != ""
  _formats  = formats
  _is_unosc = is_unosc
  _rebin_to_format = rebin_to_format

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
  _param_map = add_entries_for_internal_formats(_param_map, "p-theta", _ptheta_internal_format_suffixes)
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
    _mode_map = add_entries_for_internal_formats(_mode_map, "p-theta", _ptheta_internal_format_suffixes)
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
  print_param_name = param_name_dict[_formats[0]] if param_name_dict else ""
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

  set_empty_bins_to_value(hratio, 1.0)

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
def get_ratio_of_nominal_hists(format_pair, mode_name_dict, style_index, plot_nominal):
  hratio = None

  if _ratio_already_taken[format_pair[0]] or _ratio_already_taken[format_pair[1]]:
    print "ERROR: Cannot get ratio of nominal hists between groups", format_pair[0], "and", format_pair[1]
    sys.exit()

  # Make the ratio hist
  hnumerator   = get_nominal(mode_name_dict, format_pair[0])
  hdenominator = get_nominal(mode_name_dict, format_pair[1])
  if not hnumerator or not hdenominator:
    return None
  hratio = ROOT.TH1D(hnumerator)
  hratio.Divide(hdenominator)

  if plot_nominal:
    set_empty_bins_to_value(hratio, 0.0)
  else:
    set_empty_bins_to_value(hratio, 1.0)

  # Set the line style and colour differently for each pair of groups
  # col_offset = style_index
  # # The kGray colour offset values in TColor range from 0 to +3 (inclusive), so keep within range
  # if col_offset > 3 or col_offset < 0:
  #   print "WARNING: Setting the kGray TColor offset to", col_offset, ", which may not be a valid colour."
  hratio.SetLineColor(_default_nom_ratio_col)
  hratio.SetLineStyle(len(_formats) + 1 + style_index)  # To Make sure that we don't get the same style as any of the sigma var hists
  hratio.SetLineWidth(1 + style_index)

  ROOT.SetOwnership(hratio, False)       # This stops the object being deleted from memory when it goes out of scope
  return hratio

#_______________________________________________________________________
def get_nominal_for_plot(mode_name_dict, iformat):
  format = _formats[iformat]
  hnom   = get_nominal(mode_name_dict, format)
  if not hnom:
    return None

  hnom.SetLineStyle(1 + iformat)
  hnom.SetLineWidth(1 + iformat)

  ROOT.SetOwnership(hnom, False)       # This stops the object being deleted from memory when it goes out of scope
  return hnom

#_______________________________________________________________________
def set_empty_bins_to_value(h, value):
  for ibin in range(1, h.GetNbinsX()+1):
    if h.GetBinContent(ibin) < _a_small_num:
      h.SetBinContent(ibin, value)

#_______________________________________________________________________
def get_nominal(mode_name_dict, format):
  return get_hist(False, None, None, mode_name_dict, format)

#_______________________________________________________________________
def get_tweaked(sigma_var, param_name_dict, mode_name_dict, format):
  return get_hist(True, sigma_var, param_name_dict, mode_name_dict, format)

#_______________________________________________________________________
def get_hist(get_tweaked_hist, sigma_var, param_name_dict, mode_name_dict, format):
  global _binning_ref_axis
  global _rebin_to_format
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
    h = get_ptheta_hist(hname, format)
  elif format == "VALOR_legacy":
    h = get_t2k_legacy_hist(hname)
  else:
    for suffix in _ptheta_internal_format_suffixes:
      if format == "p-theta"+suffix:
        h = get_ptheta_hist(hname, format)
    if h is None:
      print "ERROR: Unrecognised format \'{0}\'".format(format)
      sys.exit()

  if _rebin_to_format:
    if format == _rebin_to_format:
      if not _binning_ref_axis:
        # Store the desired binning
        buffer = h.GetXaxis().GetXbins().GetArray()
        _binning_ref_axis = [buffer[i] for i in range(0, h.GetNbinsX()+1)]
        assert _binning_ref_axis
    else:
      h = rebin(h)

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
def get_ptheta_hist(hname, format):
  # Put in any special cases related to getting p-theta hists here
  fin = _input_files[format]
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
def add_entries_for_internal_formats(tuples, matchFormat, suffixes):
  formats = tuples[0]
  newtuples = []
  for i,f in enumerate(formats):
    if f == matchFormat:
      iorg = i
      for row,tup in enumerate(tuples):
        newtup = tup
        for suffix in suffixes:
          if row == 0:
            newtup += (matchFormat+suffix,)
          else:
            newtup += (tup[iorg],) # copy the value for matchFormat
        newtuples.append(newtup)
  return newtuples

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
def format_drawn_hist(h, x_max, is_binno, is_etrue, plot_nominal):
  if is_binno:
    h.GetXaxis().SetTitle("Bin Number")
  elif is_etrue:
    h.GetXaxis().SetTitle("E_{true} (GeV)")
  else:
    h.GetXaxis().SetTitle("E_{rec.} (GeV)")
  if plot_nominal:
    h.GetYaxis().SetTitle("Nominal Event Rate")
  else:
    h.GetYaxis().SetTitle("Ratio to Nominal")
  h.GetYaxis().SetLabelSize(_label_size)
  h.GetXaxis().SetLabelSize(_label_size)
  h.GetYaxis().SetTitleSize(_title_size)
  h.GetXaxis().SetTitleSize(_title_size)
  h.GetYaxis().SetTitleOffset(_title_offset)
  h.GetXaxis().SetTitleOffset(_title_offset)
  h.GetXaxis().SetRangeUser(h.GetXaxis().GetBinLowEdge(1), x_max)

#_______________________________________________________________________
def get_max_var_in_range(h, x_max, is_ratio_hist = True):
  max_var = -999
  for i in range(1, h.GetNbinsX()):
    if h.GetXaxis().GetBinUpEdge(i) <= x_max:
      if is_ratio_hist:
        max_var = max(max_var, ROOT.TMath.Abs(h.GetBinContent(i)-1))  # -1 to make the plot symmetric around ratio = 1
      else:
        max_var = max(max_var, h.GetBinContent(i))
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
def get_output_name(output_filename, selec, ip, im, is_binno, is_etrue):
  save_name = output_filename
  save_name += "_" + selec
  if is_binno:
    save_name += "_binno"
  elif is_etrue:
    save_name += "_etrue"
  else:
    save_name += "_ereco"
  if _by_mode:
    save_name += "_mode"
  save_name += ".pdf"
  start_pdf = (ip == -1) and (im == 0)
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
def compare_sigma_between_groups():
  maxlen_format = max([len(f) for f in _formats])
  for ip in range(0, len(_param_map)):
    param_name_dict = get_param_name_dict(ip)
    print param_name_dict[_formats[0]]
    for i, f in enumerate(_formats):
      if _compare_param_values[f] and _param_value_array_patterns[f]:
        name = _param_value_array_patterns[f].replace("[PARAM]", param_name_dict[f])
        param_values_vector = _input_files[f].Get(name)
        if not param_values_vector:
          print "ERROR: Could not find an object named", name, "in file", _input_files[f].GetName()
          sys.exit()
        nom = param_values_vector[2]    # We're not comparing the nominal values between groups, because some normalise this value to the generated value while other groups don't.
        sigma = param_values_vector[3] - param_values_vector[2]
        #row = "\t" + f + ": " + " "*(maxlen_format-len(f)) + " nom = " + str(nom) + ", sigma = " + str(sigma)
        row = "\t" + f + ": " + " "*(maxlen_format-len(f)) + ", sigma = " + str(sigma)
        print row

#_______________________________________________________________________
def rescale_around_y_eq_one(h, scale_factor):
  for i in range(1, h.GetNbinsX()):
    scaled_content = 1 + (h.GetBinContent(i) - 1.0) * scale_factor
    h.SetBinContent(i, scaled_content)

#_______________________________________________________________________
def rebin(h):
  global _binning_ref_axis
  global _rebin_to_format

  # Make sure we have the desired binning
  if not _binning_ref_axis:
    get_hist(False, None, None, None, _rebin_to_format)

  # The rebin is only possible if the desired bin edges lie on the current bin edges
  buffer = h.GetXaxis().GetXbins().GetArray()
  curr_binning = [buffer[i] for i in range(0, h.GetNbinsX()+1)]
  for edge in _binning_ref_axis:
    if edge not in curr_binning:
      print "WARNING: Tried to rebin hist \'{}\'".format(h.GetName()), "but the desired bin edges", \
            _binning_ref_axis, "don't align with the current binning", curr_binning, ". Will not rebin."
      return

  new_binning_array = array('d', _binning_ref_axis)
  hnew = ROOT.TH1D(h)
  hnew.Reset()
  hnew.SetBins(len(new_binning_array)-1, new_binning_array)

  for ibin_old in range(1, h.GetNbinsX()+1):
    centre  = h.GetBinCenter (ibin_old)
    content = h.GetBinContent(ibin_old)
    hnew.Fill(centre, content)

  h = copy.deepcopy(hnew)
  del hnew
  return h

#_______________________________________________________________________
if __name__=="__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Plot single-sample +/- 1 and 3 sigma variations and compare between groups")

  parser.add_argument('-o', '--out-file'         , dest='output_filename', default=""     , help="Name of the output image file")
  parser.add_argument('-i', '--in-files'         , dest='in_files'       , nargs="+"      , help="Input file names containing spectra to plot. One file per group")
  parser.add_argument('-f', '--formats'          , dest='formats'        , nargs="+"      , help="Name used to identify the format used for hist naming. Choose from: \'VALOR\', \'MaCh3\', \'p-theta\', \'VALOR_legacy\'")
  parser.add_argument('-s', '--selection'        , dest='selec'          , default=""     , help="The selection name. Choose from: \'FHC1Rmu\', \'FHC1Re\', \'RHC1Rmu\', \RHC1Re\', \'FHCMultiRe\', \'FHCMultiRmu1de\', \'FHCMultiRmu2de\'")
  parser.add_argument('-p', '--params'           , dest='params'         , default=""     , help="File containing a map between param names used by each group. Also used to decide which plots to produce.")
  parser.add_argument('-m', '--modes'            , dest='modes'          , default=""     , help="File containing a map between mode names used by each group. Also used to decide which plots to produce.")
  parser.add_argument(      '--rebin-to-format'  , dest='rebin_to_format', default=""     , help="If one group has coarser binning than the others, supply the group name here to rebin the hists from the other groups, if possible.")
  parser.add_argument('-x', '--x-max'            , dest='x_max'          , type=float     , default=3.0, help="Maximum Ereco value")
  parser.add_argument('-w', '--suppress-warnings', dest='suppress_warn'  , default=False  , action="store_true", help="Usually warnings are printed if any of the ratio hist bin contents are sufficiently different from that of other groups. Use this option to disable these warnings.")
  parser.add_argument('-c', '--compare-sigma'    , dest='compare_sigma'  , default=False  , action="store_true", help="Print a comparison between the sigma values used by each group for each param.")
  parser.add_argument('-n', '--plot-nom-ratio'   , dest='plot_nom_ratio' , default=False  , action="store_true", help="Plot the ratios of the nominal spectra for each group superimposed onto the systematic variations")
  parser.add_argument(      '--is-etrue'         , dest='is_etrue'       , default=False  , action="store_true", help="Set to true if the histograms use etrue binning, and thus the x label should be E_{true}")
  parser.add_argument(      '--is-unosc'         , dest='is_unosc'       , default=False  , action="store_true", help="The p-theta hists names contain the string \'osc\' or \'unosc\', this flag switches between them. ")
  parser.add_argument(      '--sum-t2k-cc1pi'    , dest='sum_t2k_cc1pi'  , default=False  , action="store_true", help="Sum VALOR T2K CC1piC and CC1pi0 modes when comparing to VALOR T2K legacy")
  parser.add_argument(      '--sum-t2k-ccoth'    , dest='sum_t2k_ccoth'  , default=False  , action="store_true", help="Sum VALOR T2K CCmisc, CCMpi and CCDIS modes when comparing to VALOR T2K legacy")

  o = parser.parse_args()

  main(o.output_filename, o.in_files, o.formats, o.selec, o.params, o.modes, o.x_max, o.suppress_warn, o.is_etrue, o.compare_sigma, o.plot_nom_ratio, o.is_unosc, o.sum_t2k_cc1pi, o.sum_t2k_ccoth, o.rebin_to_format)
