#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys
import copy
from copy import deepcopy

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

_hatch_fill_style = 3003    # For filling in FC intervals

#________________________________________________________________________________________________________
def plot_contours_1d(dchisq_objects, labels, pdf_filename, is_t2k_only, colours, styles, prelim, title_x, title_y,
                     line_width, text_size, conf_levels, leg_coords, titlesize,
                     fc_files, fc_graph_names, fc_line_styles, fc_line_colours, fc_draw_critvals, fc_draw_boxes,
                     is_inverted=False, force_max_y = -1):
    max_val = 0.0
    min_val = 999
    for dchisq_obj in dchisq_objects:
        max_val = max(max_val, dchisq_obj.GetMaximum())
        min_val = min(min_val, dchisq_obj.GetMinimum())
    max_val = np.ceil(max_val)
    max_val = min(max_val, 25.0)    # Not meaningful above 5 sigma, so truncate the axis
    if force_max_y > 0:
        max_val = force_max_y
    min_val = np.floor(min_val)
    for dchisq_obj in dchisq_objects:
        dchisq_obj.GetYaxis().SetRangeUser(0, max_val)

    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)

    if "mixing13" in pdf_filename:
        sinsqth13_nom = 0.0212
        sinsqth13_err = 0.0008
        X1 = sinsqth13_nom - sinsqth13_err
        X2 = sinsqth13_nom + sinsqth13_err
        print "*** test:", X1,min_val,X2,max_val
        breactors = ROOT.TBox(X1,min_val,X2,max_val)
        breactors.SetFillColor(ROOT.kYellow+1)
        dchisq_objects[0].Draw("axis")
        breactors.Draw("same")

    options = ""
    if "mixing13" in pdf_filename:
        options = "same"
    if len(dchisq_objects) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1
    for i, dchisq_obj in enumerate(dchisq_objects):
        if not colours or not styles:
            dchisq_obj.SetLineColor(line_color)
        else:
            # Custom line colours and styles have been specified
            dchisq_obj.SetLineColor(colours[i])
            dchisq_obj.SetLineStyle(styles[i])
        if line_width:
            dchisq_obj.SetLineWidth(line_width)
        if title_x:
            dchisq_obj.GetXaxis().SetTitle(title_x)
            #dchisq_obj.GetXaxis().SetTitleOffset(0.84)
            dchisq_obj.GetXaxis().SetTitleOffset(0)
            dchisq_obj.GetXaxis().SetTitleSize(titlesize)
        if title_y:
            dchisq_obj.GetYaxis().SetTitle(title_y)
            dchisq_obj.GetYaxis().SetTitleSize(titlesize)
        dchisq_obj.GetYaxis().SetTitleOffset(0.8)
        dchisq_obj.Draw("l" + options)
        options = "same"
        line_color += 1

    if not text_size:
        text_size = 0.03
    lymax = 0.86
    lymin = lymax - (1 + len(dchisq_objects)) * text_size
    lxmin = 0.15
    lxmax = 0.3
    if "mixing13" in pdf_filename:
        lxmin = 0.65
        lxmax = 0.8

    if len(leg_coords) != 0:
        # The user has specified legend coordinates, so use them instead of the defaults
        lxmin = leg_coords[0]
        lymin = leg_coords[1]
        lxmax = leg_coords[2]
        lymax = leg_coords[3]

    legend = ROOT.TLegend(lxmin, lymin, lxmax, lymax)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(132)
    legend.SetTextSize(text_size)
    if len(labels) > 0:
        for i, dchisq_obj in enumerate(dchisq_objects):
            legend.AddEntry(dchisq_obj, labels[i], "l")

    if not fc_files:
        # Draw constant Delta chisq critical value line
        xlow = dchisq_objects[0].GetXaxis().GetBinLowEdge(1)
        xhigh = dchisq_objects[0].GetXaxis().GetBinUpEdge(dchisq_objects[0].GetNbinsX())
        xpos_frac = 0.8
        if "mixing13" in pdf_filename:
            xpos_frac = 0.1
        xpos = dchisq_objects[0].GetXaxis().GetBinLowEdge( int(xpos_frac * dchisq_objects[0].GetNbinsX()) )

        critval_objects = []
        critval_draw_opt = ""
        for conf_level in conf_levels:
            critval = util.DeltaChisquare(conf_level, 1)    # The 2nd arg is the number of dimensions
            critical = ROOT.TLine(xlow, critval, xhigh, critval)
            critical.SetLineColor(ROOT.kBlue)
            critical.SetLineStyle(1)
            critical.SetLineWidth(2)
            critval_objects.append(critical)
            critical.Draw(critval_draw_opt)
            critval_draw_opt = "same"

            latLine = ROOT.TLatex(xpos, critval + 0.15, util.ConfLevelAsStringRooTex(conf_level))
            latLine.SetTextColor(ROOT.kBlue)
            latLine.SetTextSize(0.05)
            critval_objects.append(latLine)
            latLine.Draw("same")
    else:
        # Draw boxes to represent the Feldman-Cousins allowed regions
        fc_tgraphs  = []
        hatch_boxes = []
        for i, dchisq_obj in enumerate(dchisq_objects):
            fc_tfile = ROOT.TFile(fc_files[i])
            fc_critval_tgraph = fc_tfile.Get(fc_graph_names[i])
            if not fc_critval_tgraph:
                print "ERROR: Could not find a TGraph named \'{}\'".format(fc_graph_names[i]), \
                      "in the file", fc_files[i]
                sys.exit()
            fc_tgraphs.append(fc_critval_tgraph)

            if fc_draw_boxes:
                [left_boundary, right_boundary] = calc_1d_fc_xrange(dchisq_obj, fc_critval_tgraph)
                hist_left_edge  = dchisq_obj.GetXaxis().GetBinLowEdge(1)
                hist_right_edge = dchisq_obj.GetXaxis().GetBinUpEdge(dchisq_obj.GetNbinsX())
                hatch_box_max_y = 0.66 * max_val
                if left_boundary > right_boundary:
                    hatch_interval  = ROOT.TBox(left_boundary , 0, hist_right_edge, hatch_box_max_y)
                    hatch_interval.SetLineWidth(0)
                    hatch_interval.SetFillStyle(_hatch_fill_style)
                    hatch_interval.SetFillColor(dchisq_obj.GetLineColor())
                    hatch_boxes.append(hatch_interval)
                    hatch_interval.Draw("same")

                    hatch_interval2 = ROOT.TBox(hist_left_edge, 0, right_boundary , hatch_box_max_y)
                    hatch_interval2.SetLineWidth(0)
                    hatch_interval2.SetFillStyle(_hatch_fill_style)
                    hatch_interval2.SetFillColor(dchisq_obj.GetLineColor())
                    hatch_boxes.append(hatch_interval2)
                    hatch_interval2.Draw("same")
                else:
                    hatch_interval  = ROOT.TBox(left_boundary , 0, right_boundary , hatch_box_max_y)
                    hatch_interval.SetLineWidth(0)
                    hatch_interval.SetFillStyle(_hatch_fill_style)
                    hatch_interval.SetFillColor(dchisq_obj.GetLineColor())
                    hatch_boxes.append(hatch_interval)
                    hatch_interval.Draw("same")

        # Draw the FC critical values in a separate loop to the FC boxes to ensure that the critical value
        # lines are always in front of the boxes.
        if fc_draw_critvals:
            for i, dchisq_obj in enumerate(dchisq_objects):
                fc_tgraphs[i].SetLineColor(fc_line_colours[i])
                fc_tgraphs[i].SetLineStyle(fc_line_styles[i])
                fc_tgraphs[i].Draw("L,same")

                if i == 0:
                    fc_tgraph_leg_label  = "-2 #Delta lnL_{crit.} " # Have to split this up so that the RooTex subscript/superscript doesn't mess up format()
                    fc_tgraph_leg_label += "({} CL)".format(util.ConfLevelAsStringRooTex(conf_levels[0]))
                    legend.AddEntry(fc_tgraphs[i], fc_tgraph_leg_label, "l")

    if prelim:
        title = ROOT.TLatex(0.5, 0.91, "T2K Run 1-9 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")

    if "mixing13" in pdf_filename:
        legend.AddEntry(breactors, "PDG 2018", "f")

    legend.Draw()

    canvas.Print(pdf_filename + ".pdf")
    canvas.Print(pdf_filename + ".root")

#________________________________________________________________________________________________________
def plot_contours_2d(dchisq_objs, labels, theta23, dm23, pdf_filename, colours, styles, best_fit_styles,
                     best_fit_labels, prelim, title_x, title_y, line_width, marker_size, text_size,
                     conf_levels, leg_coords, bf_leg_coords, titlesize, draw_reactor=False, is_inverted=False):

    lmargin = 0.13
    rmargin = 0.08
    tmargin = 0.05
    bmargin = 0.13

    canvas = ROOT.TCanvas("canvas2d", "canvas2d", 700, 500)
    canvas.SetLeftMargin  (lmargin)
    canvas.SetRightMargin (rmargin)
    canvas.SetTopMargin   (tmargin)
    canvas.SetBottomMargin(bmargin)

    if not text_size:
        text_size = 0.03

    lymax = 0.97 - tmargin
    lymin = lymax - (len(dchisq_objs)*len(conf_levels) + 1) * text_size
    lxmax = 0.90 - rmargin
    lxmin = 0.57 + lmargin
    lx_tsize_scale = 0.3*(1 - text_size/0.03) + 0.07*(1 - float(max([len(i) for i in labels]))/15.0)  # Adjust the legend left edge of the legend depending on the text size and the size of the labels
    use_lxmin = lxmin + lx_tsize_scale
    if len(leg_coords) != 0:
        # Use user-specified legend coordinates
        lxmin = leg_coords[0]
        lymin = leg_coords[1]
        lxmax = leg_coords[2]
        lymax = leg_coords[3]
        use_lxmin = lxmin

    maxlen_bf_labels = 0
    len_bf_empty = 0
    if best_fit_labels:
        maxlen_bf_labels = float(max([len(i) for i in best_fit_labels]))
        len_bf_empty = len(filter(None, best_fit_labels))
    bf_lxmin = lxmin - 0.15 + lx_tsize_scale + 0.07*(1 - maxlen_bf_labels/15.0) - 0.03
    bf_lxmax = lxmin + lx_tsize_scale - 0.03
    bf_lymin = lymax - (len_bf_empty + int(draw_reactor) + 1) * text_size
    bf_lymax = lymax
    if len(bf_leg_coords) != 0:
        # Use user-specified legend coordinates
        bf_lxmin = bf_leg_coords[0]
        bf_lymin = bf_leg_coords[1]
        bf_lxmax = bf_leg_coords[2]
        bf_lymax = bf_leg_coords[3]

    legend = ROOT.TLegend(use_lxmin, lymin, lxmax, lymax)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(132)
    legend.SetTextSize(text_size)
    if "_{32}" in dchisq_objs[-1].GetYaxis().GetTitle():
        if is_inverted:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{13} (eV^{2}c^{-4})")
        else:
            dchisq_objs[-1].GetYaxis().SetTitle("#Deltam^{2}_{32} (eV^{2}c^{-4})")
    if title_x:
            dchisq_objs[-1].GetXaxis().SetTitle(title_x)
            dchisq_objs[-1].GetXaxis().SetTitleSize(titlesize)
    if title_y:
        dchisq_objs[-1].GetYaxis().SetTitle(title_y)
        dchisq_objs[-1].GetYaxis().SetTitleSize(titlesize)
        
    dchisq_objs[-1].GetYaxis().SetTitleOffset(0) #1.3
    dchisq_objs[-1].GetXaxis().SetTitleOffset(0)  #0.9
    dchisq_objs[-1].SetTitle("")
    dchisq_objs[-1].Draw("axis")

    # Draw the reactor constraint used in T2K, as a vertical yellow band
    sinsqth13_nom = 0.0212
    sinsqth13_err = 0.0008
    breactors_label = "PDG 2018"
    X1 = sinsqth13_nom + sinsqth13_err
    X2 = sinsqth13_nom - sinsqth13_err
    Y1 = dchisq_objs[-1].GetYaxis().GetXmin()
    Y2 = dchisq_objs[-1].GetYaxis().GetXmax()
    breactors = ROOT.TBox(X1, Y1, X2, Y2)
    breactors.SetFillColor(ROOT.kYellow+1)
    breactors.SetLineColor(ROOT.kYellow+1)
    if draw_reactor:
        breactors.Draw("same")

    critvals = [util.DeltaChisquare(conf, 2) for conf in conf_levels]    # The 2nd arg is the number of dimensions

    multigraph = ROOT.TMultiGraph()
    if len(dchisq_objs) == 3 and not "nova" in pdf_filename:
        line_color = 2
    else:
        line_color = 1

    i_dchisq_critval = 0
    for i_dchisq, dchisq_obj in enumerate(dchisq_objs):
        for i_critval, critval in enumerate(critvals):
            # Get the contour(s) for each critical value one set at a time, so we can apply the correct plot styles and labels
            dchisq_obj.SetContour(1, np.array([critval]))   # The first argument is the number of critical values
            contours = get_contour(dchisq_obj)
    
            if len(labels) > 0:
                label = labels[i_dchisq]
                if len(conf_levels) > 1:
                    label += " - " + util.ConfLevelAsStringRooTex(conf_levels[i_critval])
                legend.AddEntry(contours[0], label, "l")    # Only add the first contour from this set of contours to the legend
            for contour in contours:
                # Apply the same colour/style to all contours of the same critical value from the same file
                if not colours or not styles:
                    contour.SetLineColor(line_color)
                    if len(contours) < 2:
                        # No 68% contour, set 90% CL style to solid
                        contour.SetLineStyle(1)
                else:
                    # Custom line colours and styles have been specified
                    contour.SetLineColor(colours[i_dchisq_critval])
                    contour.SetLineStyle(styles[i_dchisq_critval])
                if line_width:
                    contour.SetLineWidth(line_width)
                multigraph.Add(contour)
            i_dchisq_critval += 1
        line_color += 1
    multigraph.Draw()

    # Plot best-fit point marker
    bf_legend = ROOT.TLegend(bf_lxmin, bf_lymin, bf_lxmax , bf_lymax)
    bf_legend.SetFillStyle(0)
    bf_legend.SetBorderSize(0)
    bf_legend.SetTextFont(132)
    bf_legend.SetTextSize(text_size)
    if best_fit_labels and best_fit_styles:
        bf_points = []
        for i, dchisq_obj in enumerate(dchisq_objs):
            if best_fit_styles[i] != 0:
                minx = ROOT.Long(0)
                miny = ROOT.Long(0)
                minz = ROOT.Long(0)
                dchisq_obj.GetMinimumBin(minx, miny, minz)
                bf_point = ROOT.TMarker(dchisq_obj.GetXaxis().GetBinCenter(minx), dchisq_obj.GetYaxis().GetBinCenter(miny), best_fit_styles[i])
                bf_point.SetMarkerColor(colours[i*len(critvals)])
                if marker_size:
                    bf_point.SetMarkerSize(marker_size)
                bf_points.append(copy.deepcopy(bf_point))   # Each iteration of the loop would overwrite the last if we didn't copy
                bf_points[len(bf_points)-1].Draw("same")
                if best_fit_labels[i].replace(" ", "") != "":
                    bf_legend.AddEntry(bf_points[len(bf_points)-1], best_fit_labels[i], "p")

    if draw_reactor:
        bf_legend.AddEntry(breactors, breactors_label, "f")

    if bf_legend.GetNColumns() > 0:
        bf_legend.Draw()

    if "atmos" in pdf_filename:
        input = ROOT.TMarker(theta23, dm23, 23)
        input.SetMarkerColor(12)
        input.Draw("SAME")
        legend.AddEntry(input, "Truth", "p")
    legend.Draw()

    if prelim:
        title = ROOT.TLatex(0.6, 0.96, "T2K Run 1-9 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")
    
    multigraph.SetName("contours_mulitgraph")
    root_filename = pdf_filename + "_contours.root"
    output_file = ROOT.TFile(root_filename,"recreate")
    print "ROOT file", root_filename, "has been created"
    multigraph.Write()

    canvas.Update()
    canvas.Write()
    canvas.Print(pdf_filename + ".pdf")
    output_file.Close()
    del output_file

    

#________________________________________________________________________________________________________
def get_contour(dchisq_obj):
    canvas = ROOT.TCanvas("c", "c", 0, 0, 700, 500)
    dchisq_obj.Draw("cont3same")
    dchisq_obj.Draw("CONT LIST same")
    canvas.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    graphs = []
    lst = contours.At(0)
    next = ROOT.TIter(lst)
    for g in next:
        g.SetLineStyle(2)
        graphs.append(g.Clone())
    if contours.GetEntriesFast() > 1:    
        lst = contours.At(1)
        next = ROOT.TIter(lst)
        for g in next:
            g.SetLineStyle(1)
            graphs.append(g.Clone())
    return graphs

#________________________________________________________________________________________________________
def get_dchisq_histos(filenames, chisq_calc_type, dchisq_objects):
    potential_chisq_hist_names = ["chi2", "h_chisq", "chisq_map_2D"]
    ndim = 0

    # check the chisq calc type
    if not (chisq_calc_type == "MHpair" or chisq_calc_type == "global" or chisq_calc_type == "local"):
        print "ERROR: chisq calculation type [", chisq_calc_type, "] not understood. Exiting."
        sys.exit()

    # Read in the chisq histos
    for filename in filenames:
        file = ROOT.TFile(filename, "READ")
        histo = None
        for hname in potential_chisq_hist_names:
            if file.GetListOfKeys().Contains(hname):
                histo = file.Get(hname)
                histo.SetDirectory(0)

                # Check the dimensionality of the histogram
                ndim_this_histo = get_ndims(histo)
                if ndim == 0:
                    ndim = ndim_this_histo
                elif ndim_this_histo != ndim:
                    print "Error: Not all contours have the same dimensionality"
                    sys.exit()

                break

        if chisq_calc_type == "MHpair":
            # Get the chisq hist for the other hierarchy and convert both to Delta chisq
            filename_other_MH = ""
            if "normal" in filename:
                filename_other_MH = filename.replace("normal", "inverted")
            else:
                filename_other_MH = filename.replace("inverted", "normal")
            file_other_MH = ROOT.TFile(filename_other_MH, "READ")
            histo_other_MH = None
            for hname in potential_chisq_hist_names:
                if file_other_MH.GetListOfKeys().Contains(hname):
                    histo_other_MH = file_other_MH.Get(hname)
                    histo_other_MH.SetDirectory(0)

                    # Check the dimensionality of the histogram
                    ndim_this_histo = get_ndims(histo)
                    if ndim == 0:
                        ndim = ndim_this_histo
                    elif ndim_this_histo != ndim:
                        print "Error: Not all contours have the same dimensionality"
                        sys.exit()

                    break

            chisq_min_pair = get_chisq_min([histo, histo_other_MH])
            convert_to_delta_chisq([histo, histo_other_MH], ndim, chisq_min_pair)
            dchisq_objects.append(histo)
            dchisq_objects.append(histo_other_MH)
        else:
            # Just append to the list as chisq hists, then convert them to Delta chisq later
            dchisq_objects.append(histo)

    # Convert the chisq hists into Delta chisq hists
    if chisq_calc_type == "global":
        convert_to_delta_chisq(dchisq_objects, ndim)
    elif chisq_calc_type == "local":
        chisq_min_all = get_chisq_min(dchisq_objects)
        convert_to_delta_chisq(dchisq_objects, ndim, chisq_min_all)

    return ndim

#________________________________________________________________________________________________________
def get_ndims(histo):
    if type(histo) == type(ROOT.TH1D()):
        return 1
    elif type(histo) == type(ROOT.TH2D()):
        return 2
    else:
        print "ERROR: Could not get the number of histogram dimensions. Exiting."
        sys.exit()

#________________________________________________________________________________________________________
def get_chisq_min(chisq_hists):
    local_min_vals = [h.GetMinimum() for h in chisq_hists]
    chisq_min = min(local_min_vals)
    if chisq_min < 0:
        print "ERROR: Found min chisq < 0! Exiting."
        sys.exit()
    return chisq_min

#________________________________________________________________________________________________________
def convert_to_delta_chisq(chisq_histos, ndim, chisq_min = -1):
    for histo in chisq_histos:
        use_chisq_min = chisq_min
        if use_chisq_min == -1:
            use_chisq_min = histo.GetMinimum()
        for i in range(1, histo.GetNbinsX() + 1):
            if ndim == 2:
                for j in range(1, histo.GetNbinsY() + 1):
                    histo.SetBinContent(i, j, histo.GetBinContent(i, j) - use_chisq_min)
            else:
                histo.SetBinContent(i, histo.GetBinContent(i) - use_chisq_min)

#________________________________________________________________________________________________________
# Search for the boundaries of the 1D allowed FC region.
# Returns [left boundary, right boundary].
# NOTE: This assumes the free parameter is deltacp and that the input Delta chisq
#       hist spans the full range of [-pi, +pi]
def calc_1d_fc_xrange(histo, fc_critval_tgraph):
    lim_x_min = None
    lim_x_max = None

    # This function assumes the free parameter is deltacp and that the input Delta chisq
    # hist spans the full range of [-pi, +pi]. Test this assumption.
    first_bin_centre = histo.GetBinCenter(1)
    last_bin_centre  = histo.GetBinCenter(histo.GetNbinsX())
    if ROOT.TMath.Abs(first_bin_centre + ROOT.TMath.Pi()) > conv.kASmallNum or \
       ROOT.TMath.Abs(last_bin_centre  - ROOT.TMath.Pi()) > conv.kASmallNum:
        print "ERROR: calc_1d_fc_xrange() assumes the input Delta chisq hist spans the full range of [-pi, +pi].", \
              "This is not the case. Instead, the hist bins centres span the range [", first_bin_centre, \
              ",", last_bin_centre, "]"
        sys.exit()

    # Shift the input Delta chisq hist by the critical values in each bin centre.
    # This allows us to calculate the range by comparing to a constant (shifted) critical value of 0.
    # This is equivalent to comparing to a variable critical value, but is easier to handle
    critval_shifted = 0.0
    hist_shifted = ROOT.TH1D(histo)
    for ibin in range(1, hist_shifted.GetNbinsX()+1):
        bincenter = histo.GetBinCenter(ibin)
        bincont = histo.GetBinContent(ibin)
        bincont_new = bincont - fc_critval_tgraph.Eval(bincenter)
        hist_shifted.SetBinContent(ibin, bincont_new)

    bin_bf = hist_shifted.GetMinimumBin()
    val_bf = hist_shifted.GetBinCenter(bin_bf)
    
    # Bins before (b) and after (a) the intersection with critical value
    bin_left_search_a  = conv.kBadResult
    bin_left_search_b  = conv.kBadResult
    bin_right_search_a = conv.kBadResult
    bin_right_search_b = conv.kBadResult
    
    # Search left of the chisq min for the boundary of the allowed FC region,
    # wrapping around the deltacp = +/- pi boundary, as needed
    found = False
    for ibin in range(bin_bf, 0, -1):
        bincont = hist_shifted.GetBinContent(ibin)
        if bincont > critval_shifted:
            # bins before (b) and after (a) DeltaChisq = critical value (following x direction)
            bin_left_search_b = ibin
            bin_left_search_a = bin_left_search_b+1
            found = True
            break
    if not found:
        # The interval hit the boundary, wrap around
        for ibin in range(hist_shifted.GetNbinsX(), bin_bf, -1):
            bincont = hist_shifted.GetBinContent(ibin)
            if bincont > critval_shifted:
                # bins before (b) and after (a) DeltaChisq = critical value (following x direction)
                bin_left_search_b = ibin
                bin_left_search_a = bin_left_search_b+1
                found = True
                break
    if not found:
        # If no bin is found it means the DeltaChisq = critical value bin is outside the physical region.
        # If so, set the left extreme to 1 (first bin of the histogram)
        print "WARNING: Could not form an allowed FC range when looking left of the chisq min"
        bin_left_search_a = 1
        bin_left_search_b = 1

    # Search right of the chisq min for the boundary of the allowed FC region,
    # wrapping around the deltacp = +/- pi boundary, as needed
    found = False
    for ibin in range(bin_bf, hist_shifted.GetNbinsX()+1):
        bincont = hist_shifted.GetBinContent(ibin)
        if bincont > critval_shifted:
            # bins before (b) and after (a) DeltaChisq = critical value
            bin_right_search_a = ibin
            bin_right_search_b = bin_right_search_a-1
            found = True
            break
    if not found:
        # The interval hit the boundary, wrap around
        for ibin in range(1, bin_bf):
            bincont = hist_shifted.GetBinContent(ibin)
            if bincont > critval_shifted:
                # bins before (b) and after (a) DeltaChisq = critical value (following x direction)
                bin_right_search_b = ibin
                bin_right_search_a = bin_right_search_b+1
                found = True
                break
    if not found:
        # If no bin is found it means the DeltaChisq = critical value bin is outside the physical region.
        # If so, set the right extreme to Nbins (last bin of the histogram)
        print "WARNING: Could not form an allowed FC range when looking right of the chisq min"
        bin_right_search_a = hist_shifted.GetNbinsX()
        bin_right_search_b = hist_shifted.GetNbinsX()

    # Extrapolation for the left extreme.
    # Extrapolate between bins to get the x-axis value corresponding to DeltaChisq = critical value.
    # Do it only if DeltaChisq = critical value bin is within the dchisq histogram
    if bin_left_search_a > 1:
        X1 = hist_shifted.GetXaxis().GetBinCenter(bin_left_search_b)
        X2 = hist_shifted.GetXaxis().GetBinCenter(bin_left_search_a)
        Y1 = hist_shifted.GetBinContent(bin_left_search_b)
        Y2 = hist_shifted.GetBinContent(bin_left_search_a)
        lim_x_min = get_interpolation_x(X1,Y1,X2,Y2,critval_shifted)
    else:
        lim_x_min = hist_shifted.GetBinCenter(bin_left_search_a)
    errneg = ROOT.TMath.Abs(val_bf - lim_x_min)

    # Extrapolation for the right extreme.
    # Extrapolate between bins to get the x-axis value corresponding to DeltaChisq = critical value.
    # Do it only if DeltaChisq = critical value bin is within the dchisq histogram
    if bin_right_search_b < hist_shifted.GetNbinsX():
        X1 = hist_shifted.GetXaxis().GetBinCenter(bin_right_search_b)
        X2 = hist_shifted.GetXaxis().GetBinCenter(bin_right_search_a)
        Y1 = hist_shifted.GetBinContent(bin_right_search_b)
        Y2 = hist_shifted.GetBinContent(bin_right_search_a)
        lim_x_max = get_interpolation_x(X1,Y1,X2,Y2,critval_shifted)
    else:
        lim_x_max = hist_shifted.GetBinCenter(bin_right_search_b)
    errpos = ROOT.TMath.Abs(val_bf - lim_x_max)

    # TODO: We use this function for multiple hists, so input a name so that the following printout is clearer
    print "\nDelta chisq min +/- allowed FC range =", val_bf, "+", errpos, "-", errneg
    print "\nThe FC confidence interval is: [", val_bf - errneg, ";", val_bf + errpos, "]"

    return [lim_x_min, lim_x_max]

#________________________________________________________________________________________________________
# Linearly interpolate between two points to get the x-value corresponding to an input y-value
def get_interpolation_x(X1, Y1, X2, Y2, Yin):
    if Y1 <= Y2:
        if Yin < Y1 or Yin > Y2:
            print "WARNING: Trying to linearly interpolate y =", Yin, \
                  "outside of the valid y range: [", Y1, ",", Y2, "]"
    else:
        if Yin < Y2 or Yin > Y1:
            print "WARNING: Trying to linearly interpolate y =", Yin, \
                "outside of the valid y range: [", Y2, ",", Y1, "]"
    slope = (Y2-Y1)/(X2-X1)
    Xout = (Yin-Y1)/slope + X1
    return Xout

#________________________________________________________________________________________________________
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename"       , dest="in_filenames"     , nargs="+"              , help="Input filename(s) for contour(s).")
    parser.add_argument("-l", "--label"          , dest="labels"           , nargs="+"              , help="Legend label(s) for contours(s). This should have either 0 length or the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument("-c", "--colour"         , dest="colours"          , nargs="+"              , help="Custom line colour(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the colours is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument("-s", "--style"          , dest="styles"           , nargs="+"              , help="Custom line style(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename  (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the styles is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--best-fit-label" , dest="best_fit_labels"  , nargs="+"              , help="Labels for best-fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--best-fit-style" , dest="best_fit_styles"  , nargs="+"              , help="Marker styles for best fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--conf-levels"    , dest="conf_levels_str"  , nargs="+"              , help="The confidence levels to plot const chisq lines/contours for. See ConfLevel::AsString() for the choices")
    parser.add_argument(      "--title-x"        , dest="title_x"                                   , help="Override default x-axis title")
    parser.add_argument(      "--title-y"        , dest="title_y"          , default="-2#Deltaln(L)", help="Override default y-axis title")
    parser.add_argument(      "--line-width"     , dest="lwidth"           , type=int               , help="Contour line width")
    parser.add_argument(      "--marker-size"    , dest="msize"            , type=float             , help="Best-fit marker size")
    parser.add_argument(      "--text-size"      , dest="tsize"            , type=float             , help="Legend text size")
    parser.add_argument(      "--title-size"     , dest="titlesize"        , type=float,default=0.05, help="Title text size")
    parser.add_argument("-o", "--output"         , dest="out_filename"                              , help="Output PDF filename")
    parser.add_argument(      "--chisq-calc"     , dest="chisq_calc"                                , help="If \'global\', the min. Dchisq of each contour is set to 0. If \'local\', Dchisq is calculated relative to the min. of all contours. If \'MHpair\', the contour for the opposite MH is found and the Dchisq of each contour set is relative to the min. of pair of NH/IH contours")
    parser.add_argument(      "--fc-files"       , dest="fc_files"         , nargs="+"              , help="ROOT files containing TGraphs of the FC critical values as a function of the parameter of interest and the mass hierarchy. For 1D fits only. If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--fc-graph-names" , dest="fc_graph_names"   , nargs="+"              , help="Names of the TGraphs containing the FC critical values. If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--fc-line-styles" , dest="fc_line_styles"   , nargs="+"              , help="Line styles used for plotting the FC critical values with \'--fc-draw-vals\'. If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--fc-line-colours", dest="fc_line_colours"  , nargs="+"              , help="Line colours used for plotting the FC critical values with \'--fc-draw-vals\'. If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--fc-draw-vals"   , dest="fc_draw_vals"     , action="store_true"    , help="If true, draws the FC critical values on top of the Delta chisq hists")
    parser.add_argument(      "--fc-draw-boxes"  , dest="fc_draw_boxes"    , action="store_true"    , help="If true, draws shaded boxes representing the FC allowed region on top of the Delta chisq hists")
    parser.add_argument(      "--t2k"            , dest="is_t2k_only"      , action="store_true"    )
    parser.add_argument(      "--prelim"         , dest="prelim"           , action="store_true"    )
    parser.add_argument(      "--theta23"        , dest="theta23"          , type=float             , help="Input value of theta23.")
    parser.add_argument(      "--dm23"           , dest="dm23"             , type=float             , help="Input value of dm23.")
    parser.add_argument(      "--draw-reactor"   , dest="draw_reactor"     , action="store_true"    , help="Draw a vertical band indicating the PDG 2018 reactor constraint on sinsq(theta13) used in the T2K oscillation analysis")
    parser.add_argument(      "--force-max-y"    , dest="force_max_y"      , type=float             , help="Set the maximum y-axis value in 1D plots to this.")
    parser.add_argument(      "--max-ndigits-y"  , dest="max_ndigits_y"    , type=int, default=3    , help="Max number of digits on the y-axis labels.")
    parser.add_argument(      "--leg-coords"     , dest="leg_coords_str"   , nargs="+"              , help="Legend corner coordinates [xmin, ymin, xmax, ymax].")
    parser.add_argument(      "--bf-leg-coords"  , dest="bf_leg_coords_str", nargs="+"              , help="Best-fit marker legend corner coordinates [xmin, ymin, xmax, ymax] (for 2D plots only).")
    args = parser.parse_args()

    # Check the validity of the input arguments

    labels = []
    if args.labels:
        labels = [i for i in args.labels]
        if not len(args.in_filenames) == len(labels) and args.chisq_calc != "MHpair":
            print "Error: The number of filenames and labels should match"
            sys.exit()
        elif not 2 * len(args.in_filenames) == len(labels) and args.chisq_calc == "MHpair":
            print "Error: The number of labels should be twice the number of filenames"
            print "filenames =", args.in_filenames
            print "labels =", labels
            sys.exit()


    dchisq_objects = []
    num_var = get_dchisq_histos(args.in_filenames, args.chisq_calc, dchisq_objects)
    print "Number of chisq hists to plot =", len(dchisq_objects)

    conf_levels = []
    if not args.conf_levels_str:
        conf_levels = [conv.kCL1Sigma, conv.kCL90]
        print "Setting default confidence levels to [1sigma, 90%]"
    else:
        conf_levels = [util.ConfLevelFromString(conf_str) for conf_str in args.conf_levels_str]

    colours = []
    if args.colours:
        colours = [int(i) for i in args.colours]
        if num_var == 1:
            if not len(colours) == len(args.in_filenames) and args.chisq_calc != "MHpair":
                print "Error: The number of filenames and colours should match"
                sys.exit()
            elif not len(colours) == 2 * len(args.in_filenames) and args.chisq_calc == "MHpair":
                print "Error: The number of colours should be twice the number of filenames"
                print "filenames =", args.in_filenames
                print "colours =", colours
                sys.exit()
        elif num_var == 2:
            if not len(colours) == len(args.in_filenames) * len(conf_levels) and args.chisq_calc != "MHpair":
                print "Error: The number of colours should = number of filenames * number of confidence levels"
                sys.exit()
            elif not len(colours) == 2 * len(args.in_filenames) * len(conf_levels) and args.chisq_calc == "MHpair":
                print "Error: The number of colours should = 2 * number of filenames * number of confidence levels"
                print "filenames =", args.in_filenames
                print "confidence levels =", args.conf_levels_str
                print "colours =", colours
                sys.exit()
        else:
            print "Error: Only 1D or 2D contours are supported."
            sys.exit()


    styles = []
    if args.styles:
        styles = [int(i) for i in args.styles]
        if num_var == 1:
            if not len(styles) == len(args.in_filenames) and args.chisq_calc != "MHpair":
                print "Error: The number of filenames and styles should match"
                sys.exit()
            elif not len(styles) == 2 * len(args.in_filenames) and args.chisq_calc == "MHpair":
                print "Error: The number of styles should be twice the number of filenames"
                print "filenames =", args.in_filenames
                print "styles =", styles
                sys.exit()
        elif num_var == 2:
            if not len(styles) == len(args.in_filenames) * len(conf_levels) and args.chisq_calc != "MHpair":
                print "Error: The number of styles should = number of filenames * number of confidence levels"
                sys.exit()
            elif not len(styles) == 2 * len(args.in_filenames) * len(conf_levels) and args.chisq_calc == "MHpair":
                print "Error: The number of styles should = 2 * number of filenames * number of confidence levels"
                print "filenames =", args.in_filenames
                print "confidence levels =", args.conf_levels_str
                print "styles =", styles
                sys.exit()
        else:
            print "Error: Only 1D or 2D contours are supported."
            sys.exit()
        
        
    best_fit_styles = []
    if args.best_fit_styles:
        best_fit_styles = [int(i) for i in args.best_fit_styles]
        if num_var == 2 and len(best_fit_styles) != len(args.in_filenames):
            print "ERROR: The number of best-fit marker styles must equal the number of input files"
            print "filenames =", args.in_filenames
            print "best-fit styles =", args.best_fit_styles
            sys.exit()

    leg_coords = []
    if args.leg_coords_str:
        leg_coords = [float(i) for i in args.leg_coords_str]
        if len(leg_coords) != 4:
            print "ERROR: The number of legend coordinates must equal 4"

    bf_leg_coords = []
    if args.bf_leg_coords_str:
        bf_leg_coords = [float(i) for i in args.bf_leg_coords_str]
        if len(leg_coords) != 4:
            print "ERROR: The number of best-fit point legend coordinates must equal 4"

    fc_line_styles  = []
    fc_line_colours = []
    if args.fc_files:
        if args.chisq_calc != "MHpair":
            if len(args.fc_files) != len(args.in_filenames):
                print "ERROR: The number of input FC files must match the number of input chisq hist files."
                sys.exit()
        else:
            if len(args.fc_files) != 2 * len(args.in_filenames):
                print "ERROR: The number of input FC files must be double the number of input chisq hist files."
                sys.exit()
        if len(args.fc_graph_names) != len(args.fc_files):
            print "ERROR: The number of input FC TGraph names must match the number of input FC files."
            sys.exit()
        if len(args.fc_line_styles) != len(args.fc_files):
            print "ERROR: The number of input FC line styles must match the number of input FC files."
            sys.exit()
        if len(args.fc_line_colours) != len(args.fc_files):
            print "ERROR: The number of input FC line colours must match the number of input FC files."
            sys.exit()
        if len(conf_levels) != 1:
            print "ERROR: Tried to plot 1D Feldman-Cousins with != 1 input confidence levels. Plotting FC is only supported for a single confidence level."
            sys.exit()
        if num_var != 1:
            print "ERROR: Only 1D Feldman-Cousins is supported."
            sys.exit()
        if "cp" not in args.title_x.lower():
            print "ERROR: Currently plotting Feldman-Cousins is only supported for deltacp"
            sys.exit()
        fc_line_styles  = [int(i) for i in args.fc_line_styles ]
        fc_line_colours = [int(i) for i in args.fc_line_colours]

    ROOT.TGaxis.SetMaxDigits(args.max_ndigits_y)

    is_inverted = False     # TODO: deal with this. Should be specified individually?

    if num_var == 1:
        plot_contours_1d(dchisq_objects, labels, args.out_filename, args.is_t2k_only, colours, styles, args.prelim,
                         args.title_x, args.title_y, args.lwidth, args.tsize, conf_levels, leg_coords, args.titlesize,
                         args.fc_files, args.fc_graph_names, fc_line_styles, fc_line_colours, args.fc_draw_vals,
                         args.fc_draw_boxes, is_inverted, args.force_max_y)
    elif num_var == 2:
        plot_contours_2d(dchisq_objects, labels, args.theta23, args.dm23, args.out_filename, colours, styles,
                         best_fit_styles, args.best_fit_labels, args.prelim, args.title_x, args.title_y, args.lwidth,
                         args.msize, args.tsize, conf_levels, leg_coords, bf_leg_coords, args.titlesize, args.draw_reactor)
    else:
        print "Error: Not all contours have the same dimensionality"
        sys.exit()