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

_hatch_fill_style   = 1001    # For filling in FC intervals
_hatch_fill_opacity = 0.15    # For filling in FC intervals

#________________________________________________________________________________________________________
def plot_contours_1d(dchisq_objects, labels, pdf_filename, is_t2k_only, colours, styles, prelim, title_x, title_y,
                     line_width, text_size, conf_levels, leg_coords, titlesize,
                     fc_files, fc_graph_pattern, fc_line_styles, fc_line_colours, fc_draw_critvals,
                     fc_draw_boxes, ptheta_objects = None, mach3_objects = None,
                     x_range = None, y_range = None, chisq_calc = "global", hier_leg_coords = None):
    max_val = 0.0
    min_val = 999
    for dchisq_obj in dchisq_objects:
        max_val = max(max_val, dchisq_obj.GetMaximum())
        min_val = min(min_val, dchisq_obj.GetMinimum())
    if ptheta_objects:
        for dchisq_obj in ptheta_objects:
            max_val = max(max_val, dchisq_obj.GetMaximum())
            min_val = min(min_val, dchisq_obj.GetMinimum())
    if mach3_objects:
        for dchisq_obj in mach3_objects:
            max_val = max(max_val, dchisq_obj.GetMaximum())
            min_val = min(min_val, dchisq_obj.GetMinimum())
    max_val = np.ceil(max_val)
    max_val = min(max_val, 25.0)    # Not meaningful above 5 sigma, so truncate the axis
    if len(y_range) == 2:
        max_val = y_range[1]
    min_val = np.floor(min_val)
    for dchisq_obj in dchisq_objects:
        if len(x_range) == 2:
            dchisq_obj.GetXaxis().SetRangeUser(x_range[0], x_range[1])
        if len(y_range) == 2:
            dchisq_obj.GetYaxis().SetRangeUser(y_range[0], y_range[1])
        else:
            dchisq_obj.GetYaxis().SetRangeUser(min_val, max_val)

    # When drawing group comparisons, the draw order is 1) VALOR, 2) p-theta, & 3) MaCh3, so make the VALOR and p-theta
    # lines and markers thicker so they can be seen more easily
    mach3_line_width   = line_width
    ptheta_line_width  = line_width
    if mach3_objects:
        if line_width:
            line_width  += 1
    if mach3_objects:
        if line_width:
            line_width        += 1
            ptheta_line_width += 1

    canvas = ROOT.TCanvas("canvas1d", "canvas1d", 700, 500)

    breactors = None
    breactors_label = None
    if "mixing13" in pdf_filename:
        sinsqth13_nom = util.DoubleAngleToSingleAngle(0.0853)
        sinsqth13_err = util.DoubleAngleSigmaToSingleAngleSigma(0.0853, 0.0027)
        breactors_label = "PDG 2019"
        X1 = sinsqth13_nom - sinsqth13_err
        X2 = sinsqth13_nom + sinsqth13_err
        breactors = ROOT.TBox(X1,min_val,X2,max_val)
        breactors.SetFillColor(ROOT.kYellow+1)
        dchisq_objects[0].Draw("axis")
        breactors.Draw("same")

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
        line_color += 1

    num_leg_entries = len(dchisq_objects)
    if ptheta_objects: num_leg_entries += len(ptheta_objects)
    if mach3_objects : num_leg_entries += len(mach3_objects)

    if not text_size:
        text_size = 0.03
    lymax = 0.86
    lymin = lymax - (1 + num_leg_entries) * text_size
    lxmin = 0.15
    lxmax = 0.3

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
        i_label = 0
        for i, dchisq_obj in enumerate(dchisq_objects):
            if chisq_calc == "MHpair" and i % 2 != 0:
                # If we're plotting both mass hierarchies, then each input file has a pair corresponding to
                # the other hierarchy. So only add the first file to the legend, and have a separate legend
                # for the hierarchy
                continue
            if labels[i_label] != "":
                legend.AddEntry(dchisq_obj, labels[i_label], "l")
                i_label += 1

    legend_hierarchy = None
    if chisq_calc == "MHpair" and hier_leg_coords:
        legend_hierarchy = ROOT.TLegend(hier_leg_coords[0], hier_leg_coords[1], hier_leg_coords[2], hier_leg_coords[3])
        legend_hierarchy.SetFillStyle(0)
        legend_hierarchy.SetBorderSize(0)
        legend_hierarchy.SetTextFont(132)
        legend_hierarchy.SetTextSize(text_size)
        legend_hierarchy.AddEntry(dchisq_objects[0], "Normal"  , "l")
        legend_hierarchy.AddEntry(dchisq_objects[1], "Inverted", "l")

    # If we have contours from other groups, format and draw them here
    if ptheta_objects:
        for i, obj in enumerate(ptheta_objects):
            obj.SetLineWidth(ptheta_line_width)
            obj.Draw("l,same")
            if obj.GetTitle() != "" and obj.GetTitle() != "-":
                if chisq_calc == "MHpair" and i % 2 != 0:
                    continue
                legend.AddEntry(obj, obj.GetTitle(), "l")
    if mach3_objects:
        for i, obj in enumerate(mach3_objects):
            obj.SetLineWidth(mach3_line_width)
            obj.Draw("l,same")
            if obj.GetTitle() != "" and obj.GetTitle() != "-":
                if chisq_calc == "MHpair" and i % 2 != 0:
                    continue
                legend.AddEntry(obj, obj.GetTitle(), "l")

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

        print_best_fit_and_range_1d_const_chisq(dchisq_objects, conf_levels)

    else:
        # Draw boxes to represent the Feldman-Cousins allowed regions
        fc_tgraphs        = []
        hatch_boxes = []
        for i, dchisq_obj in enumerate(dchisq_objects):
            if fc_graph_pattern[i] == "None" or fc_graph_pattern[i] == "none" or fc_graph_pattern[i] == "":
                # Not every contour needs to have an associated FC hatch box or graph of critical values
                continue
            for j, conf_level in enumerate(conf_levels):
                fc_tfile = ROOT.TFile(fc_files[i])
                fc_graph_name = (fc_graph_pattern[i]).replace("[CL]",util.ConfLevelAsString(conf_level))
                fc_critval_tgraph = fc_tfile.Get(fc_graph_name)
                if not fc_critval_tgraph:
                    # Try the legacy analysis naming convention before we give up
                    if "pct" in fc_graph_name:
                        fc_critval_tgraph = fc_tfile.Get(fc_graph_name.replace("pct","CL"))
                    if not fc_critval_tgraph:
                        print "ERROR: Could not find a TGraph named \'{}\'".format(fc_graph_name), \
                              " or \'{}\'".format(fc_graph_name.replace("pct","CL")), \
                              "in the file", fc_files[i]
                        sys.exit()

                fc_critval_tgraph.SetLineStyle(fc_line_styles [i*(len(conf_levels)) +j])
                fc_critval_tgraph.SetLineColor(fc_line_colours[i*(len(conf_levels)) +j])
                fc_tgraph_leg_label  = "FC" + " " + title_y + "_{crit.} " # Have to split this up so that the RooTex subscript/superscript doesn't mess up format()
                fc_tgraph_leg_label += "({} CL)".format(util.ConfLevelAsStringRooTex(conf_level))
                fc_critval_tgraph.SetName(fc_tgraph_leg_label)
                fc_tgraphs.append(fc_critval_tgraph)

                if fc_draw_boxes:
                    # This only gets done for a single confidence level. This condition is is enforced
                    # when processing command line args.
                    [left_boundary, right_boundary] = calc_1d_fc_xrange(dchisq_obj, fc_critval_tgraph)
                    hist_left_edge  = dchisq_obj.GetXaxis().GetBinLowEdge(1)
                    hist_right_edge = dchisq_obj.GetXaxis().GetBinUpEdge(dchisq_obj.GetNbinsX())
                    if not (left_boundary == conv.kBadResult and right_boundary == conv.kBadResult):

                        bin_bf = dchisq_obj.GetMinimumBin()
                        val_bf = dchisq_obj.GetBinCenter(bin_bf)
                        errneg = str(ROOT.TMath.Abs(val_bf - left_boundary))
                        left_boundary_str = str(left_boundary)
                        if left_boundary == conv.kBadResult:
                            errneg = "not found"
                            left_boundary_str = errneg
                        errpos = str(ROOT.TMath.Abs(val_bf - right_boundary))
                        right_boundary_str = str(right_boundary)
                        if right_boundary == conv.kBadResult:
                            errpos = "not found"
                            right_boundary_str = errpos

                        print "======================================="
                        print "For the fit from file \'{}\'".format(dchisq_obj.GetTitle()), \
                              "and the", util.ConfLevelAsString(conf_level), "CL:"
                        print "\nDelta chisq min +/- allowed FC range =", val_bf, "+", errpos, "-", errneg
                        print "\nThe FC confidence interval is: [", left_boundary_str, ";", right_boundary_str, "]"
                        print "======================================="

                        if left_boundary == conv.kBadResult:
                            left_boundary = hist_left_edge
                        if right_boundary == conv.kBadResult:
                            right_boundary = hist_right_edge
                        hatch_box_max_y = 0.66 * max_val
                        if left_boundary > right_boundary:
                            hatch_interval  = ROOT.TBox(left_boundary , 0, hist_right_edge, hatch_box_max_y)
                            hatch_interval.SetLineWidth(0)
                            hatch_interval.SetFillStyle(_hatch_fill_style)
                            hatch_interval.SetFillColorAlpha(dchisq_obj.GetLineColor(), _hatch_fill_opacity)
                            hatch_boxes.append(hatch_interval)
                            hatch_interval.Draw("same")

                            hatch_interval2 = ROOT.TBox(hist_left_edge, 0, right_boundary , hatch_box_max_y)
                            hatch_interval2.SetLineWidth(0)
                            hatch_interval2.SetFillStyle(_hatch_fill_style)
                            hatch_interval2.SetFillColorAlpha(dchisq_obj.GetLineColor(), _hatch_fill_opacity)
                            hatch_boxes.append(hatch_interval2)
                            hatch_interval2.Draw("same")
                        else:
                            hatch_interval  = ROOT.TBox(left_boundary , 0, right_boundary , hatch_box_max_y)
                            hatch_interval.SetLineWidth(0)
                            hatch_interval.SetFillStyle(_hatch_fill_style)
                            hatch_interval.SetFillColorAlpha(dchisq_obj.GetLineColor(), _hatch_fill_opacity)
                            hatch_boxes.append(hatch_interval)
                            hatch_interval.Draw("same")
                fc_tfile.Close()

        # Draw the FC critical values in a separate loop to the FC boxes to ensure that the critical value
        # lines are always in front of the boxes.
        if fc_draw_critvals:
            for i, fc_tgraph in enumerate(fc_tgraphs):
                fc_tgraph.Draw("L,same")
                legend.AddEntry(fc_tgraph, fc_tgraph.GetName(), "l")

        if fc_draw_boxes:
            box_name = "FC "+util.ConfLevelAsStringRooTex(conf_levels[0])+" conf. region"
            legend.AddEntry(hatch_boxes[0],box_name , "f")

    if prelim:
        title = ROOT.TLatex(0.5, 0.96, "T2K Run 1-10 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")

    if "mixing13" in pdf_filename:
        legend.AddEntry(breactors, breactors_label, "f")

    legend.Draw()
    if legend_hierarchy:
        legend_hierarchy.Draw()

    canvas.RedrawAxis()     # to stop objects from obscuring the axes
    canvas.Print(pdf_filename + ".pdf")
    canvas.Print(pdf_filename + ".root")

#________________________________________________________________________________________________________
def plot_contours_2d(dchisq_objs, labels, theta23, dm23, pdf_filename, colours, styles, best_fit_styles,
                     best_fit_labels, prelim, title_x, title_y, line_width, marker_size, text_size,
                     conf_levels, leg_coords, leg2_coords, titlesize, draw_reactor=False,
                     ptheta_objects = None, mach3_objects = None, x_range = [], y_range = []):

    lmargin = 0.13
    rmargin = 0.08
    tmargin = 0.05
    bmargin = 0.13

    # When drawing group comparisons, the draw order is 1) VALOR, 2) MaCh3, & 3) p-theta, so make the VALOR and MaCh3
    # lines and markers thicker so they can be seen more easily
    mach3_line_width   = line_width
    ptheta_line_width  = line_width
    mach3_marker_size  = marker_size
    ptheta_marker_size = marker_size
    if mach3_objects:
        if line_width:
            line_width  += 1
        if marker_size:
            marker_size += 1
    if ptheta_objects:
        if line_width:
            line_width        += 1
            mach3_line_width  += 1
        if marker_size:
            marker_size       += 1
            mach3_marker_size += 1


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
    if len(leg2_coords) != 0:
        # Use user-specified legend coordinates
        bf_lxmin = leg2_coords[0]
        bf_lymin = leg2_coords[1]
        bf_lxmax = leg2_coords[2]
        bf_lymax = leg2_coords[3]

    # Legend showing the styles of each separate fit (e.g. NH/IH or VALOR/MaCh3/ptheta etc.)
    legend = ROOT.TLegend(use_lxmin, lymin, lxmax, lymax)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextFont(132)
    legend.SetTextSize(text_size)

    # legend showing the styles of each conf level and the best fit point
    conf_legend = ROOT.TLegend(bf_lxmin, bf_lymin, bf_lxmax , bf_lymax)
    conf_legend.SetFillStyle(0)
    conf_legend.SetBorderSize(0)
    conf_legend.SetTextFont(132)
    conf_legend.SetTextSize(text_size)

    if title_x:
        dchisq_objs[-1].GetXaxis().SetTitle(title_x)
        dchisq_objs[-1].GetXaxis().SetTitleSize(titlesize)
    if title_y:
        dchisq_objs[-1].GetYaxis().SetTitle(title_y)
        dchisq_objs[-1].GetYaxis().SetTitleSize(titlesize)

    dchisq_objs[-1].GetYaxis().SetTitleOffset(0) #1.3
    dchisq_objs[-1].GetXaxis().SetTitleOffset(0)  #0.9
    if len(x_range) == 2:
        dchisq_objs[-1].GetXaxis().SetRangeUser(x_range[0], x_range[1])
    if len(y_range) == 2:
        dchisq_objs[-1].GetYaxis().SetRangeUser(y_range[0], y_range[1])
    dchisq_objs[-1].SetTitle("")
    dchisq_objs[-1].Draw("axis")

    # Draw the reactor constraint used in T2K, as a vertical yellow band
    sinsqth13_nom = util.DoubleAngleToSingleAngle(0.0853)
    sinsqth13_err = util.DoubleAngleSigmaToSingleAngleSigma(0.0853, 0.0027)
    breactors_label = "PDG 2019"
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

    # Add contours
    add_2d_contours(multigraph, legend, conf_legend, dchisq_objects, critvals, labels, conf_levels, colours, styles, line_width, line_color)
    if mach3_objects:
        add_2d_contours(multigraph, legend, conf_legend, mach3_objects, critvals, None, conf_levels, mach3_colours, mach3_styles, mach3_line_width, line_color)
    if ptheta_objects:
        add_2d_contours(multigraph, legend, conf_legend, ptheta_objects, critvals, None, conf_levels, ptheta_colours, ptheta_styles, ptheta_line_width, line_color)
    multigraph.Draw()

    # Plot best-fit point marker
    bf_points = []
    if best_fit_labels and best_fit_styles:
        draw_bf_markers(bf_points, conf_legend, dchisq_objs, best_fit_labels, best_fit_styles, colours, critvals, marker_size)
        if mach3_objects:
            draw_bf_markers(bf_points, conf_legend, mach3_objects, None, best_fit_styles, mach3_colours, None, mach3_marker_size)
        if ptheta_objects:
            draw_bf_markers(bf_points, conf_legend, ptheta_objects, None, best_fit_styles, ptheta_colours, None, ptheta_marker_size)

    if draw_reactor:
        conf_legend.AddEntry(breactors, breactors_label, "f")

    if conf_legend.GetNRows() > 0:
        conf_legend.Draw()

    if "atmos" in pdf_filename:
        input = ROOT.TMarker(theta23, dm23, 23)
        input.SetMarkerColor(12)
        input.Draw("SAME")
        legend.AddEntry(input, "Truth", "p")
    legend.Draw()

    if prelim:
        title = ROOT.TLatex(0.6, 0.96, "T2K Run 1-10 Preliminary")
        title.SetTextSize(0.05)
        title.SetNDC(True)
        title.Draw("same")

    multigraph.SetName("contours_mulitgraph")
    root_filename = pdf_filename + "_contours.root"
    output_file = ROOT.TFile(root_filename,"recreate")
    print "ROOT file", root_filename, "has been created"
    multigraph.Write()

    canvas.RedrawAxis()     # to stop objects from obscuring the axes
    canvas.Update()
    canvas.Write()
    canvas.Print(pdf_filename + ".pdf")
    output_file.Close()
    del output_file



#________________________________________________________________________________________________________
def add_2d_contours(multigraph, legend, conf_legend, dchisq_objs, critvals, leg_labels, conf_levels, colours, styles, line_width, line_color):
    already_have_contours = (type(dchisq_objs[0]) == list)  # i.e. for MaCh3 and p-theta

    # For the MaCh3 and p-theta TGraphs, we're storing the labels in the hist titles
    if not leg_labels:
        leg_labels = []
        for dchisq_obj in dchisq_objs:
            if not already_have_contours:
                leg_labels.append(dchisq_obj.GetTitle())
            elif type(dchisq_obj) == list:
                leg_labels.append(dchisq_obj[0].GetTitle())

    if not already_have_contours:
        i_dchisq_critval = 0
        for i_dchisq, dchisq_obj in enumerate(dchisq_objs):
            for i_critval, critval in enumerate(critvals):
                # Get the contour(s) for each critical value one set at a time, so we can apply the correct plot styles and labels
                dchisq_obj.SetContour(1, np.array([critval]))   # The first argument is the number of critical values
                contours = get_contour(dchisq_obj)

                if len(leg_labels) > 0 and i_critval == 0:
                    label = leg_labels[i_dchisq]
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

                if i_dchisq == 0 and len(conf_levels) > 1 and len(contours) > 0:
                    conf_legend.AddEntry(contours[0], util.ConfLevelAsStringRooTex(conf_levels[i_critval]) + " CL", "l")

                i_dchisq_critval += 1
            line_color += 1
    else:
        i_dchisq_critval = 0
        for i_cl, cl_objs in enumerate(dchisq_objs):
            if type(cl_objs) != list:
                # In dchisq_objs we've shoved both BF point markers and lists of contours, but here we only want to deal with the contours
                break
            if len(leg_labels) > 0 and i_cl == 0:
                label = leg_labels[i_cl]
                legend.AddEntry(cl_objs[0], label, "l")    # Only add the first contour from this set of contours to the legend

            for i_cont, contour in enumerate(cl_objs):
                # Apply the same colour/style to all contours of the same critical value from the same file
                if not colours or not styles:
                    contour.SetLineColor(line_color)
                    if len(cl_objs) < 2:
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
def draw_bf_markers(bf_points, conf_legend, dchisq_objs, best_fit_labels, best_fit_styles, colours, critvals, marker_size):
    have_marker_graphs = (type(dchisq_objs[0]) == list)

    if not have_marker_graphs:
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
                if i < len(best_fit_labels):
                    if best_fit_labels[i].replace(" ", "") != "":
                        conf_legend.AddEntry(bf_points[-1], best_fit_labels[i], "p")
    else:
        # TODO: This only allows us to get one set of BF markers for each of MaCh3 and p-theta, so need to get this in a better way. e.g. having their markers in a separate list would make more sense than this.
        marker_graph = dchisq_objs[-1]
        assert(marker_graph.GetN() == 1)
        x = ROOT.Double(-999)
        y = ROOT.Double(-999)
        marker_graph.GetPoint(0,x,y)
        bf_point = ROOT.TMarker(x, y, best_fit_styles[0])
        bf_point.SetMarkerColor(colours[0])
        if marker_size:
            bf_point.SetMarkerSize(marker_size)
        bf_points.append(copy.deepcopy(bf_point))   # Each iteration of the loop would overwrite the last if we didn't copy
        bf_points[len(bf_points)-1].Draw("same")
        #conf_legend.AddEntry(bf_points[-1], marker_graph.GetTitle(), "p")


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
        if not file or not file.IsOpen():
            print "ERROR: Could not open the file \'{}\'".format(filename)
            sys.exit()
        histo = None
        for hname in potential_chisq_hist_names:
            if file.GetListOfKeys().Contains(hname):
                histo = file.Get(hname)
                histo.SetDirectory(0)
                histo.SetTitle(filename)    # So we know what fit this hist corresponds to when getting the BF point and ranges

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
            if not file_other_MH or not file_other_MH.IsOpen():
                print "ERROR: Could not open the file \'{}\'".format(filename_other_MH)
                sys.exit()
            histo_other_MH = None
            for hname in potential_chisq_hist_names:
                if file_other_MH.GetListOfKeys().Contains(hname):
                    histo_other_MH = file_other_MH.Get(hname)
                    histo_other_MH.SetDirectory(0)
                    histo_other_MH.SetTitle(filename_other_MH)    # So we know what fit this hist corresponds to when getting the BF point and ranges

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
def print_best_fit_and_range_1d_const_chisq(dchisq_objects, conf_levels):
    for hist in dchisq_objects:
        hist_left_edge  = hist.GetXaxis().GetBinLowEdge(1)
        hist_right_edge = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
        vec_x = ROOT.vector('double')()
        vec_x.push_back(hist_left_edge)
        vec_x.push_back(hist_right_edge)
        bin_bf = hist.GetMinimumBin()
        val_bf = hist.GetBinCenter(bin_bf)

        print "======================================="
        print "For the fit from file \'{}\':".format(hist.GetTitle())
        print "Best fit = {:.4f}".format(val_bf)
        print "Confidence regions:"

        for conf_level in conf_levels:
            critval = util.DeltaChisquare(conf_level, 1)    # The 2nd arg is the number of dimensions
            vec_y = ROOT.vector('double')()
            vec_y.push_back(critval)
            vec_y.push_back(critval)
            critval_tgraph = ROOT.TGraph(vec_x.size(), vec_x.data(), vec_y.data())
            [left_boundary, right_boundary] = calc_1d_fc_xrange(hist, critval_tgraph)    # This function was written for FC, but works for const chisq too

            # Do we have a cyclic variable such as deltacp? (where e.g. if we don't find a lower limit by searching left of
            # the best fit point, then we can search left of the hist's right-hand edge.)
            # If so, does the hist span the full range of this parameter?
            cyclic_limits = False
            if "cpv" in hist.GetTitle() or "deltacp" in hist.GetTitle():
                first_bin_centre = hist.GetBinCenter(1)
                last_bin_centre  = hist.GetBinCenter(hist.GetNbinsX())
                if ROOT.TMath.Abs(first_bin_centre + ROOT.TMath.Pi()) < conv.kASmallNum and \
                        ROOT.TMath.Abs(last_bin_centre  - ROOT.TMath.Pi()) < conv.kASmallNum:
                    cyclic_limits = True

            left_boundary_for_err  = left_boundary
            right_boundary_for_err = right_boundary
            if left_boundary > right_boundary_for_err and cyclic_limits:
                left_boundary_for_err = left_boundary - 2.0*ROOT.TMath.Pi()
            if right_boundary < left_boundary_for_err and cyclic_limits:
                right_boundary_for_err = right_boundary + 2.0*ROOT.TMath.Pi()

            left_boundary_str = "{:.4f}".format(left_boundary)
            err_low_str       = "{:.4f}".format(val_bf - left_boundary_for_err)
            if left_boundary == conv.kBadResult:
                left_boundary_str = "not found"
                err_low_str       = "not found"
            right_boundary_str = "{:.4f}".format(right_boundary)
            err_high_str       = "{:.4f}".format(right_boundary_for_err - val_bf)
            if right_boundary == conv.kBadResult:
                right_boundary_str = "not found"
                err_high_str       = "not found"

            print "   ", util.ConfLevelAsString(conf_level), "CR: [", left_boundary_str, ",", right_boundary_str, \
                  "]  -->  {:.4f}".format(val_bf), "+", err_high_str, "-", err_low_str
            del critval_tgraph
            del vec_y

        print "======================================="
        del vec_x

#________________________________________________________________________________________________________
def get_ptheta_objects(ptheta_filename, num_var, ptheta_labels, ptheta_colours, ptheta_styles, chisq_calc_type, valor_filenames, conf_levels):
    objects = []

    ptheta_potential_chisq_hist_names = ["SA_dcp", "SA_dcp_IH"]
    ptheta_potential_contour_graph_names = "cont_like_all_conv[IS_IH]_[CL]_[ICONT]"   # e.g. cont_like_all_conv_68_0, cont_like_all_conv_90_1, or cont_like_all_conv_IH_997_0
    ptheta_bf_pt_name = "BestFit_Point_All_conv[IS_IH]" # e.g. BestFit_Point_All_conv or BestFit_Point_All_conv_IH

    cl_valor_to_ptheta = {
        conv.kCL68     : "68",
        conv.kCL90     : "90",
        conv.kCL95     : "95",
        conv.kCL99     : "99",
        conv.kCL1Sigma : "68",
        conv.kCL2Sigma : "954",
        conv.kCL3Sigma : "997"
    }

    # Are we exclusively plotting file files of a single mass hierarchy?
    # We need to know which MH we're plotting for VALOR because both MHs are included in each p-theta file
    _found_valor_nh = False
    _found_valor_ih = False
    for fn in valor_filenames:
        if "normal" in fn:
            _found_valor_nh = True
        elif "inverted" in fn:
            _found_valor_ih = True
        else:
            print "ERROR: Could not determine the true mass hierarchy from VALOR filename \'{}\'".format(fn)
            sys.exit()

    # p-theta 1D hists have deltachisq = chisq - chisq_min_NH (i.e. same as the MHpair option),
    # so for simplicity, only allow for this option to be used when comparing to 1D p-theta hists.
    if num_var == 1:
        if chisq_calc_type == "MHpair":
            _found_valor_nh = True
            _found_valor_ih = True
        else:
            print "ERROR: Tried to plot a mix of NH and IH hists with a comparison to p-theta, but the chisq calculation type was not set to \'MHpair\'"
            sys.exit()

    file_ptheta = ROOT.TFile(ptheta_filename, "read")
    if not file_ptheta or not file_ptheta.IsOpen():
        print "ERROR: Could not open the p-theta comparison file \'{}\'".format(ptheta_filename)
        sys.exit()

    if num_var == 1:
        # Get 1D deltachisq hists
        for i, hname in enumerate(ptheta_potential_chisq_hist_names):
            if file_ptheta.GetListOfKeys().Contains(hname):
                if (_found_valor_nh and "IH" not in hname) or (_found_valor_ih and "IH" in hname):
                    histo = file_ptheta.Get(hname)
                    histo.SetDirectory(0)
                    if get_ndims(histo) != 1:
                        print "ERROR: Read p-theta hist with", get_ndims(histo), "dimensions. Expected 1."
                        sys.exit()
                    if chisq_calc_type == "MHpair":
                        histo.SetTitle(ptheta_labels[0])     # We'll use this title for the legend later
                    else:
                        histo.SetTitle(ptheta_labels[i])     # We'll use this title for the legend later
                    histo.SetLineColor(ptheta_colours[i])
                    histo.SetLineStyle(ptheta_styles [i])
                    objects.append(histo)
    else:
        # Get contour TGraphs
        iterate_hier_label = []
        if _found_valor_nh:
            iterate_hier_label.append("")
        if _found_valor_ih:
            iterate_hier_label.append("_IH")

        for iconf, conf_level in enumerate(conf_levels):
            ptheta_conf_level_str = cl_valor_to_ptheta[conf_level]
            for ih, h in enumerate(iterate_hier_label):
                conf_regions_at_this_cl = []
                for icont in range(0, 10):
                    ptheta_contour_name = ptheta_potential_contour_graph_names
                    ptheta_contour_name = ptheta_contour_name.replace("[CL]"   , ptheta_conf_level_str)
                    ptheta_contour_name = ptheta_contour_name.replace("[IS_IH]", h)   # There can be multiple confidence regions at a single confidence level
                    ptheta_contour_name = ptheta_contour_name.replace("[ICONT]", str(icont))   # There can be multiple confidence regions at a single confidence level
                    if file_ptheta.GetListOfKeys().Contains(ptheta_contour_name):
                        if (_found_valor_nh and "IH" not in h) or (_found_valor_ih and "IH" in h):
                            if not (_found_valor_nh and _found_valor_ih):
                                ih = 0
                            isty = ih
                            contour = file_ptheta.Get(ptheta_contour_name)
                            MH_str = "NO"
                            if "IH" in h:
                                MH_str = "IO"
                            label = ptheta_labels[isty]
                            contour.SetTitle(label)     # We'll use this title for the legend later
                            contour.SetLineColor(ptheta_colours[isty])
                            contour.SetLineStyle(ptheta_styles [isty])
                            conf_regions_at_this_cl.append(contour)
                objects.append(conf_regions_at_this_cl)

        # Get BF point markers
        for h in iterate_hier_label:
            ptheta_bf_name = ptheta_bf_pt_name
            ptheta_bf_name = ptheta_bf_name.replace("[IS_IH]", h)
            if file_ptheta.GetListOfKeys().Contains(ptheta_bf_name):
                if (_found_valor_nh and "IH" not in h) or (_found_valor_ih and "IH" in h):
                    bf = file_ptheta.Get(ptheta_bf_name)
                    MH_str = " NO"
                    if "IH" in h:
                        MH_str = " IO"
                    if (_found_valor_nh or _found_valor_ih) and not (_found_valor_nh and _found_valor_ih):
                        MH_str = ""
                    bf.SetTitle("p-theta" + MH_str)     # We'll use this title for the legend later
                    objects.append(bf)

    return objects

#________________________________________________________________________________________________________
def get_mach3_objects(mach3_filename, num_var, conf_levels, mach3_labels, mach3_colours, mach3_styles, chisq_calc_type):
    objects = []

    # The contours drawn in MaCh3's graphs are just the same 2D posterior probability hist drawn multiple times with different CLs
    mach3_potential_chisq_hist_names = ["h_dcp", "h_dm32", "h_th13", "h_th23"]
    mach3_potential_contour_graph_names = "g[CL]_conf_[ICONT]"   # e.g. g68_conf_0, g90_conf_0, or g_90_conf_1 etc.
    mach3_bf_pt_name = "Bestfit"

    # TODO: check these, not sure about the 2 and 3 sigma strings
    cl_valor_to_mach3 = {
        conv.kCL68     : "68",
        conv.kCL90     : "90",
        conv.kCL95     : "95",
        conv.kCL99     : "99",
        conv.kCL1Sigma : "68",
        conv.kCL2Sigma : "954",
        conv.kCL3Sigma : "997"
    }


    file_mach3 = ROOT.TFile(mach3_filename, "read")
    if not file_mach3 or not file_mach3.IsOpen():
        print "ERROR: Could not get MaCh3 input file \'{}\'".format(mach3_filename)
        sys.exit()

    for i, hname in enumerate(mach3_potential_chisq_hist_names):
        if chisq_calc_type == "MHpair":
            # Get the posterior prob hist for the other hierarchy
            filename_other_MH = ""
            if "nh" in mach3_filename or "NH" in mach3_filename:
                filename_other_MH = mach3_filename   .replace("nh", "ih")
                filename_other_MH = filename_other_MH.replace("NH", "IH")
            else:
                filename_other_MH = mach3_filename   .replace("ih", "nh")
                filename_other_MH = filename_other_MH.replace("IH", "NH")
            file_other_MH = ROOT.TFile(filename_other_MH, "READ")
            if not file_other_MH or not file_other_MH.IsOpen():
                print "ERROR: Could not get MaCh3 input file \'{}\'".format(filename_other_MH)
                sys.exit()

            filename_both_MH = ""
            if "nh" in mach3_filename or "NH" in mach3_filename:
                filename_both_MH = mach3_filename  .replace("nh", "both")
                filename_both_MH = filename_both_MH.replace("NH", "both")
            else:
                filename_both_MH = mach3_filename   .replace("ih", "both")
                filename_both_MH = filename_both_MH.replace("IH", "both")
            file_both_MH = ROOT.TFile(filename_both_MH, "READ")
            if not file_both_MH or not file_both_MH.IsOpen():
                print "ERROR: Could not get MaCh3 input file \'{}\'".format(filename_both_MH)
                sys.exit()

            if num_var == 1:
                hist_posterior = get_single_mach3_hist(file_mach3, hname, num_var)
                if not hist_posterior:
                    continue

                hist_posterior_other_MH = get_single_mach3_hist(file_other_MH, hname, num_var)
                if not hist_posterior_other_MH:
                    print "ERROR: Could not get MaCh3 hist named \'{}\'".format(hname),  "from file \'{}\'".format(filename_other_MH)
                    sys.exit()

                if "dm23" in mach3_filename or "th23" in mach3_filename:
                    # The single MH posterior prob hists have been normalised by the posterior prob of each MO, so we can't
                    # directly extract the NO-IO chisq separation from them. To get this, we need to calculate the
                    # posterior prob of each MO from the given information (and that PNO + PIO = 1).
                    # This can be calculated using a single bin, but the bin content must be non-zero
                    # Try a bin in the middle of the posterior prob hist, as the middle is the least likely region to be empty
                    hist_posterior_both_MH = get_single_mach3_hist(file_both_MH, hname, num_var)
                    if not hist_posterior_both_MH:
                        print "ERROR: Could not get MaCh3 hist named \'{}\'".format(hname),  "from file \'{}\'".format(filename_both_MH)
                        sys.exit()

                    test_bin = hist_posterior.GetMaximumBin()
                    test_bin_IO = test_bin
                    if "dm23" in mach3_filename:
                        dm12sq = 7.53e-5
                        IO_cent = hist_posterior.GetBinCenter(test_bin) + dm12sq
                        test_bin_IO = hist_posterior_other_MH.FindBin(IO_cent)
                    assert(hist_posterior         .GetBinContent(test_bin) > 0.0)
                    assert(hist_posterior_other_MH.GetBinContent(test_bin_IO) > 0.0)
                    assert(hist_posterior_both_MH .GetBinContent(test_bin) > 0.0)

                    # TODO: I don't think this will work if hist_posterior corresponds to IO - check this
                    assert("NH" in mach3_filename or "nh" in mach3_filename)
                    PNO = (hist_posterior_both_MH.GetBinContent(test_bin) - hist_posterior_other_MH.GetBinContent(test_bin_IO)) / (hist_posterior.GetBinContent(test_bin) - hist_posterior_other_MH.GetBinContent(test_bin_IO))
                    hist_posterior         .Scale(PNO)
                    hist_posterior_other_MH.Scale(1.0 - PNO)

                hist_chisq = convert_posterior_prob_to_chisq(hist_posterior, num_var, file_mach3.GetName())
                hist_chisq.SetLineColor(mach3_colours[0])
                hist_chisq.SetLineStyle(mach3_styles [0])
                hist_chisq.SetFillColorAlpha(0, 0.0)
                hist_chisq.SetTitle(mach3_labels[0])
                min_chisq1 = hist_chisq.GetMinimum()

                hist_chisq_other_MH = convert_posterior_prob_to_chisq(hist_posterior_other_MH, num_var, filename_other_MH)
                hist_chisq_other_MH.SetLineColor(mach3_colours[1])
                hist_chisq_other_MH.SetLineStyle(mach3_styles [1])
                hist_chisq_other_MH.SetFillColorAlpha(0, 0.0)
                hist_chisq_other_MH.SetTitle(mach3_labels[0])
                min_chisq2 = hist_chisq_other_MH.GetMinimum()

                if min_chisq1 <= min_chisq2:
                    convert_to_delta_chisq([hist_chisq_other_MH], num_var, min_chisq1)  # negative chisq_min because we want to add it to the IH
                    convert_to_delta_chisq([hist_chisq], num_var)
                else:
                    convert_to_delta_chisq([hist_chisq], num_var, min_chisq2)
                    convert_to_delta_chisq([hist_chisq_other_MH], num_var)

                objects.append(hist_chisq)
                objects.append(hist_chisq_other_MH)
            else:
                print "ERROR: Plotting 2D MaCh3 contours with the chisq calculation set to \'MHpair\' is currently unsupported"
                sys.exit()
        else:
            if num_var == 1:
                hist_posterior = get_single_mach3_hist(file_mach3, hname, num_var)
                if not hist_posterior:
                    continue

                hist_chisq = convert_posterior_prob_to_chisq(hist_posterior, num_var, file_mach3.GetName())
                hist_chisq.SetLineColor(mach3_colours[0])
                hist_chisq.SetLineStyle(mach3_styles [0])
                hist_chisq.SetFillColorAlpha(0, 0.0)
                hist_chisq.SetTitle(mach3_labels[0])
                convert_to_delta_chisq([hist_chisq], num_var)
                objects.append(hist_chisq)
            else:
                # Get contour TGraphs
                for iconf, conf_level in enumerate(conf_levels):
                    mach3_conf_level_str = cl_valor_to_mach3[conf_level]
                    conf_regions_at_this_cl = []
                    for icont in range(0, 10):
                        mach3_contour_name = mach3_potential_contour_graph_names
                        mach3_contour_name = mach3_contour_name.replace("[CL]"   , mach3_conf_level_str)
                        mach3_contour_name = mach3_contour_name.replace("[ICONT]", str(icont))   # There can be multiple confidence regions at a single confidence level
                        if file_mach3.GetListOfKeys().Contains(mach3_contour_name):
                            contour = file_mach3.Get(mach3_contour_name)

                            if "dm23" in file_mach3.GetName() and ("IH" in file_mach3.GetName() or "ih" in file_mach3.GetName()):
                                # Then this is really Dm31sq, which is negative. We need this to be directly comparable to Dm23sq,
                                # so we want to convert to |Dm31sq|
                                contour = convert_mach3_Dm31sq_to_absDm31sq(contour)

                            label = mach3_labels[0]
                            contour.SetTitle(label)     # We'll use this title for the legend later
                            contour.SetLineColor(mach3_colours[iconf])
                            contour.SetLineStyle(mach3_styles [iconf])
                            conf_regions_at_this_cl.append(contour)
                    if len(conf_regions_at_this_cl) == 0:
                        print "ERROR: Could not find a MaCh3 contour at the", util.ConfLevelAsString(conf_level), "CL"
                        sys.exit()
                    objects.append(conf_regions_at_this_cl)

                # Get BF point markers
                if file_mach3.GetListOfKeys().Contains(mach3_bf_pt_name):
                    bf = file_mach3.Get(mach3_bf_pt_name)
                    bf.SetTitle("MaCh3")     # We'll use this title for the legend later
                    objects.append(bf)
        break

    # Convert the chisq hists into Delta chisq hists
    if num_var == 1:
        if chisq_calc_type == "global":
            convert_to_delta_chisq(objects, num_var)
        elif chisq_calc_type == "local":
            chisq_min_all = get_chisq_min(dchisq_objects)
            convert_to_delta_chisq(objects, num_var, chisq_min_all)

    return objects

#________________________________________________________________________________________________________
def get_single_mach3_hist(file, hist_name, ndim):
    mach3_potential_canvas_names = ["c1", "c"]
    canv = None
    for canv_name in mach3_potential_canvas_names:
        if file.GetListOfKeys().Contains(canv_name):
            canv = file.Get(canv_name)
            break
    if not canv:
        print "ERROR: Could not get the MaCh3 canvas from file \'{}\'".format(file.GetName())
        sys.exit()

    # The MaCh3 ROOT files store the canvas of some plot, rather than storing the drawn objects individually,
    # so we must extract the objects from the canvas
    plot_objects = canv.GetListOfPrimitives()
    if not plot_objects:
        print "ERROR: Could not get the list of drawn objects from the MaCh3 canvas in file \'{}\'".format(file.GetName())
        sys.exit()

    if not plot_objects.FindObject(hist_name):
        return None

    # MaCh3 is Bayesian, so uses posterior probability instead of chisq.
    hist_posterior = plot_objects.FindObject(hist_name)
    if not hist_posterior:
        print "ERROR: Could not get the hist named file \'{}\'".format(file.GetName())
        sys.exit()
    hist_posterior.SetDirectory(0)
    if get_ndims(hist_posterior) != ndim:
        print "ERROR: Read MaCh3 hist with", get_ndims(hist_posterior), "dimensions. Expected", ndim
        sys.exit()

    if "dm23" in file.GetName() and ("IH" in file.GetName() or "ih" in file.GetName()):
        # Then this is really Dm31sq, which is negative. We need this to be directly comparable to Dm23sq,
        # so we want to convert to |Dm31sq|
        hist_posterior = convert_mach3_Dm31sq_to_absDm31sq(hist_posterior)

    return hist_posterior

#________________________________________________________________________________________________________
def convert_mach3_Dm31sq_to_absDm31sq(obj_posterior):
    dm12sq = 7.53e-5    # PDG 2019 central value, I think
    if type(obj_posterior) == ROOT.TH1D:
        nbins = obj_posterior.GetNbinsX()
        low = obj_posterior.GetBinLowEdge(1)
        high = obj_posterior.GetBinLowEdge(nbins) + obj_posterior.GetBinWidth(nbins)
        hist_posterior_new = ROOT.TH1D(obj_posterior.GetName()+"_abs", obj_posterior.GetTitle()+"_abs", nbins, low, high)
        for ib in range(1, nbins+1):
            fill = obj_posterior.GetBinContent(ib)
            cent = obj_posterior.GetBinCenter (ib)
            if cent < 0.0:
                hist_posterior_new.Fill(ROOT.TMath.Abs(cent)-dm12sq, fill)
        obj_posterior = copy.deepcopy(hist_posterior_new)
        del hist_posterior_new
    elif type(obj_posterior) == ROOT.TGraph:
        cont_new = ROOT.TGraph()
        for ip in range(0, obj_posterior.GetN()):
            x = ROOT.Double(-999)
            y = ROOT.Double(-999)
            obj_posterior.GetPoint(ip, x, y)
            y_new = ROOT.TMath.Abs(y) - dm12sq
            cont_new.SetPoint(cont_new.GetN(), x, y_new)
        obj_posterior = copy.deepcopy(cont_new)
        del cont_new
    else:
        print "ERROR: unrecognised MaCh3 object type:", type(obj_posterior)
        sys.exit()

    return obj_posterior
    

#________________________________________________________________________________________________________
def convert_posterior_prob_to_chisq(hist_posterior, ndim, filename):
    apply_penalty_to_mach3_hist = False # This shouldn't be needed, as I'm told that the MaCh3 posterior prob should be directly equivalent to our likelihood... somehow...
    hist_chisq = None

    # Convert from posterior probability to Delta chisq, with min = min of this chisq hist only
    # (not chisq as the conversion from posterior probability to chisq results in negative chisq values)
    if ndim == 1:
        hist_chisq = ROOT.TH1D(hist_posterior)
    else:
        hist_chisq = ROOT.TH2D(hist_posterior)
    hist_chisq.Reset()

    # Fill chisq hist
    for ib in range(1, hist_posterior.GetNbinsX()+1):
        if ndim == 2:
            for jb in range(1, hist_posterior.GetNbinsY() + 1):
                if (hist_posterior.GetBinContent(ib, jb) > 0):
                    fill = -2.0*ROOT.TMath.Log(hist_posterior.GetBinContent(ib, jb))
                    if apply_penalty_to_mach3_hist:
                        fill += get_chisq_penalty(filename, ndim, hist_posterior, ib, jb)
                else: 
                    fill = 999.0 # as MaCh3 MaCh3 has a finite number of steps in parameter space, so where the posterior prob is really small, 
                                 #they'll get no steps, so the input MaCh3 hists will have bins with posterior prob = 0.  ln(0) = inf.
                                 #This give us weird plots. Here is a workaround to fix it. 
                hist_chisq.SetBinContent(ib, jb, fill)
        else:
            if (hist_posterior.GetBinContent(ib) >0 ):
                fill = -2.0*ROOT.TMath.Log(hist_posterior.GetBinContent(ib))
                if apply_penalty_to_mach3_hist:
                    fill += get_chisq_penalty(filename, ndim, hist_posterior, ib, None)
            else:
                fill = 990.0 #Same trick as for 2D contours 
            hist_chisq.SetBinContent(ib, fill) 
            

    hist_chisq.SetDirectory(0)

    return hist_chisq

#________________________________________________________________________________________________________
def get_chisq_penalty(filename, ndim, hist_posterior, xbin, ybin):
    chisq_penalty = 0.0

    # The mixing13 prior is not flat, so we must account for it in the posterior prob to
    # likelihood calculation (because of Bayes' theorem).
    if "th13" in filename:
        prior_ss2th13_mean  = 0.0853
        prior_ss2th13_sigma = 0.0027

        ssth13 = -999
        if "13" in hist_posterior.GetXaxis().GetTitle():
            ssth13 = hist_posterior.GetXaxis().GetBinCenter(xbin)
        elif ndim == 2:
            ssth13 = hist_posterior.GetYaxis().GetBinCenter(ybin)
        else:
            return 0.0
        ss2th13 = util.SingleAngleToDoubleAngle(ssth13)

        chisq_penalty += ROOT.TMath.Power(ss2th13 - prior_ss2th13_mean,2) / ROOT.TMath.Power(prior_ss2th13_sigma,2)

    return chisq_penalty

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
# NOTE: Although this function was written with FC in mind, it can still be used for
#       calculating const dchisq ranges, as long as a flat, two-point TGraph is input.
def calc_1d_fc_xrange(histo, fc_critval_tgraph):
    lim_x_min = None
    lim_x_max = None

    # Do we have a cyclic variable such as deltacp? (where e.g. if we don't find a lower limit by searching left of
    # the best fit point, then we can search left of the hist's right-hand edge.)
    # If so, does the hist span the full range of this parameter?
    cyclic_limits = False
    if "cpv" in histo.GetTitle() or "deltacp" in histo.GetTitle():
        first_bin_centre = histo.GetBinCenter(1)
        last_bin_centre  = histo.GetBinCenter(histo.GetNbinsX())
        if ROOT.TMath.Abs(first_bin_centre + ROOT.TMath.Pi()) < conv.kASmallNum and \
                ROOT.TMath.Abs(last_bin_centre  - ROOT.TMath.Pi()) < conv.kASmallNum:
            cyclic_limits = True

    # Shift the input Delta chisq hist by the critical values in each bin centre.
    # This allows us to calculate the range by comparing to a constant (shifted) critical value of 0.
    # This is equivalent to comparing to a variable critical value, but is easier to handle
    critval_shifted = 0.0
    hist_shifted = ROOT.TH1D(histo)
    min_shifted =  999.0
    max_shifted = -999.0
    for ibin in range(1, hist_shifted.GetNbinsX()+1):
        bincenter = histo.GetBinCenter(ibin)
        bincont = histo.GetBinContent(ibin)
        bincont_new = bincont - fc_critval_tgraph.Eval(bincenter)
        hist_shifted.SetBinContent(ibin, bincont_new)
        min_shifted = min(min_shifted, bincont_new)
        max_shifted = max(max_shifted, bincont_new)

    # Does the allowed region exist at all for this hist? If not, return some nonsense numbers
    if max_shifted < critval_shifted or min_shifted > critval_shifted:
        return [conv.kBadResult, conv.kBadResult]

    bin_bf = hist_shifted.GetMinimumBin()

    # Bins before (b) and after (a) the intersection with critical value
    bin_left_search_a  = conv.kBadResult
    bin_left_search_b  = conv.kBadResult
    bin_right_search_a = conv.kBadResult
    bin_right_search_b = conv.kBadResult

    # Search left of the chisq min for the boundary of the allowed FC region,
    # wrapping around, if needed and if possible
    found = False
    for ibin in range(bin_bf, 0, -1):
        bincont = hist_shifted.GetBinContent(ibin)
        if bincont > critval_shifted:
            # bins before (b) and after (a) DeltaChisq = critical value (following x direction)
            bin_left_search_b = ibin
            bin_left_search_a = bin_left_search_b+1
            found = True
            break
    if cyclic_limits and not found:
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
        # If no bin is found it means the DeltaChisq = critical value bin is outside the scanned region.
        # If so, set the left extreme to 1 (first bin of the histogram)
        print "WARNING: Could not form an allowed FC range when looking left of the chisq min"
        bin_left_search_a = 1
        bin_left_search_b = 1

    # Search right of the chisq min for the boundary of the allowed FC region,
    # wrapping around, if needed and if possible
    found = False
    for ibin in range(bin_bf, hist_shifted.GetNbinsX()+1):
        bincont = hist_shifted.GetBinContent(ibin)
        if bincont > critval_shifted:
            # bins before (b) and after (a) DeltaChisq = critical value
            bin_right_search_a = ibin
            bin_right_search_b = bin_right_search_a-1
            found = True
            break
    if cyclic_limits and not found:
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
        if cyclic_limits:
            lim_x_min = hist_shifted.GetBinCenter(bin_left_search_a)
        else:
            lim_x_min = conv.kBadResult

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
        if cyclic_limits:
            lim_x_max = hist_shifted.GetBinCenter(bin_right_search_b)
        else:
            lim_x_max = conv.kBadResult

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
    parser.add_argument("-f", "--filename"        , dest="in_filenames"        , nargs="+"              , help="Input filename(s) for contour(s).")
    parser.add_argument("-l", "--label"           , dest="labels"              , nargs="+"              , help="Legend label(s) for contours(s). This should have either 0 length or the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument("-c", "--colour"          , dest="colours"             , nargs="+"              , help="Custom line colour(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the colours is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument("-s", "--style"           , dest="styles"              , nargs="+"              , help="Custom line style(s) for contours(s) (optional). If used, this should have the same number of arguments as provided to --filename  (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the styles is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--best-fit-label"  , dest="best_fit_labels"     , nargs="+"              , help="Labels for best-fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--best-fit-style"  , dest="best_fit_styles"     , nargs="+"              , help="Marker styles for best fit point(s) for 2D contours(s) (optional). If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--conf-levels"     , dest="conf_levels_str"     , nargs="+"              , help="The confidence levels to plot const chisq lines/contours for. See ConfLevel::AsString() for the choices")
    parser.add_argument(      "--title-x"         , dest="title_x"                                      , help="Override default x-axis title")
    parser.add_argument(      "--title-y"         , dest="title_y"             , default="-2#Deltaln(L)", help="Override default y-axis title")
    parser.add_argument(      "--line-width"      , dest="lwidth"              , type=int               , help="Contour line width")
    parser.add_argument(      "--marker-size"     , dest="msize"               , type=float             , help="Best-fit marker size")
    parser.add_argument(      "--text-size"       , dest="tsize"               , type=float             , help="Legend text size")
    parser.add_argument(      "--title-size"      , dest="titlesize"           , type=float,default=0.05, help="Title text size")
    parser.add_argument("-o", "--output"          , dest="out_filename"                                 , help="Output PDF filename")
    parser.add_argument(      "--chisq-calc"      , dest="chisq_calc"                                   , help="If \'global\', the min. Dchisq of each contour is set to 0. If \'local\', Dchisq is calculated relative to the min. of all contours. If \'MHpair\', the contour for the opposite MH is found and the Dchisq of each contour set is relative to the min. of pair of NH/IH contours")
    parser.add_argument(      "--fc-files"        , dest="fc_files"            , nargs="+"              , help="ROOT files containing TGraphs of the FC critical values as a function of the parameter of interest and the mass hierarchy. For 1D fits only. If used, this should have the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used, but leaving some args as \'None\' is fine if not every dchisq surface/contour file has an associated FC file.")
    parser.add_argument(      "--fc-graph-pattern", dest="fc_graph_pattern"    , nargs="+"              , help="Naming pattern of the TGraphs containing the FC critical values with the confidence level of each TGraph replaced with \'[CL]\', e.g. \'gFC_1sigma\'. If used, this should have the same number of arguments as provided to --filename.")
    parser.add_argument(      "--fc-line-styles"  , dest="fc_line_styles"      , nargs="+"              , help="Line styles used for plotting the FC critical values with \'--fc-draw-vals\'. If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels), or double if \'--chisq-calc MHpair\' is used, the order of the styles is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--fc-line-colours" , dest="fc_line_colours"     , nargs="+"              , help="Line colours used for plotting the FC critical values with \'--fc-draw-vals\'. If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels), or double if \'--chisq-calc MHpair\' is used, the order of the colours is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--fc-draw-vals"    , dest="fc_draw_vals"        , action="store_true"    , help="If true, draws the FC critical values on top of the Delta chisq hists")
    parser.add_argument(      "--fc-draw-boxes"   , dest="fc_draw_boxes"       , action="store_true"    , help="If true, draws shaded boxes representing the FC allowed region on top of the Delta chisq hists")
    parser.add_argument(      "--t2k"             , dest="is_t2k_only"         , action="store_true"    )
    parser.add_argument(      "--prelim"          , dest="prelim"              , action="store_true"    )
    parser.add_argument(      "--theta23"         , dest="theta23"             , type=float             , help="Input value of theta23.")
    parser.add_argument(      "--dm23"            , dest="dm23"                , type=float             , help="Input value of dm23.")
    parser.add_argument(      "--draw-reactor"    , dest="draw_reactor"        , action="store_true"    , help="Draw a vertical band indicating the PDG 2018 reactor constraint on sinsq(theta13) used in the T2K oscillation analysis")
    parser.add_argument(      "--x-range"         , dest="x_range_str"         , nargs="+"              , help="x-axis range in user coordinates. Specify as min, max")
    parser.add_argument(      "--y-range"         , dest="y_range_str"         , nargs="+"              , help="y-axis range in user coordinates. Specify as min, max")
    parser.add_argument(      "--max-ndigits-y"   , dest="max_ndigits_y"       , type=int, default=3    , help="Max number of digits on the y-axis labels.")
    parser.add_argument(      "--leg-coords"      , dest="leg_coords_str"      , nargs="+"              , help="Legend corner coordinates [xmin, ymin, xmax, ymax]. This legend displays the styles of each fit, e.g. NH/IH or VALOR/MaCh3/p-theta.")
    parser.add_argument(      "--leg2-coords"     , dest="leg2_coords_str"     , nargs="+"              , help="2nd Legend corner coordinates [xmin, ymin, xmax, ymax]. This legend displays the styles of each confidence level and of the best fit point markers (2D plots only), or the style of each MH (1D plots only)")

    # For comparisons to other groups
    parser.add_argument(      "--compare-ptheta"  , dest="compare_ptheta"                            , help="Input filename for plotting comparisons to p-theta. Only file allowed (as multiple hierarchies in the same file).")
    parser.add_argument(      "--ptheta-label"    , dest="ptheta_labels"    , nargs="+"              , help="Legend label(s) for p-theta contours(s). This should have either 0 length or the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--ptheta-colour"   , dest="ptheta_colours"   , nargs="+"              , help="Custom line colour(s) for p-theta contours(s) (optional). If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the colours is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--ptheta-style"    , dest="ptheta_styles"    , nargs="+"              , help="Custom line style(s) for p-theta contours(s) (optional). If used, this should have the same number of arguments as provided to --filename  (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the styles is taken to be: file0+CL0, file0+CL1, file1+CL0, file1+CL1, etc.")
    parser.add_argument(      "--compare-mach3"   , dest="compare_mach3"                             , help="Input filename for plotting comparisons to MaCh3.")
    parser.add_argument(      "--mach3-label"     , dest="mach3_labels"     , nargs="+"              , help="Legend label(s) for MaCh3 contours(s). This should have either 0 length or the same number of arguments as provided to --filename, or double if \'--chisq-calc MHpair\' is used.")
    parser.add_argument(      "--mach3-colour"    , dest="mach3_colours"    , nargs="+"              , help="Custom line colour(s) for MaCh3 contours(s) (optional). If used, this should have the same number of arguments as provided to --filename (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the colours is taken to be: NH+CL0, NH+CL1, IH+CL0, IH+CL1, etc. (one MH used if VALOR is only using one)")
    parser.add_argument(      "--mach3-style"     , dest="mach3_styles"     , nargs="+"              , help="Custom line style(s) for MaCh3 contours(s) (optional). If used, this should have the same number of arguments as provided to --filename  (times the number of confidence levels for 2D plots), or double if \'--chisq-calc MHpair\' is used. e.g. for 2D plots, the order of the styles is taken to be:  NH+CL0, NH+CL1, IH+CL0, IH+CL1, etc. (one MH used if VALOR is only using one)")
    args = parser.parse_args()

    # Check the validity of the input arguments

    dchisq_objects = []
    num_var = get_dchisq_histos(args.in_filenames, args.chisq_calc, dchisq_objects)
    print "Number of chisq hists to plot =", len(dchisq_objects)

    labels = []
    if args.labels:
        labels = [i for i in args.labels]
        if num_var == 1 and len(args.leg2_coords_str) > 0:
            if not len(args.in_filenames) == len(labels):
                print "Error: The number of filenames and labels should match"
                sys.exit()
        else:
            if not len(args.in_filenames) == len(labels) and args.chisq_calc != "MHpair":
                print "Error: The number of filenames and labels should match"
                sys.exit()
            elif not 2 * len(args.in_filenames) == len(labels) and args.chisq_calc == "MHpair":
                print "Error: The number of labels should be twice the number of filenames"
                print "filenames =", args.in_filenames
                print "labels =", labels
                sys.exit()

    conf_levels = []
    if not args.conf_levels_str:
        conf_levels = [conv.kCL1Sigma, conv.kCL90]
        print "Setting default confidence levels to [1sigma, 90%]"
    else:
        conf_levels = [util.ConfLevelFromString(conf_str) for conf_str in args.conf_levels_str]

    ptheta_objects = None
    if args.compare_ptheta:
        ptheta_colours = [int(c) for c in args.ptheta_colours]
        ptheta_styles  = [int(c) for c in args.ptheta_styles]
        ptheta_objects = get_ptheta_objects(args.compare_ptheta, num_var, args.ptheta_labels, ptheta_colours,
                                            ptheta_styles, args.chisq_calc, args.in_filenames, conf_levels)

    mach3_objects = None
    if args.compare_mach3:
        mach3_colours = [int(c) for c in args.mach3_colours]
        mach3_styles  = [int(c) for c in args.mach3_styles]
        mach3_objects  = get_mach3_objects(args.compare_mach3, num_var, conf_levels, args.mach3_labels, mach3_colours,
                                           mach3_styles, args.chisq_calc)

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
        if num_var == 2:
            if args.chisq_calc != "MHpair":
                if len(best_fit_styles) != len(args.in_filenames):
                    print "ERROR: The number of best-fit marker styles must equal the number of input files"
                    print "filenames =", args.in_filenames
                    print "best-fit styles =", args.best_fit_styles
                    sys.exit()
            else:
                if len(best_fit_styles) != 2 * len(args.in_filenames):
                    print "ERROR: The number of best-fit marker styles must be double the number of input files"
                    print "filenames =", args.in_filenames
                    print "best-fit styles =", args.best_fit_styles
                    sys.exit()

    leg_coords = []
    if args.leg_coords_str:
        leg_coords = [float(i) for i in args.leg_coords_str]
        if len(leg_coords) != 4:
            print "ERROR: The number of legend coordinates must equal 4"

    leg2_coords = []
    if args.leg2_coords_str:
        leg2_coords = [float(i) for i in args.leg2_coords_str]
        if len(leg2_coords) != 4:
            print "ERROR: The number of coordinates for the 2nd legend must equal 4. Instead found:", leg2_coords

    x_range = []
    if args.x_range_str:
        x_range = [float(i) for i in args.x_range_str]
        assert(len(x_range) == 2)

    y_range = []
    if args.y_range_str:
        y_range = [float(i) for i in args.y_range_str]
        assert(len(y_range) == 2)

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
        if len(args.fc_graph_pattern) != len(args.fc_files):
            print "ERROR: The number of input FC TGraph names must match the number of input FC files."
            sys.exit()
        if len(args.fc_line_styles) != len(args.fc_files) * len(conf_levels):
            print "ERROR: The number of input FC line styles must match the number of input FC files * number of confidence levels."
            sys.exit()
        if len(args.fc_line_colours) != len(args.fc_files) * len(conf_levels):
            print "ERROR: The number of input FC line colours must match the number of input FC files * number of confidence levels."
            sys.exit()
        if num_var != 1:
            print "ERROR: Only 1D Feldman-Cousins is supported."
            sys.exit()
        if "cp" not in args.title_x.lower():
            print "ERROR: Currently plotting Feldman-Cousins is only supported for deltacp"
            sys.exit()
        if args.fc_draw_boxes and len(conf_levels) != 1:
            print "ERROR: Drawing FC hatch boxes is not supported for != 1 confidence levels."
            sys.exit()
        fc_line_styles  = [int(i) for i in args.fc_line_styles ]
        fc_line_colours = [int(i) for i in args.fc_line_colours]

    ROOT.TGaxis.SetMaxDigits(args.max_ndigits_y)

    if num_var == 1:
        plot_contours_1d(dchisq_objects, labels, args.out_filename, args.is_t2k_only, colours, styles, args.prelim,
                         args.title_x, args.title_y, args.lwidth, args.tsize, conf_levels, leg_coords, args.titlesize,
                         args.fc_files, args.fc_graph_pattern, fc_line_styles, fc_line_colours,
                         args.fc_draw_vals, args.fc_draw_boxes, ptheta_objects, mach3_objects,
                         x_range, y_range, args.chisq_calc, leg2_coords)
    elif num_var == 2:
        plot_contours_2d(dchisq_objects, labels, args.theta23, args.dm23, args.out_filename, colours, styles,
                         best_fit_styles, args.best_fit_labels, args.prelim, args.title_x, args.title_y, args.lwidth,
                         args.msize, args.tsize, conf_levels, leg_coords, leg2_coords, args.titlesize, args.draw_reactor,
                         ptheta_objects, mach3_objects, x_range, y_range)
    else:
        print "Error: Not all contours have the same dimensionality"
        sys.exit()