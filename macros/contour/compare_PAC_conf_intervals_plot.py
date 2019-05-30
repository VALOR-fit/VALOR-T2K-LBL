#!/usr/bin/env python

import collections
import sys
import numpy as np
from array import array
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

#_______________________________________________________________________________________________________________________
# User defined variables
_data_type  = collections.OrderedDict()
_allocation = collections.OrderedDict()
_hierarchy  = collections.OrderedDict()
_POT        = collections.OrderedDict()

# key = name in file, value = name in Tex file
_data_type["pure_AsimovA"] = "Asimov"
_data_type["hybrid"]       = "hybrid"

_allocation["FHC"]     = "all FHC"
_allocation["FHC_RHC"] = "equal"

_hierarchy["normal"]   = "NH"
_hierarchy["inverted"] = "IH"

_POT["7.8E21"] = "7.8"
_POT["10E21"]  = "10"
_POT["15E21"]  = "15"

_box_half_width   = 0.2
_title_width      = 1.5
_title_offset     = 0.015
_half_title_hight = 0.25


#_______________________________________________________________________________________________________________________
def main(file_pattern, outfilename, nbatch, critval_str, var, do_global):

    crit_dchisq = get_crit_dchisq(critval_str)
    print "Input critical value string of", critval_str, "gives critical Delta chisq =", crit_dchisq

    max_lines = len(_hierarchy) * len(_data_type) * len(_allocation) * len(_POT)

    xaxis_range_low  = 0
    xaxis_range_high = 0
    x_title = "Range"
    if var == "dmsq23":
        x_title = "#Delta m^2"
        xaxis_range_low  = 2.3e-3
        xaxis_range_high = 2.7e-3
    elif var == "mixing23":
        x_title = "#sin^{2}(#theta_{23})"
        xaxis_range_low  = 0.3
        xaxis_range_high = 0.7
    elif var == "deltacp":
        x_title = "#delta_{CP} (radians)"
        xaxis_range_low  = -conv.kpi
        xaxis_range_high = conv.kpi

    mg = ROOT.TMultiGraph()
    box_list = []
    title_list = []

    # For deltacp only
    mg_dist_from_pi = ROOT.TMultiGraph()
    mg_dist_from_0  = ROOT.TMultiGraph()
    max_dist_to_pi  = 0
    max_dist_to_0   = 0
    mg_ddchisq_pi = ROOT.TMultiGraph()
    mg_ddchisq_0  = ROOT.TMultiGraph()
    max_ddchisq_pi  = 0
    max_ddchisq_0   = 0

    # force the whole range to be drawn by adding two almost invisible points at each end of the axis
    x_arr1, y_arr1, xerr_arr1, yerr_arr1 = array('d'), array('d'), array('d'), array('d')
    x_arr1.append(xaxis_range_low)
    x_arr1.append(xaxis_range_high)
    y_arr1.append(1)
    y_arr1.append(max_lines)
    xerr_arr1.append(0.0)
    xerr_arr1.append(0.0)
    yerr_arr1.append(0.0)
    yerr_arr1.append(0.0)
    tiny_points = ROOT.TGraphErrors(2, x_arr1, y_arr1, xerr_arr1, yerr_arr1)
    tiny_points.SetMarkerStyle(1)  # Sets the marker to a single do so it can't be seen (not drawing the marker gets rid of the axis for some reason)
    mg.Add(tiny_points,"P")

    irow = max_lines+1
    for hier in _hierarchy:
        hier_tex = _hierarchy[hier]
        
        for type in _data_type:
            type_tex = _data_type[type]
            
            for alloc in _allocation:
                alloc_tex = _allocation[alloc]
                
                for pot in _POT:
                    irow -= 1
                    pot_tex = _POT[pot]

                    #if type == "pure_AsimovA" and alloc == "FHC" and hier == "normal" and pot == "10E21" and var == "deltacp":
                        #print "\n\n\n *** WARNING: Skipping scenario pure_AsimovA, FHC, normal, 10E21, deltacp - it needs fixing! *** \n\n\n"
                        #continue

                    bf_point   = []
                    range_low  = []
                    range_high = []
                    wrap_range = False
                    wrap_range_batch = [False for i in range(0, nbatch)]
                    dchisq_at_pi = []   # for deltacp only
                    dchisq_at_0  = []   # for deltacp only
                    for ibatch in range(0, nbatch):
                        file_in = file_pattern.replace("[BATCH]"     , str(ibatch))
                        file_in = file_in.replace     ("[DATATYPE]"  , type)
                        file_in = file_in.replace     ("[ALLOCATION]", alloc)
                        file_in = file_in.replace     ("[MH]"        , hier)
                        file_in = file_in.replace     ("[POT]"       , pot)
                        file_in = file_in.replace     ("[VAR]"       , var)
                        print "Reading from file", file_in
                        f = ROOT.TFile(file_in, "read")
                        hchisq = None

                        if do_global and hier == "inverted":
                            hchisq = f.Get("h_chisq")
                            # Open up the NH chisq contour to get the minimum
                            filename_NH = file_in.replace("inverted", "normal")
                            fNH = ROOT.TFile(filename_NH, "read")
                            hchisq_NH = fNH.Get("h_chisq")
                            bf_chisq_NH = hchisq_NH.GetMinimum()
                            bf_chisq_IH = hchisq.GetMinimum()
                            bf_chisq_glob = min(bf_chisq_NH, bf_chisq_IH)
                            for i in range(1, hchisq.GetNbinsX() + 1):
                                hchisq.SetBinContent(i, hchisq.GetBinContent(i) - bf_chisq_glob)
                        else:
                            hchisq = f.Get("h_deltachisq")

                        # Get parameter best fit value
                        bf_bin = hchisq.GetMinimumBin()
                        bf_bin_x = ROOT.Long()
                        dummy = ROOT.Long()
                        hchisq.GetBinXYZ(bf_bin, bf_bin_x, dummy, dummy)
                        bf_point.append(hchisq.GetXaxis().GetBinCenter(bf_bin_x))

                        # For deltacp get the chisq values at +/- pi and 0
                        if var == "deltacp":
                            dchisq_at_pi.append(hchisq.GetBinContent(1))    # pi is at the centre of this bin, so no need to extrapolate
                            deltacp_0_bin_glob = hchisq.FindBin(0.0)
                            deltacp_0_bin_x = ROOT.Long()
                            hchisq.GetBinXYZ(deltacp_0_bin_glob, deltacp_0_bin_x, dummy, dummy)
                            dchisq_at_0.append(hchisq.GetBinContent(deltacp_0_bin_x))    # 0 is at the centre of this bin, so no need to extrapolate

                        # TODO: If there's no dchisq value in the hist > crit_dchisq, then just leave this row of the plot blank instead of dying
                        # Search forwards from left for the two bins directly greater than and less than the critical value:
                        ibin_high = -1  # chisq hist bin with content > critical value
                        ibin_low  = -1  # chisq hist bin with content < critical value
                        found = False
                        for ibin in range(1, hchisq.GetNbinsX()+1):
                            if hchisq.GetBinContent(ibin) > crit_dchisq:
                                ibin_high = ibin
                            elif hchisq.GetBinContent(ibin) <= crit_dchisq:
                                ibin_low = ibin
                                if ibin_high != -1:
                                    # Interpolate between bins
                                    range_low.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_high, ibin_low, var))
                                    found = True
                                    break
                                else:
                                    if var == "deltacp":
                                        # Try to search backwards from right instead
                                        break
                                    else:
                                        print "ERROR: No bins > critval, when searching from left. Exiting."
                                        sys.exit()
                        if var == "deltacp" and not found:
                            # Search backwards from right with swapped criteria:
                            ibin_high = -1
                            ibin_low  = -1
                            for ibin in range(hchisq.GetNbinsX(), 0, -1):
                                if hchisq.GetBinContent(ibin) <= crit_dchisq:
                                    ibin_low = ibin
                                elif hchisq.GetBinContent(ibin) > crit_dchisq:
                                    ibin_high = ibin
                                    if ibin_low != -1:
                                        # Interpolate between bins
                                        range_low.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_high, ibin_low, var))
                                        found = True
                                        wrap_range = True
                                        wrap_range_batch[ibatch] = True
                                        break
                                    else:
                                        # Maybe the bins we want are either side of the +/- pi boundary, so try that
                                        if hchisq.GetBinContent(1) <= crit_dchisq and hchisq.GetBinContent(hchisq.GetNbinsX()) > crit_dchisq:
                                            print "*** WARNING: dchisq hist crosses the crit dchisq at the +/- pi boundary. Behaviour may be strange as this is untested. ***"
                                            ibin_low  = 1
                                            ibin_high = hchisq.GetNbinsX()
                                            x1 = hchisq.GetBinCenter(ibin_high)    # high because treating it like moving the 1st bin onto the end of the last
                                            x2 = hchisq.GetBinCenter(ibin_high) + hchisq.GetBinWidth(ibin_low)
                                            y1 = hchisq.GetBinContent(ibin_high)
                                            y2 = hchisq.GetBinContent(ibin_low)
                                            range_low.append(interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2, var))
                                        else:
                                            print "ERROR: No bins < critval, when searching from right. Exiting."
                                            sys.exit()

                                    
                        # Search backwards from right for the two bins directly greater than and less than the critical value:
                        ibin_high = -1
                        ibin_low  = -1
                        found = False
                        for ibin in range(hchisq.GetNbinsX(), 0, -1):
                            if hchisq.GetBinContent(ibin) > crit_dchisq:
                                ibin_high = ibin
                            elif hchisq.GetBinContent(ibin) <= crit_dchisq:
                                ibin_low = ibin
                                if ibin_high != -1:
                                    # Interpolate between bins
                                    range_high.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_low, ibin_high, var))
                                    found = True
                                    break
                                else:
                                    if var == "deltacp":
                                        # Try to search forwards from left instead
                                        break
                                    else:
                                        print "ERROR: No bins > critval, when searching from right. Exiting."
                                        sys.exit()
                        if var == "deltacp" and not found:
                            # Search forwards from left with swapped criteria
                            if var == "deltacp":
                                ibin_high = -1
                                ibin_low  = -1
                                for ibin in range(1, hchisq.GetNbinsX()+1):
                                    if hchisq.GetBinContent(ibin) <= crit_dchisq:
                                        ibin_low = ibin
                                    elif hchisq.GetBinContent(ibin) > crit_dchisq:
                                        ibin_high = ibin
                                        if ibin_low != -1:
                                            # Interpolate between bins
                                            range_high.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_low, ibin_high, var))
                                            found = True
                                            wrap_range = True
                                            wrap_range_batch[ibatch] = True
                                            break
                                        else:
                                            # Maybe the bins we want are either side of the +/- pi boundary, so try that
                                            if hchisq.GetBinContent(1) > crit_dchisq and hchisq.GetBinContent(hchisq.GetNbinsX()) <= crit_dchisq:
                                                print "*** WARNING: dchisq hist crosses the crit dchisq at the +/- pi boundary. Behaviour may be strange as this is untested. ***"
                                                ibin_low  = hchisq.GetNbinsX()
                                                ibin_high = 1
                                                #x1 = hchisq.GetBinLowEdge(ibin_high) - hchisq.GetBinWidth(ibin_high)  # high - width because treating it like moving the last bin before the end of the 1st bin
                                                #x2 = hchisq.GetBinLowEdge(ibin_high) + hchisq.GetBinWidth(ibin_high)
                                                x1 = hchisq.GetBinCenter(ibin_high) - hchisq.GetBinWidth(ibin_high)  # high - width because treating it like moving the last bin before the end of the 1st bin
                                                x2 = hchisq.GetBinCenter(ibin_high)
                                                y1 = hchisq.GetBinContent(ibin_low)
                                                y2 = hchisq.GetBinContent(ibin_high)
                                                range_high.append(interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2, var))
                                                found = True
                                                wrap_range = True
                                                wrap_range_batch[ibatch] = True
                                            else:
                                                print "ERROR: No bins < critval, when searching from left. Exiting."
                                                sys.exit()
                            else:
                                print "ERROR: should not be here."
                                sys.exit()

                        if not found:
                            print "ERROR: Failed to find a confidence interval. Exiting."
                            sys.exit()

                    mixed_wrap = False
                    for tval in wrap_range_batch:
                        if tval != wrap_range:
                            mixed_wrap = True

                    # Print formatting
                    # if var == "dmsq23":
                    #     range_low  = [1000.0 * i for i in range_low]
                    #     range_high = [1000.0 * i for i in range_high]


                    # correct individual ranges
                    if var == "deltacp" and mixed_wrap:
                        # The confidence interval can also be wrong, so correct it here
                        # TODO: put all of these into functions to tidy up
                        sum_dist_to_left_edge = 0
                        sum_dist_to_right_edge = 0
                        for val in range_low:
                            sum_dist_to_left_edge  += ROOT.TMath.Abs(val + conv.kpi)
                            sum_dist_to_right_edge += ROOT.TMath.Abs(conv.kpi - val)
                        print "sum l =", sum_dist_to_left_edge
                        print "sum r =", sum_dist_to_right_edge
                        if sum_dist_to_left_edge < sum_dist_to_right_edge:
                            # Subtract 2pi from the values that are closer to the right edge
                            for iv in range(0, len(range_low)):
                                print "For val =", range_low[iv], ":"
                                dist_to_left_edge  = ROOT.TMath.Abs(range_low[iv] + conv.kpi)
                                dist_to_right_edge = ROOT.TMath.Abs(conv.kpi - range_low[iv])
                                print "  dist_to_left_edge =", dist_to_left_edge
                                print "  dist_to_right_edge =", dist_to_right_edge
                                if dist_to_right_edge < dist_to_left_edge:
                                    range_low[iv] = -conv.kpi + dist_to_right_edge
                                    print "  new val =", -conv.kpi + dist_to_right_edge
                        else:
                            # Add 2pi to the values that are closer to the left edge
                            for iv in range(0, len(range_low)):
                                print "For val =", range_low[iv], ":"
                                dist_to_left_edge  = ROOT.TMath.Abs(range_low[iv] + conv.kpi)
                                dist_to_right_edge = ROOT.TMath.Abs(conv.kpi - range_low[iv])
                                print "  dist_to_left_edge =", dist_to_left_edge
                                print "  dist_to_right_edge =", dist_to_right_edge
                                if dist_to_left_edge < dist_to_right_edge:
                                    range_low[iv] = conv.kpi + dist_to_left_edge
                                    print "  new val =", conv.kpi + dist_to_left_edge

                        sum_dist_to_left_edge = 0
                        sum_dist_to_right_edge = 0
                        for val in range_high:
                            sum_dist_to_left_edge  += ROOT.TMath.Abs(val + conv.kpi)
                            sum_dist_to_right_edge += ROOT.TMath.Abs(conv.kpi - val)
                        if sum_dist_to_left_edge < sum_dist_to_right_edge:
                            # Subtract 2pi from the values that are closer to the right edge
                            for iv in range(0, len(range_high)):
                                dist_to_left_edge  = ROOT.TMath.Abs(range_high[iv] + conv.kpi)
                                dist_to_right_edge = ROOT.TMath.Abs(conv.kpi - range_high[iv])
                                if dist_to_right_edge < dist_to_left_edge:
                                    range_high[iv] = -conv.kpi + dist_to_right_edge
                        else:
                            # Add 2pi to the values that are closer to the left edge
                            for iv in range(0, len(range_high)):
                                dist_to_left_edge  = ROOT.TMath.Abs(range_high[iv] + conv.kpi)
                                dist_to_right_edge = ROOT.TMath.Abs(conv.kpi - range_high[iv])
                                if dist_to_left_edge < dist_to_right_edge:
                                    range_high[iv] = conv.kpi + dist_to_left_edge

                    # Take the mean range of each batch rather than the range of the marginalised chisq surface, because
                    # we're interested in the stability of the confidence interval range at the number of toys used to
                    # produce each batch rather than at Nbatches * Ntoys_per_batch toys.
                    xrange_low_mean  = np.mean(range_low)
                    xrange_high_mean = np.mean(range_high)

                    # Now that we have the x range for each fit, calculate error
                    xrange_low_err_up   = max(range_low) - xrange_low_mean
                    xrange_low_err_down = xrange_low_mean - min(range_low)
                    xrange_high_err_up   = max(range_high) - xrange_high_mean
                    xrange_high_err_down = xrange_high_mean - min(range_high)
                    # xrange_low_err   = max(range_low) - min(range_low)
                    # xrange_high_err  = max(range_high) - min(range_high)
                    # xrange_low_err  = np.std(range_low, ddof=len(range_low)-1) / np.sqrt(len(range_low))
                    # xrange_high_err = np.std(range_high, ddof=len(range_high)-1) / np.sqrt(len(range_high))

                    # correct errors
                    if var == "deltacp" and mixed_wrap:
                        # Sometimes the different batches of fits have confidence intervals spread across
                        # the +/- pi boundary, so the error comes out as huge, this fixes it
                        if xrange_low_err_up > conv.kpi:
                            xrange_low_err_up = ROOT.TMath.Abs(xrange_low_err_up - 2*conv.kpi)
                        if xrange_low_err_down > conv.kpi:
                            xrange_low_err_down = ROOT.TMath.Abs(xrange_low_err_down - 2*conv.kpi)
                        if xrange_high_err_up > conv.kpi:
                            xrange_high_err_up = ROOT.TMath.Abs(xrange_high_err_up - 2*conv.kpi)
                        if xrange_high_err_down > conv.kpi:
                            xrange_high_err_down = ROOT.TMath.Abs(xrange_high_err_down - 2*conv.kpi)

                    # Now we have the mean high and low edges, for deltacp this can be slightly outside of the
                    # [-pi, +pi] range (due to corrections to each individual range above), so now correct the mean range
                    # NOTE: Must do this after calculating the errors
                    if var == "deltacp":
                        if xrange_low_mean < -conv.kpi:
                            xrange_low_mean += 2*conv.kpi
                        if xrange_low_mean > conv.kpi:
                            xrange_low_mean -= 2*conv.kpi
                        if xrange_high_mean < -conv.kpi:
                            xrange_high_mean += 2*conv.kpi
                        if xrange_high_mean > conv.kpi:
                            xrange_high_mean -= 2*conv.kpi

                    print "CI low edges:", range_low
                    print "CI high edges:", range_high
                    print "CI mean low edge:", xrange_low_mean
                    print "CI mean high edge:", xrange_high_mean
                    print "CI mean low edge err down =", xrange_low_err_down, ", up =", xrange_low_err_up
                    print "CI mean high edge err down =", xrange_high_err_down, ", up =", xrange_high_err_up
                    print "\n\n"

                    # Make a box for the mean range with error bars on each side. Collect all objects for now, draw later
                    x_arr, y_arr, xerrup_arr, xerrdown_arr, yerrup_arr, yerrdown_arr= array('d'), array('d'), array('d'), array('d'), array('d'), array('d')
                    x_arr.append(xrange_low_mean)
                    x_arr.append(xrange_high_mean)
                    y_arr.append(irow)
                    y_arr.append(irow)
                    xerrup_arr.append(xrange_low_err_up)
                    xerrup_arr.append(xrange_high_err_up)
                    xerrdown_arr.append(xrange_low_err_down)
                    xerrdown_arr.append(xrange_high_err_down)
                    yerrup_arr.append(0.0)
                    yerrup_arr.append(0.0)
                    yerrdown_arr.append(0.0)
                    yerrdown_arr.append(0.0)
                    graph_errors = ROOT.TGraphAsymmErrors(2, x_arr, y_arr, xerrdown_arr, xerrup_arr, yerrdown_arr, yerrup_arr)
                    graph_errors.SetMarkerStyle(1)  # Sets the marker to a single do so it can't be seen (not drawing the marker gets rid of the axis for some reason)
                    ROOT.SetOwnership(graph_errors, False)  # This stops the object being deleted from memory when it goes out of scope
                    mg.Add(graph_errors,"AP")

                    # Sometimes the errors go across the +/- pi boundary for deltacp, so more error bars need to be drawn to show this
                    if var == "deltacp":
                        new_err_low = -1
                        new_err_high = -1
                        if xrange_low_mean < xrange_high_mean:
                            if xrange_low_mean - xrange_low_err_down < -conv.kpi:
                                new_err_low = ROOT.TMath.Abs(xrange_low_mean - xrange_low_err_down + conv.kpi)
                            if xrange_high_mean + xrange_high_err_up > conv.kpi:
                                new_err_high = ROOT.TMath.Abs(xrange_high_mean + xrange_high_err_up - conv.kpi)
                        else:
                            if xrange_low_mean + xrange_low_err_up > conv.kpi:
                                new_err_low = ROOT.TMath.Abs(xrange_low_mean + xrange_low_err_up - conv.kpi)
                            if xrange_high_mean - xrange_high_err_down < -conv.kpi:
                                new_err_high = ROOT.TMath.Abs(xrange_high_mean - xrange_high_err_down + conv.kpi)
                        if new_err_low > 0 or new_err_high > 0:
                            npoints = 1
                            if new_err_low > 0 and new_err_high > 0:
                                npoints = 2
                            x_arr2, y_arr2, xerrup_arr2, xerrdown_arr2, yerrup_arr2, yerrdown_arr2 = array('d'), array('d'), array('d'), array('d'), array('d'), array('d')
                            if xrange_low_mean < xrange_high_mean:
                                if new_err_low > 0:
                                    x_arr2.append(conv.kpi)
                                    xerrup_arr2.append(0)
                                    xerrdown_arr2.append(new_err_low)
                                if new_err_high > 0:
                                    x_arr2.append(-conv.kpi)
                                    xerrup_arr2.append(new_err_high)
                                    xerrdown_arr2.append(0)
                            else:
                                if new_err_low > 0:
                                    x_arr2.append(-conv.kpi)
                                    xerrup_arr2.append(new_err_low)
                                    xerrdown_arr2.append(0)
                                if new_err_high > 0:
                                    x_arr2.append(+conv.kpi)
                                    xerrup_arr2.append(0)
                                    xerrdown_arr2.append(new_err_high)
                            y_arr2.append(irow)
                            yerrup_arr2.append(0.0)
                            yerrdown_arr2.append(0.0)
                            if new_err_low > 0 and new_err_high > 0:
                                y_arr2.append(irow)
                                yerrup_arr2.append(0.0)
                                yerrdown_arr2.append(0.0)
                            graph_errors2 = ROOT.TGraphAsymmErrors(npoints, x_arr2, y_arr2, xerrdown_arr2, xerrup_arr2, yerrdown_arr2, yerrup_arr2)
                            graph_errors2.SetMarkerStyle(1)  # Sets the marker to a single do so it can't be seen (not drawing the marker gets rid of the axis for some reason)
                            ROOT.SetOwnership(graph_errors2, False)  # This stops the object being deleted from memory when it goes out of scope
                            mg.Add(graph_errors2,"AP")


                    box_sublist = []
                    if xrange_low_mean < xrange_high_mean:
                        # Only need one box
                        b = ROOT.TBox(xrange_low_mean, irow-_box_half_width, xrange_high_mean, irow+_box_half_width)
                        b.SetFillColor(0)
                        b.SetLineStyle(1)
                        b.SetLineColor(1)
                        b.SetFillStyle(0)
                        ROOT.SetOwnership(b, False)  # This stops the object being deleted from memory when it goes out of scope
                        box_sublist.append(b)
                    else:
                        # Need two boxes
                        b1 = ROOT.TBox(xaxis_range_low, irow-_box_half_width, xrange_high_mean, irow+_box_half_width)
                        b1.SetFillColor(0)
                        b1.SetLineStyle(1)
                        b1.SetLineColor(1)
                        b1.SetFillStyle(0)
                        ROOT.SetOwnership(b1, False)  # This stops the object being deleted from memory when it goes out of scope
                        box_sublist.append(b1)
                        b2 = ROOT.TBox(xrange_low_mean, irow-_box_half_width, xaxis_range_high, irow+_box_half_width)
                        b2.SetFillColor(0)
                        b2.SetLineStyle(1)
                        b2.SetLineColor(1)
                        b2.SetFillStyle(0)
                        ROOT.SetOwnership(b2, False)  # This stops the object being deleted from memory when it goes out of scope
                        box_sublist.append(b2)
                    box_list.append(box_sublist)

                    title_box = ROOT.TPaveText(xaxis_range_low - _title_width - _title_offset, irow - _half_title_hight,
                                               xaxis_range_low - _title_offset               , irow + _half_title_hight)
                    title_str = '{:>2}'.format(hier_tex)  + ", " + '{:>6}'.format(type_tex) + ", " + \
                                '{:>7}'.format(alloc_tex) + ", " + '{:>3}'.format(pot_tex) + "E21 POT"
                    title_box.AddText(title_str)
                    title_box.SetFillColor(0)
                    title_box.SetFillStyle(0)
                    title_box.SetLineStyle(0)
                    title_box.SetLineColor(0)
                    ROOT.SetOwnership(title_box, False)  # This stops the object being deleted from memory when it goes out of scope
                    title_list.append(title_box)

                    # Add mean best fit point and error
                    # TODO: Apply all of the above relevant deltacp range and error checks to the BF point and BF error
                    bf_point_mean     = np.mean(bf_point)
                    bf_point_err_up   = max(bf_point) - bf_point_mean
                    bf_point_err_down = bf_point_mean - min(bf_point)
                    x_arr3, y_arr3, xerrup_arr3, xerrdown_arr3, yerrup_arr3, yerrdown_arr3= array('d'), array('d'), array('d'), array('d'), array('d'), array('d')
                    x_arr3.append(bf_point_mean)
                    y_arr3.append(irow)
                    xerrup_arr3.append(bf_point_err_up)
                    xerrdown_arr3.append(bf_point_err_down)
                    yerrup_arr3.append(0.0)
                    yerrdown_arr3.append(0.0)
                    graph_errors3 = ROOT.TGraphAsymmErrors(1, x_arr3, y_arr3, xerrdown_arr3, xerrup_arr3, yerrdown_arr3, yerrup_arr3)
                    graph_errors3.SetMarkerStyle(23)    # TODO: marker size may need automatically adjusting based on the number of rows
                    ROOT.SetOwnership(graph_errors3, False)  # This stops the object being deleted from memory when it goes out of scope
                    mg.Add(graph_errors3,"AP")

                    # For deltacp we need plots of the distance between the confidence interval edges and the 0 and +/- pi boundaries
                    if var == "deltacp":
                        if not (bf_point_mean < 0 and bf_point_mean > -conv.kpi):
                            print "WARNING: When drawing the distance to the deltacp boundaries, I've assumed that the " \
                                  "mean BF point lies between -pi and 0, but this is not the case for this scenario."
                        dist_to_pi          = None
                        dist_to_pi_err_down = None
                        dist_to_pi_err_up   = None
                        if ROOT.TMath.Abs(xrange_low_mean + conv.kpi) < ROOT.TMath.Abs(conv.kpi - xrange_low_mean):
                            dist_to_pi          = xrange_low_mean + conv.kpi
                            dist_to_pi_err_down = xrange_low_err_down
                            dist_to_pi_err_up   = xrange_low_err_up
                        else:
                            # Has crossed boundary
                            dist_to_pi          = xrange_low_mean - conv.kpi
                            dist_to_pi_err_down = xrange_low_err_up
                            dist_to_pi_err_up   = xrange_low_err_down
                        max_dist_to_pi = max(max_dist_to_pi, 1.05 * ROOT.TMath.Abs(dist_to_pi))

                        dist_to_0          = -xrange_high_mean
                        dist_to_0_err_down = None
                        dist_to_0_err_up   = None
                        if dist_to_0 > 0:
                            dist_to_0_err_down = xrange_high_err_up
                            dist_to_0_err_up   = xrange_high_err_down
                        else:
                            # Has crossed boundary
                            dist_to_0_err_down = xrange_high_err_down
                            dist_to_0_err_up   = xrange_high_err_up
                        max_dist_to_0 = max(max_dist_to_0, 1.05 * ROOT.TMath.Abs(dist_to_0))
                            
                        # Dist to pi graph
                        x_arr_pi, y_arr_pi, xerrdown_arr_pi, xerrup_arr_pi, yerrdown_arr_pi, yerrup_arr_pi = array('d'), array('d'), array('d'),  array('d'), array('d'), array('d')
                        x_arr_pi.append(dist_to_pi)
                        y_arr_pi.append(irow)
                        xerrdown_arr_pi.append(dist_to_pi_err_down)
                        xerrup_arr_pi.append(dist_to_pi_err_up)
                        yerrdown_arr_pi.append(0.0)
                        yerrup_arr_pi.append(0.0)
                        graph_errors_pi = ROOT.TGraphAsymmErrors(1, x_arr_pi, y_arr_pi, xerrdown_arr_pi, xerrup_arr_pi, yerrdown_arr_pi, yerrup_arr_pi)
                        graph_errors_pi.SetMarkerStyle(23)    # TODO: marker size may need automatically adjusting based on the number of rows
                        ROOT.SetOwnership(graph_errors_pi, False)  # This stops the object being deleted from memory when it goes out of scope
                        mg_dist_from_pi.Add(graph_errors_pi,"AP")
                        
                        # Dist to 0 graph
                        x_arr_0, y_arr_0, xerrdown_arr_0, xerrup_arr_0, yerrdown_arr_0, yerrup_arr_0 = array('d'), array('d'), array('d'),  array('d'), array('d'), array('d')
                        x_arr_0.append(dist_to_0)
                        y_arr_0.append(irow)
                        xerrdown_arr_0.append(dist_to_0_err_down)
                        xerrup_arr_0.append(dist_to_0_err_up)
                        yerrdown_arr_0.append(0.0)
                        yerrup_arr_0.append(0.0)
                        graph_errors_0 = ROOT.TGraphAsymmErrors(1, x_arr_0, y_arr_0, xerrdown_arr_0, xerrup_arr_0, yerrdown_arr_0, yerrup_arr_0)
                        graph_errors_0.SetMarkerStyle(23)    # TODO: marker size may need automatically adjusting based on the number of rows
                        ROOT.SetOwnership(graph_errors_0, False)  # This stops the object being deleted from memory when it goes out of scope
                        mg_dist_from_0.Add(graph_errors_0,"AP")


                    # For deltacp we need plots of the mean difference in Delta chisq between the critical value and the Delta chisq at +/- pi and 0
                    if var == "deltacp":
                        # Difference at deltacp = +/- pi
                        mean_dchisq_at_pi = np.mean(dchisq_at_pi)
                        mean_dchisq_at_pi_err_down = mean_dchisq_at_pi - min(dchisq_at_pi)  # This is the same as the error on the difference, as the critical dchisq has no error
                        mean_dchisq_at_pi_err_up   = max(dchisq_at_pi) - mean_dchisq_at_pi
                        diff_dchisq_at_pi = mean_dchisq_at_pi - crit_dchisq
                        max_ddchisq_pi = max(max_ddchisq_pi, 1.05 * diff_dchisq_at_pi)

                        x_arr_pi, y_arr_pi, xerrdown_arr_pi, xerrup_arr_pi, yerrdown_arr_pi, yerrup_arr_pi = array('d'), array('d'), array('d'),  array('d'), array('d'), array('d')
                        x_arr_pi.append(diff_dchisq_at_pi)
                        y_arr_pi.append(irow)
                        xerrdown_arr_pi.append(mean_dchisq_at_pi_err_down)
                        xerrup_arr_pi.append(mean_dchisq_at_pi_err_up)
                        yerrdown_arr_pi.append(0.0)
                        yerrup_arr_pi.append(0.0)
                        graph_errors_pi = ROOT.TGraphAsymmErrors(1, x_arr_pi, y_arr_pi, xerrdown_arr_pi, xerrup_arr_pi, yerrdown_arr_pi, yerrup_arr_pi)
                        graph_errors_pi.SetMarkerStyle(23)    # TODO: marker size may need automatically adjusting based on the number of rows
                        ROOT.SetOwnership(graph_errors_pi, False)  # This stops the object being deleted from memory when it goes out of scope
                        mg_ddchisq_pi.Add(graph_errors_pi,"AP")
                        
                        # Difference at deltacp = 0
                        mean_dchisq_at_0  = np.mean(dchisq_at_0)
                        mean_dchisq_at_0_err_down = mean_dchisq_at_0 - min(dchisq_at_0)
                        mean_dchisq_at_0_err_up   = max(dchisq_at_0) - mean_dchisq_at_0
                        diff_dchisq_at_0 = mean_dchisq_at_0 - crit_dchisq
                        max_ddchisq_0 = max(max_ddchisq_0, 1.05 * diff_dchisq_at_0)
                        
                        x_arr_0, y_arr_0, xerrdown_arr_0, xerrup_arr_0, yerrdown_arr_0, yerrup_arr_0 = array('d'), array('d'), array('d'),  array('d'), array('d'), array('d')
                        x_arr_0.append(diff_dchisq_at_0)
                        y_arr_0.append(irow)
                        xerrdown_arr_0.append(mean_dchisq_at_0_err_down)
                        xerrup_arr_0.append(mean_dchisq_at_0_err_up)
                        yerrdown_arr_0.append(0.0)
                        yerrup_arr_0.append(0.0)
                        graph_errors_0 = ROOT.TGraphAsymmErrors(1, x_arr_0, y_arr_0, xerrdown_arr_0, xerrup_arr_0, yerrdown_arr_0, yerrup_arr_0)
                        graph_errors_0.SetMarkerStyle(23)    # TODO: marker size may need automatically adjusting based on the number of rows
                        ROOT.SetOwnership(graph_errors_0, False)  # This stops the object being deleted from memory when it goes out of scope
                        mg_ddchisq_0.Add(graph_errors_0,"AP")


    # Plot the confidence intervals
    canvas = ROOT.TCanvas("ranges", "ranges", 1024, 768)
    canvas.SetLeftMargin(1.05*(_title_width - _title_offset)/(xaxis_range_high - xaxis_range_low + _title_width))
    canvas.SetRightMargin(0.01)
    mg.Draw("AP")
    canvas.Range(xaxis_range_low - _title_width - _title_offset, 0, xaxis_range_high + _title_offset, max_lines+1)

    #axis = ROOT.TGaxis(xaxis_range_low - _title_width - _title_offset, 0.5, xaxis_range_high, 0.5,
    #                   xaxis_range_low, xaxis_range_high, 7+100*5,"")
    #axis.SetTitle(x_title)
    #axis.Draw()
    canvas.Update()
    mg.GetYaxis().SetLabelOffset(999)   # Get rid of y-axis labels, as we've made our own
    mg.GetYaxis().SetLabelSize(0)
    mg.GetYaxis().SetTickLength(0)
    mg.GetYaxis().SetLimits(0, max_lines+1)
    mg.GetXaxis().SetLimits(xaxis_range_low, xaxis_range_high)
    mg.GetXaxis().SetTitle(x_title)
    canvas.Update()
    for i in range(0, len(title_list)):
        for box in box_list[i]:
            box.Draw("l,same")
        title_list[i].Draw("NB,same")
    canvas.Update()

    print "saving image to file", outfilename
    outfilename_use = outfilename
    if var == "deltacp":
        outfilename_use += "("
    canvas.SaveAs(outfilename_use, "recreate")


    # For deltacp, plot graphs of the distance of the CI edges to the 0 and +/- pi boundaries
    if var == "deltacp":
        # deltacp = +/- pi boundary plot
        canvas_pi = ROOT.TCanvas("dist_to_pi", "dist_to_pi", 1024, 768)
        canvas_pi.SetLeftMargin(1.05*(_title_width - _title_offset)/(xaxis_range_high - xaxis_range_low + _title_width))
        canvas_pi.SetRightMargin(0.01)

        mg_dist_from_pi.Draw("AP")
        canvas_pi.Range(-0.5 - _title_width - _title_offset, 0, max_dist_to_pi + _title_offset, max_lines+1)
        canvas_pi.Update()
        mg_dist_from_pi.GetYaxis().SetLabelOffset(999)   # Get rid of y-axis labels, as we've made our own
        mg_dist_from_pi.GetYaxis().SetLabelSize(0)
        mg_dist_from_pi.GetYaxis().SetTickLength(0)
        mg_dist_from_pi.GetYaxis().SetLimits(0, max_lines+1)
        mg_dist_from_pi.GetXaxis().SetLimits(-0.5, max_dist_to_pi)
        mg_dist_from_pi.GetXaxis().SetTitle("Distance of CI edge to the #delta_{CP} = #pm #pi boundary (radians)")

        for i in range(0, len(title_list)):
            title_list[i].Draw("NB,same")

        #pave_cross_pi = ROOT.TPaveText(-0.95*max_dist_to_pi, max_lines-3, -0.45*max_dist_to_pi, max_lines+1)
        pave_cross_pi = ROOT.TPaveText(-0.48, max_lines-2.5, 0, max_lines+0.5)
        pave_cross_pi.SetShadowColor(0)
        pave_cross_pi.SetFillColor(0)
        pave_cross_pi.SetFillStyle(0)
        pave_cross_pi.SetLineColor(0)
        pave_cross_pi.AddText("CI crosses #delta_{CP} = #pm #pi")
        pave_cross_pi.SetTextSize(0.02)
        pave_cross_pi.Draw("NB,same")

        #pave_no_cross_pi = ROOT.TPaveText(0.45*max_dist_to_pi, max_lines-3, 0.95*max_dist_to_pi, max_lines+1)
        pave_no_cross_pi = ROOT.TPaveText(0, max_lines-2.5, max_dist_to_pi*0.98, max_lines+0.5)
        pave_no_cross_pi.SetShadowColor(0)
        pave_no_cross_pi.SetFillColor(0)
        pave_no_cross_pi.SetFillStyle(0)
        pave_no_cross_pi.SetLineColor(0)
        pave_no_cross_pi.AddText("CI doesn't cross #delta_{CP} = #pm #pi")
        pave_no_cross_pi.SetTextSize(0.02)
        pave_no_cross_pi.Draw("NB,same")

        line_0 = ROOT.TLine(0, 0, 0, max_lines+1)
        line_0.SetLineStyle(2)
        line_0.Draw("same")
        canvas_pi.Update()
        canvas_pi.SaveAs(outfilename, "recreate")

        # deltacp = 0 boundary plot
        canvas_0 = ROOT.TCanvas("dist_to_0", "dist_to_0", 1024, 768)
        canvas_0.SetLeftMargin(1.05*(_title_width - _title_offset)/(xaxis_range_high - xaxis_range_low + _title_width))
        canvas_0.SetRightMargin(0.01)

        mg_dist_from_0.Draw("AP")
        canvas_0.Range(-0.5 - _title_width - _title_offset, 0, max_dist_to_0 + _title_offset, max_lines+1)
        canvas_0.Update()
        mg_dist_from_0.GetYaxis().SetLabelOffset(999)   # Get rid of y-axis labels, as we've made our own
        mg_dist_from_0.GetYaxis().SetLabelSize(0)
        mg_dist_from_0.GetYaxis().SetTickLength(0)
        mg_dist_from_0.GetYaxis().SetLimits(0, max_lines+1)
        mg_dist_from_0.GetXaxis().SetLimits(-0.5, max_dist_to_0)
        mg_dist_from_0.GetXaxis().SetTitle("Distance of CI edge to the #delta_{CP} = 0 boundary (radians)")

        for i in range(0, len(title_list)):
            title_list[i].Draw("NB,same")

        #pave_cross_0 = ROOT.TPaveText(-0.95*max_dist_to_0, max_lines-3, -0.45*max_dist_to_0, max_lines+1)
        pave_cross_0 = ROOT.TPaveText(-0.48, max_lines-2.5, 0, max_lines+0.5)
        pave_cross_0.SetShadowColor(0)
        pave_cross_0.SetFillColor(0)
        pave_cross_0.SetFillStyle(0)
        pave_cross_0.SetLineColor(0)
        pave_cross_0.AddText("CI crosses #delta_{CP} = 0")
        pave_cross_0.SetTextSize(0.02)
        pave_cross_0.Draw("NB,same")

        #pave_no_cross_0 = ROOT.TPaveText(0.45*max_dist_to_0, max_lines-3, 0.95*max_dist_to_0, max_lines+1)
        pave_no_cross_0 = ROOT.TPaveText(0, max_lines-2.5, max_dist_to_0*0.98, max_lines+0.5)
        pave_no_cross_0.SetShadowColor(0)
        pave_no_cross_0.SetFillColor(0)
        pave_no_cross_0.SetFillStyle(0)
        pave_no_cross_0.SetLineColor(0)
        pave_no_cross_0.AddText("CI doesn't cross #delta_{CP} = 0")
        pave_no_cross_0.SetTextSize(0.02)
        pave_no_cross_0.Draw("NB,same")

        line_0.Draw("same")
        canvas_0.Update()
        canvas_0.SaveAs(outfilename, "recreate")

    # For deltacp plot the mean difference in Delta chisq between the critical value and the Delta chisq at +/- pi and 0
    if var == "deltacp":
        min_ddchisq = -7
        # deltacp = +/- pi boundary plot
        canvas_pi = ROOT.TCanvas("ddchisq_pi", "ddchisq_pi", 1024, 768)
        canvas_pi.SetLeftMargin(1.05*(_title_width - _title_offset)/(xaxis_range_high - xaxis_range_low + _title_width))
        canvas_pi.SetRightMargin(0.01)

        mg_ddchisq_pi.Draw("AP")
        canvas_pi.Range(min_ddchisq - _title_width - _title_offset, 0, max_ddchisq_pi + _title_offset, max_lines+1)
        canvas_pi.Update()
        mg_ddchisq_pi.GetYaxis().SetLabelOffset(999)   # Get rid of y-axis labels, as we've made our own
        mg_ddchisq_pi.GetYaxis().SetLabelSize(0)
        mg_ddchisq_pi.GetYaxis().SetTickLength(0)
        mg_ddchisq_pi.GetYaxis().SetLimits(0, max_lines+1)
        mg_ddchisq_pi.GetXaxis().SetLimits(min_ddchisq, max_ddchisq_pi)
        mg_ddchisq_pi.GetXaxis().SetTitle("#Delta #chi^{2} #left( #delta_{CP} = #pm #pi #right) - #Delta #chi^{2}_{crit.}")

        for i in range(0, len(title_list)):
            title_list[i].Draw("NB,same")

        pave_cross_pi = ROOT.TPaveText(min_ddchisq*0.98, max_lines-2.5, 0, max_lines+0.5)
        pave_cross_pi.SetShadowColor(0)
        pave_cross_pi.SetFillColor(0)
        pave_cross_pi.SetFillStyle(0)
        pave_cross_pi.SetLineColor(0)
        pave_cross_pi.AddText("CI crosses #delta_{CP} = #pm #pi")
        pave_cross_pi.SetTextSize(0.02)
        pave_cross_pi.Draw("NB,same")

        pave_no_cross_pi = ROOT.TPaveText(0, max_lines-2.5, max_ddchisq_pi*0.98, max_lines+0.5)
        pave_no_cross_pi.SetShadowColor(0)
        pave_no_cross_pi.SetFillColor(0)
        pave_no_cross_pi.SetFillStyle(0)
        pave_no_cross_pi.SetLineColor(0)
        pave_no_cross_pi.AddText("CI doesn't cross #delta_{CP} = #pm #pi")
        pave_no_cross_pi.SetTextSize(0.02)
        pave_no_cross_pi.Draw("NB,same")

        line_0 = ROOT.TLine(0, 0, 0, max_lines+1)
        line_0.SetLineStyle(2)
        line_0.Draw("same")
        canvas_pi.Update()
        canvas_pi.SaveAs(outfilename, "recreate")

        # deltacp = 0 boundary plot
        canvas_0 = ROOT.TCanvas("ddchisq_0", "ddchisq_0", 1024, 768)
        canvas_0.SetLeftMargin(1.05*(_title_width - _title_offset)/(xaxis_range_high - xaxis_range_low + _title_width))
        canvas_0.SetRightMargin(0.01)

        mg_ddchisq_0.Draw("AP")
        canvas_0.Range(min_ddchisq - _title_width - _title_offset, 0, max_ddchisq_0 + _title_offset, max_lines+1)
        canvas_0.Update()
        mg_ddchisq_0.GetYaxis().SetLabelOffset(999)   # Get rid of y-axis labels, as we've made our own
        mg_ddchisq_0.GetYaxis().SetLabelSize(0)
        mg_ddchisq_0.GetYaxis().SetTickLength(0)
        mg_ddchisq_0.GetYaxis().SetLimits(0, max_lines+1)
        mg_ddchisq_0.GetXaxis().SetLimits(min_ddchisq, max_ddchisq_0)
        mg_ddchisq_0.GetXaxis().SetTitle("#Delta #chi^{2} #left( #delta_{CP} = 0 #right) - #Delta #chi^{2}_{crit.}")

        for i in range(0, len(title_list)):
            title_list[i].Draw("NB,same")

        pave_cross_0 = ROOT.TPaveText(min_ddchisq*0.98, max_lines-2.5, 0, max_lines+0.5)
        pave_cross_0.SetShadowColor(0)
        pave_cross_0.SetFillColor(0)
        pave_cross_0.SetFillStyle(0)
        pave_cross_0.SetLineColor(0)
        pave_cross_0.AddText("CI crosses #delta_{CP} = 0")
        pave_cross_0.SetTextSize(0.02)
        pave_cross_0.Draw("NB,same")

        pave_no_cross_0 = ROOT.TPaveText(0, max_lines-2.5, max_ddchisq_0*0.98, max_lines+0.5)
        pave_no_cross_0.SetShadowColor(0)
        pave_no_cross_0.SetFillColor(0)
        pave_no_cross_0.SetFillStyle(0)
        pave_no_cross_0.SetLineColor(0)
        pave_no_cross_0.AddText("CI doesn't cross #delta_{CP} = 0")
        pave_no_cross_0.SetTextSize(0.02)
        pave_no_cross_0.Draw("NB,same")

        line_0 = ROOT.TLine(0, 0, 0, max_lines+1)
        line_0.SetLineStyle(2)
        line_0.Draw("same")
        canvas_0.Update()
        canvas_0.SaveAs(outfilename+")", "recreate")

    print "Done!"


#_______________________________________________________________________________________________________________________
def interpolate_crit_x(crit_dchisq, hchisq, ibin1, ibin2, var):
    y1 = hchisq.GetBinContent(ibin1)
    y2 = hchisq.GetBinContent(ibin2)
    # x1 = hchisq.GetBinLowEdge(ibin1)
    # x2 = hchisq.GetBinLowEdge(ibin2) + hchisq.GetBinWidth(ibin2)
    x1 = hchisq.GetBinCenter(ibin1)
    x2 = hchisq.GetBinCenter(ibin2)
    return interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2, var)

#_______________________________________________________________________________________________________________________
def interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2, var):
    # Keep values within physical limits
    lim_low  = -conv.kALargeNum
    lim_high =  conv.kALargeNum
    #if var == "deltacp":
    #    lim_low  = -conv.kpi
    #    lim_high =  conv.kpi
    if "mixing" in var:
        lim_low  = 0
        lim_high = 1
    if x1 < lim_low:
        x1 = lim_low
    if x2 < lim_low:
        x2 = lim_low
    if x1 > lim_high:
        x1 = lim_high
    if x2 > lim_high:
        x2 = lim_high
    if y1 < 0:
        y1 = 0
    if y2 < 0:
        y2 = 0

    # Interpolate
    gradient = (y2 - y1) / (x2 - x1)
    intercept = y1 - gradient * x1
    return (crit_dchisq - intercept) / gradient

#_______________________________________________________________________________________________________________________
def get_crit_dchisq(critval_str_in):
    # Assume the input can either be "Xsigma", "x%" or "x". Treat the latter as a decimal.
    critval = -1
    is_sigma = False
    if "sigma" in critval_str_in:
        critval = float(critval_str_in.replace("sigma", ""))
        is_sigma = True
    elif "%" in critval_str_in:
        critval = float(critval_str_in.replace("%", ""))/100.
    else:
        critval = float(critval_str_in)

    if is_sigma:
        return critval**2
    else:
        return sigma(critval)**2

#_______________________________________________________________________________________________________________________
def sigma(frac):
    return ROOT.TMath.Sqrt(2) * ROOT.TMath.ErfInverse(frac)

#_______________________________________________________________________________________________________________________
def likelihood_mean_from_chisq_list(chisq_list):
    likelihood_list = [chisq_to_likelihood(i) for i in chisq_list]
    return np.mean(likelihood_list)

#_______________________________________________________________________________________________________________________
def chisq_to_likelihood(chisq):
    return ROOT.TMath.Exp(-0.5*chisq)

#_______________________________________________________________________________________________________________________
def likelihood_to_chisq(likelihood):
    return -2.0*ROOT.TMath.Log(likelihood)

#_______________________________________________________________________________________________________________________
if __name__=="__main__":
  
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..",description="Plot the confidence intervals for many difference PAC scenarios on the same graph")
    parser.add_option('-p', '--file-pattern', dest='file_pattern', default="batch[BATCH]/fit_oa2018_[MH]_reactor_[VAR]_nexp_1_0.003_-999.000-job0_[ALLOCATION]_[POT]_[DATATYPE]-_0_0_contour.root", help="File naming pattern for input ROOT files.")
    parser.add_option('-o', '--out-file'    , dest='output_file' , default="PAC_conf_interval_comparison.tex"                                                                                     , help="Name for output TeX file.")
    parser.add_option('-n', '--n-batches'   , dest='nbatch'      , default=3,   type=int                                                                                                          , help="The number of batches of fits to compare.")
    parser.add_option('-c', '--crit-val'    , dest='critval'     , default="4sigma"                                                                                                               , help="Critical value (in sigma or percent) to evaluate the range at.")
    parser.add_option('-v', '--var'         , dest='var'         , default="dmsq23"                                                                                                               , help="Oscillation param used in fit. \'dmsq23\' or \'mixing23\'")
    parser.add_option('-g', '--global'      , dest='do_global'   , action='store_true'                                                                                                               , help="Shift the chisq values of the IH contour relative to NH")
    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o,other = parser.parse_args()

    main(o.file_pattern, o.output_file, o.nbatch, o.critval, o.var, o.do_global)

    if other:
        print "Warning - ignored arguments: ",other
