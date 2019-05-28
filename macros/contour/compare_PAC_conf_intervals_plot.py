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
_title_offset     = 0.02
_half_title_hight = 0.25


#_______________________________________________________________________________________________________________________
def main(file_pattern, outfilename, nbatch, critval_str, var):

    crit_dchisq = get_crit_dchisq(critval_str)
    print "*** test: crit_dchisq =", crit_dchisq

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

                    if type == "pure_AsimovA" and alloc == "FHC" and hier == "normal" and pot == "10E21" and var == "deltacp":
                        print "\n\n\n *** WARNING: Skipping scenario pure_AsimovA, FHC, normal, 10E21, deltacp - it needs fixing! *** \n\n\n"
                        continue

                    range_low  = []
                    range_high = []
                    wrap_range = False
                    wrap_range_batch = [False for i in range(0, nbatch)]
                    for ibatch in range(0, nbatch):
                        file_in = file_pattern.replace("[BATCH]"     , str(ibatch))
                        file_in = file_in.replace     ("[DATATYPE]"  , type)
                        file_in = file_in.replace     ("[ALLOCATION]", alloc)
                        file_in = file_in.replace     ("[MH]"        , hier)
                        file_in = file_in.replace     ("[POT]"       , pot)
                        file_in = file_in.replace     ("[VAR]"       , var)
                        print "Reading from file", file_in
                        f = ROOT.TFile(file_in, "read")
                        hchisq = f.Get("h_deltachisq")
                        
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
                                    range_low.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_high, ibin_low))
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
                                        range_low.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_high, ibin_low))
                                        found = True
                                        wrap_range = True
                                        wrap_range_batch[ibatch] = True
                                        break
                                    else:
                                        # Maybe the bins we want are either side of the +/- pi boundary, so try that
                                        if hchisq.GetBinContent(1) <= crit_dchisq and hchisq.GetBinContent(hchisq.GetNbinsX()) > crit_dchisq:
                                            ibin_low  = 1
                                            ibin_high = hchisq.GetNbinsX()
                                            x1 = hchisq.GetBinLowEdge(ibin_high)    # high because treating it like moving the 1st bin onto the end of the last
                                            x2 = hchisq.GetBinLowEdge(ibin_high) + 2*hchisq.GetBinWidth(ibin_low)
                                            y1 = hchisq.GetBinContent(ibin_high)
                                            y2 = hchisq.GetBinContent(ibin_low)
                                            range_low.append(interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2))
                                        else:
                                            print "ERROR: No bins < critval, when searching from right. Exiting."
                                            sys.exit()

                                    
                        # Search backwards from right for the two bins directly greater than and less than the critical value:
                        ibin_high = -1
                        ibin_low  = -1
                        found = False
                        if not wrap_range:
                            for ibin in range(hchisq.GetNbinsX(), 0, -1):
                                if hchisq.GetBinContent(ibin) > crit_dchisq:
                                    ibin_high = ibin
                                elif hchisq.GetBinContent(ibin) <= crit_dchisq:
                                    ibin_low = ibin
                                    if ibin_high != -1:
                                        # Interpolate between bins
                                        range_high.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_low, ibin_high))
                                        found = True
                                        break
                                    else:
                                        print "ERROR: No bins > critval, when searching from right. Exiting."
                                        sys.exit()
                        else:
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
                                            range_high.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_low, ibin_high))
                                            found = True
                                            break
                                        else:
                                            # Maybe the bins we want are either side of the +/- pi boundary, so try that
                                            if hchisq.GetBinContent(1) > crit_dchisq and hchisq.GetBinContent(hchisq.GetNbinsX()) <= crit_dchisq:
                                                ibin_low  = hchisq.GetNbinsX()
                                                ibin_high = 1
                                                x1 = hchisq.GetBinLowEdge(ibin_high) - hchisq.GetBinWidth(ibin_high)  # high - width because treating it like moving the last bin before the end of the 1st bin
                                                x2 = hchisq.GetBinLowEdge(ibin_high) + hchisq.GetBinWidth(ibin_high)
                                                y1 = hchisq.GetBinContent(ibin_low)
                                                y2 = hchisq.GetBinContent(ibin_high)
                                                range_low.append(interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2))
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

                    print "*** test: found low edges:", range_low
                    print "*** test: found high edges:", range_high
                    print "*** test: wrap_range =", wrap_range
                    print "*** test: wrap_range_batch =", wrap_range_batch
                    print "*** test: mixed_wrap =", mixed_wrap

                    # Print formatting
                    # if var == "dmsq23":
                    #     range_low  = [1000.0 * i for i in range_low]
                    #     range_high = [1000.0 * i for i in range_high]

                    # Now that we have the x range for each fit, calculate error (the mean is calculated after deltacp corrections)
                    xrange_low_err   = max(range_low) - min(range_low)
                    xrange_high_err  = max(range_high) - min(range_high)
                    # xrange_low_err  = np.std(range_low, ddof=len(range_low)-1) / np.sqrt(len(range_low))
                    # xrange_high_err = np.std(range_high, ddof=len(range_high)-1) / np.sqrt(len(range_high))

                    # correct ranges and errors
                    if var == "deltacp" and mixed_wrap:
                        print "*** test: HERE1"
                        # Sometimes the different batches of fits have confidence intervals spread across
                        # the +/- pi boundary, so the error comes out as huge, this fixes it
                        if xrange_low_err > conv.kpi:
                            xrange_low_err = ROOT.TMath.Abs(xrange_low_err - 2*conv.kpi)
                        if xrange_high_err > conv.kpi:
                            xrange_high_err = ROOT.TMath.Abs(xrange_high_err - 2*conv.kpi)

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

                    xrange_low_mean  = np.mean(range_low)
                    xrange_high_mean = np.mean(range_high)

                    # Make a box for the mean range with error bars on each side. Collect all objects for now, draw later
                    x_arr, y_arr, xerr_arr, yerr_arr = array('d'), array('d'), array('d'), array('d')
                    x_arr.append(xrange_low_mean)
                    x_arr.append(xrange_high_mean)
                    y_arr.append(irow)
                    y_arr.append(irow)
                    xerr_arr.append(xrange_low_err)
                    xerr_arr.append(xrange_high_err)
                    yerr_arr.append(0.0)
                    yerr_arr.append(0.0)
                    graph_errors = ROOT.TGraphErrors(2, x_arr, y_arr, xerr_arr, yerr_arr)
                    graph_errors.SetMarkerStyle(1)  # Sets the marker to a single do so it can't be seen (not drawing the marker gets rid of the axis for some reason)
                    ROOT.SetOwnership(graph_errors, False)  # This stops the object being deleted from memory when it goes out of scope
                    mg.Add(graph_errors,"AP")

                    box_sublist = []
                    if not wrap_range:
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
                    title_box.AddText(hier_tex +" "+ type_tex +" "+ alloc_tex +" "+ pot_tex+"E21 POT")
                    title_box.SetFillColor(0)
                    title_box.SetFillStyle(0)
                    title_box.SetLineStyle(0)
                    title_box.SetLineColor(0)
                    ROOT.SetOwnership(title_box, False)  # This stops the object being deleted from memory when it goes out of scope
                    title_list.append(title_box)

    # Plot everything
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
    canvas.SaveAs(outfilename, "recreate")
    print "Done!"


#_______________________________________________________________________________________________________________________
def interpolate_crit_x(crit_dchisq, hchisq, ibin1, ibin2):
    y1 = hchisq.GetBinContent(ibin1)
    y2 = hchisq.GetBinContent(ibin2)
    x1 = hchisq.GetBinLowEdge(ibin1)
    x2 = hchisq.GetBinLowEdge(ibin2) + hchisq.GetBinWidth(ibin2)
    gradient = (y2 - y1) / (x2 - x1)
    intercept = y1 - gradient * x1
    return (crit_dchisq - intercept) / gradient

#_______________________________________________________________________________________________________________________
def interpolate_crit_x_no_hist(crit_dchisq, x1, x2, y1, y2):
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
if __name__=="__main__":
  
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")
    parser.add_option('-p', '--file-pattern', dest='file_pattern', default="batch[BATCH]/fit_oa2018_[MH]_reactor_[VAR]_nexp_1_0.003_-999.000-job0_[ALLOCATION]_[POT]_[DATATYPE]-_0_0_contour.root", help="File naming pattern for input ROOT files.")
    parser.add_option('-o', '--out-file'    , dest='output_file' , default="PAC_conf_interval_comparison.tex"                                                                                      , help="Name for output TeX file.")
    parser.add_option('-n', '--n-batches'   , dest='nbatch'      , default=3,   type=int                                                                                                           , help="The number of batches of fits to compare.")
    parser.add_option('-c', '--crit-val'    , dest='critval'     , default="4sigma"                                                                                                                , help="Critical value (in sigma or percent) to evaluate the range at.")
    parser.add_option('-v', '--var'         , dest='var'         , default="dmsq23"                                                                                                                , help="Oscillation param used in fit. \'dmsq23\' or \'mixing23\'")

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o,other = parser.parse_args()

    main(o.file_pattern, o.output_file, o.nbatch, o.critval, o.var)

    if other:
        print "Warning - ignored arguments: ",other
