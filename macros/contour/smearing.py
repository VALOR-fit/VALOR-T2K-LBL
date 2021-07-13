#!/usr/bin/env python

import argparse
import ROOT
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptTitle(0)

#_______________________________________________________________________________________________________________________
def catmull_rom(x):
    """Apply Catmull-Rom cubic interpolation"""
    B = 0.0
    C = 0.5
    x = np.abs(x)
    if x < 1:
        x2 = x * x
        x3 = x2 * x
        return ((12 - 9*B -6*C)*x3 + (-18 + 12*B + 6*C)*x2 + (6 - 2*B)) / 6.0
    elif 1 <= x < 2:
        x2 = x * x
        x3 = x2 * x
        return ((-B - 6*C)*x3 + (6*B + 30*C)*x2 + (-12*B - 48*C)*x + (8*B + 24*C)) / 6.0
    else:
        return 0.0

#_______________________________________________________________________________________________________________________
def sub_sample(input, row, col):
    """Exract an 4 x 4 sample around input[row][col]"""
    sample = np.zeros((4, 4))
    for dy in range(-1, 3):
        y = row + dy
        if y < 0:
            y = 0
        elif y >= input.shape[0]:
            y = input.shape[0] - 1
        for dx in range(-1, 3):
            x = col + dx
            if x < 0:
                x = 0
            elif x >= input.shape[1]:
                x = input.shape[1] - 1
            sample[dy + 1][dx + 1] = input[y][x]
    return sample

#_______________________________________________________________________________________________________________________
def sub_sample_1d(input, idx):
    """Exract an 4 x 4 sample around input[row][col]"""
    sample = np.zeros(4)
    for dx in range(-1, 3):
        x = idx + dx
        if x < 0:
            x = 0
        elif x >= len(input):
            x = len(input) - 1
        sample[dx + 1] = input[x]
    return sample

#_______________________________________________________________________________________________________________________
def make_filter(row, col, scale_row, scale_col):
    temp = (row - scale_row * 0.5) / scale_row
    pixel_row = int(temp)
    shift_row = temp - pixel_row
    temp = (col - scale_col * 0.5) / scale_col
    pixel_col = int(temp)
    shift_col = temp - pixel_col
    filter = np.zeros((4, 4))
    for dy in range(-1, 3):
        y = dy - shift_row
        for dx in range(-1, 3):
            x = -(dx - shift_col)
            filter[dy + 1][dx + 1] = catmull_rom(y) * catmull_rom(x)
    return filter

#_______________________________________________________________________________________________________________________
def make_filter_1d(idx, scale):
    temp = (idx - scale * 0.5) / scale
    pixel = int(temp)
    shift = temp - pixel
    filter = np.zeros(4)
    for dx in range(-1, 3):
        x = -(dx - shift)
        filter[dx + 1] = catmull_rom(x)
    return filter

#_______________________________________________________________________________________________________________________
def convolve(sub_sample, filter):
    return np.sum(sub_sample * filter)

#_______________________________________________________________________________________________________________________
def resample(input, scale):
    """Apply Catmull-Rom cubic interpolation"""
    if len(input.shape) == 1:
        # 1D resampling
        width = int(scale * len(input))
        output = np.zeros(width)
        for c in range(width):
            mc = int((c - scale * 0.5) / scale)
            sample = sub_sample_1d(input, mc)
            filter = make_filter_1d(c, scale)
            output[c] = convolve(sample, filter)
        return output
    elif len(input.shape) == 2:
        # 2D resampling
        width = int(scale * input.shape[1])
        height = int(scale * input.shape[0])
        new_shape =  (height, width)
        output = np.zeros(new_shape)
        for r in range(height):
            mr = int((r - scale * 0.5) / scale)
            for c in range(width):
                mc = int((c - scale * 0.5) / scale)
                sample = sub_sample(input, mr, mc)
                filter = make_filter(r, c, scale, scale)
                output[r][c] = convolve(sample, filter)
        return output
    else:
        print "Error, input %dD. 1D or 2D resampling only" % len(input.shape)
        sys.exit()

#_______________________________________________________________________________________________________________________
def smear(unsmeared, error):
    Ny = unsmeared.GetNbinsY()
    Nx = unsmeared.GetNbinsX()
    xaxis = unsmeared.GetXaxis()
    yaxis = unsmeared.GetYaxis()
    output = unsmeared.Clone()
    output.SetDirectory(0)
    for i in range(1, Nx + 1):
        for j in range(1, Ny + 1):
            y_u = yaxis.GetBinCenter(j)
            kmin = yaxis.FindBin(y_u - 5 * error)
            kmax = yaxis.FindBin(y_u + 5 * error)
            content = 0.0
            for k in range(kmin, kmax + 1):
                mean = yaxis.GetBinCenter(k)
                content += unsmeared.GetBinContent(i, k) * ROOT.TMath.Gaus(y_u, mean, error, False)
            content = content * xaxis.GetBinWidth(i) * yaxis.GetBinWidth(j) / (2 * np.pi * error)
            output.SetBinContent(i, j, content)
    return output

#_______________________________________________________________________________________________________________________
def smear_1d(unsmeared, error):
    N = unsmeared.GetNbinsX()
    xaxis = unsmeared.GetXaxis()
    output = unsmeared.Clone()
    output.SetDirectory(0)
    for i in range(1, N + 1):
        x_u = xaxis.GetBinCenter(i)
        content = 0.0
        jmin = xaxis.FindBin(x_u - 5 * error)
        jmax = xaxis.FindBin(x_u + 5 * error)
        for j in range(jmin, jmax + 1):
            mean = xaxis.GetBinCenter(j)
            content += unsmeared.GetBinContent(j) * ROOT.TMath.Gaus(x_u, mean, error, False)
        content = content * xaxis.GetBinWidth(i) / (2 * np.pi * error)
        output.SetBinContent(i, content)
    return output

#_______________________________________________________________________________________________________________________
def chisq_to_likelihood(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output.SetBinContent(i, j, np.exp(-0.5 * input.GetBinContent(i, j)))
    return output

#_______________________________________________________________________________________________________________________
def likelihood_to_chisq(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output.SetBinContent(i, j, -2.0 * np.log(input.GetBinContent(i, j)))
    return output

#_______________________________________________________________________________________________________________________
def make_projection_x(input):
    nbins_x = input.GetNbinsX()
    nbins_y = input.GetNbinsY()
    xaxis = input.GetXaxis()
    hchisq = ROOT.TH1D("hchi2_x", "hchi2_x", nbins_x, xaxis.GetBinLowEdge(1), xaxis.GetBinUpEdge(nbins_x))
    hchisq.SetDirectory(0)

    for i in range(1, nbins_x + 1):
        Lavg = 0.0
        for j in range(1 , nbins_y + 1):
            Lavg += input.GetBinContent(i, j)
        Lavg /= nbins_y
        hchisq.SetBinContent(i, -2 * ROOT.TMath.Log(Lavg))

    return hchisq

#_______________________________________________________________________________________________________________________
def make_projection_y(input):
    nbins_x = input.GetNbinsX()
    nbins_y = input.GetNbinsY()
    yaxis = input.GetYaxis()
    hchisq = ROOT.TH1D("hchi2_y", "hchi2_y", nbins_y, yaxis.GetBinLowEdge(1), yaxis.GetBinUpEdge(nbins_y))
    hchisq.SetDirectory(0)

    for i in range(1, nbins_y + 1):
        Lavg = 0.0
        for j in range(1 , nbins_x + 1):
            Lavg += input.GetBinContent(j, i)
        Lavg /= nbins_x
        hchisq.SetBinContent(i, -2 * ROOT.TMath.Log(Lavg))

    return hchisq

#_______________________________________________________________________________________________________________________
def chisq_to_likelihood_1d(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    for i in range(1, NXbins + 1):
        output.SetBinContent(i, np.exp(-0.5 * input.GetBinContent(i)))
    return output

#_______________________________________________________________________________________________________________________
def likelihood_to_chisq_1d(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    for i in range(1, NXbins + 1):
        output.SetBinContent(i, -2.0 * np.log(input.GetBinContent(i)))
    return output

#_______________________________________________________________________________________________________________________
def histo_to_nparray(input):
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    output = np.zeros((NYbins, NXbins))
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output[j - 1][i - 1] = input.GetBinContent(i, j)
    return output

#_______________________________________________________________________________________________________________________
def copy_nparray_to_histo(input_array, output_histo):
    NXbins = output_histo.GetNbinsX()
    NYbins = output_histo.GetNbinsY()
    if NXbins != input_array.shape[1] or NYbins != input_array.shape[0]:
        print "Error: Shape mismatch between array and histogram:"
        print "Array %dx%d Histogram %d%d (row x column)" % (input_array.shape[0], input_array.shape[1], NYbins, NXbins)
        sys.exit()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output_histo.SetBinContent(i, j, input_array[j - 1][i - 1])
    return output_histo

#_______________________________________________________________________________________________________________________
def histo_to_nparray_1d(input):
    NXbins = input.GetNbinsX()
    output = np.zeros(NXbins)
    for i in range(1, NXbins + 1):
        output[i - 1] = input.GetBinContent(i)
    return output

#_______________________________________________________________________________________________________________________
def copy_nparray_to_histo_1d(input_array, output_histo):
    NXbins = output_histo.GetNbinsX()
    if NXbins != len(input_array):
        print "Error: Shape mismatch between array and histogram:"
        print "Array %d Histogram %d" % (len(input_array), NXbins)
        sys.exit()
    for i in range(1, NXbins + 1):
        output_histo.SetBinContent(i, input_array[i - 1])
    return output_histo

#_______________________________________________________________________________________________________________________
def process_2d(original, root_out_file, error):
    xaxis = original.GetXaxis()
    yaxis = original.GetYaxis()

    original_array = histo_to_nparray(original)
    scale = 9.0
    resampled_array = resample(original_array, scale)
    resampled_histogram = ROOT.TH2D("atmos", "atmos", int(scale * original.GetNbinsX()), original.GetXaxis().GetXmin(), original.GetXaxis().GetXmax(), int(scale * original.GetNbinsY()), original.GetYaxis().GetXmin(), original.GetYaxis().GetXmax())
    copy_nparray_to_histo(resampled_array, resampled_histogram)
    likelihood_histogram = chisq_to_likelihood(resampled_histogram)
    likelihood_histogram.SetTitle("Resampled likelihood surface")

    smeared_histogram = smear(likelihood_histogram, error)
    smeared_histogram.GetXaxis().SetTitle("sin^{2}#theta_{23}")
    if "nh" in root_out_file.GetName() or "normal" in root_out_file.GetName():
        smeared_histogram.GetYaxis().SetTitle("|#Deltam^{2}_{32}| (eV^{2}c^{-4})")
    else:
        smeared_histogram.GetYaxis().SetTitle("|#Deltam^{2}_{13}| (eV^{2}c^{-4})")
    smeared_histogram.GetXaxis().SetTitleOffset(1.1)
    textsize = 0.05
    smeared_histogram.GetXaxis().SetTitleSize(textsize)
    smeared_histogram.GetYaxis().SetTitleOffset(1.3)
    smeared_histogram.GetYaxis().SetTitleSize(textsize)
    smeared_histogram.SetTitle("")

    chisq_histogram = likelihood_to_chisq(smeared_histogram)

    root_out_file.cd()
    xaxis_label = xaxis.GetTitle()
    yaxis_label = yaxis.GetTitle()
    smeared_histogram.Write("likelihood")
    chisq_histogram.GetXaxis().SetTitle(xaxis_label)
    chisq_histogram.GetYaxis().SetTitle(yaxis_label)
    chisq_histogram.Write("chi2")

    # Get the 1D projections
    chisq_x_histogram = make_projection_x(smeared_histogram)
    chisq_x_histogram.GetXaxis().SetTitle(xaxis_label)
    chisq_x_histogram.GetYaxis().SetTitle("#chi^{2}")
    chisq_x_histogram.Write("chi2_x")

    chisq_y_histogram = make_projection_y(smeared_histogram)
    chisq_y_histogram.GetXaxis().SetTitle(yaxis_label)
    chisq_y_histogram.GetYaxis().SetTitle("#chi^{2}")
    chisq_y_histogram.Write("chi2_y")

#_______________________________________________________________________________________________________________________
def process_1d(original, root_out_file, error):
    xaxis = original.GetXaxis()

    original_array = histo_to_nparray_1d(original)
    scale = 5.0
    resampled_array = resample(original_array, scale)
    resampled_histogram = ROOT.TH1D("atmos", "atmos", int(scale * original.GetNbinsX()), original.GetXaxis().GetXmin(), original.GetXaxis().GetXmax())
    copy_nparray_to_histo_1d(resampled_array, resampled_histogram)
    likelihood_histogram = chisq_to_likelihood_1d(resampled_histogram)

    smeared_histogram = smear_1d(likelihood_histogram, error)

    chisq_histogram = likelihood_to_chisq_1d(smeared_histogram)

    root_out_file.cd()
    xaxis_label = xaxis.GetTitle()
    smeared_histogram.GetXaxis().SetTitle(xaxis_label)
    smeared_histogram.Write("likelihood")
    chisq_histogram.GetXaxis().SetTitle(xaxis_label)
    chisq_histogram.Write("h_chisq")

#_______________________________________________________________________________________________________________________
if __name__=="__main__":
    ROOT.gStyle.SetOptStat(0)
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input contour file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output root file.")
    parser.add_argument("-e", "--error", dest="error", type=float, help="Smearing error to apply.")
    args = parser.parse_args()


    in_file       = ROOT.TFile(args.in_filename , "READ")
    root_out_file = ROOT.TFile(args.out_filename, "RECREATE")


    if in_file.GetListOfKeys().Contains("chi2"):
        original = in_file.Get("chi2")
        original.SetDirectory(0)
        if type(original) == ROOT.TH2D:
            process_2d(original, root_out_file, args.error)
        else:
            process_1d(original, root_out_file, args.error)

    # Copy over unchanged stuff
    in_file.cd()
    vec_list = list(in_file.Get("chi2_best_fit_toy"))
    vec = ROOT.vector('double')()
    for i in vec_list: vec.push_back(i)
    root_out_file.cd()
    root_out_file.WriteObject(vec, "chi2_best_fit_toy")

    in_file.Close()
    root_out_file.Close()

    print "Done!"

