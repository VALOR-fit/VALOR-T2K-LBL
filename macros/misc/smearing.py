#!/usr/bin/env python

import argparse
import ROOT
import array
import numpy as np
import sys

ROOT.gROOT.SetBatch(1)

def make_erec_spectrum(erec_theta):
    Nerec = erec_theta.GetNbinsY()
    erec_axis = erec_theta.GetYaxis()
    erec_low = erec_axis.GetBinLowEdge(1)
    erec_high = erec_axis.GetBinUpEdge(Nerec)
    Ntheta = erec_theta.GetNbinsX()
    theta_axis = erec_theta.GetXaxis()
    theta_low = theta_axis.GetBinLowEdge(1)
    theta_high = theta_axis.GetBinUpEdge(Ntheta)
    erec_name = erec_theta.GetName() + "_erec"
    erec = ROOT.TH1D(erec_name, erec_name, Nerec, erec_low, erec_high)
    for e in range(1, Nerec + 1):
        for t in range(1, Ntheta + 1):
            erec.SetBinContent(e, erec.GetBinContent(e) + erec_theta.GetBinContent(t, e))
    return erec

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

def make_filter_1d(idx, scale):
    temp = (idx - scale * 0.5) / scale
    pixel = int(temp)
    shift = temp - pixel
    filter = np.zeros(4)
    for dx in range(-1, 3):
        x = -(dx - shift)
        filter[dx + 1] = catmull_rom(x)
    return filter


def convolve(sub_sample, filter):
    return np.sum(sub_sample * filter)

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

def smear(input, error=3.556768e-5):
    # Smear the input TH2D - assumes smear variable is Y
    NXbins = input.GetNbinsX()
    xaxis = input.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)

    NYbins = input.GetNbinsY()
    yaxis = input.GetYaxis()
    y_low = yaxis.GetBinLowEdge(1)
    y_high = yaxis.GetBinUpEdge(NYbins)
    bin_width = yaxis.GetBinWidth(1)

    norm = 1.0 / (error * np.sqrt(2 * np.pi))
    max = float(-np.inf)
    output = input.Clone()
    output.SetDirectory(0)
    for i in range(0, NXbins):
        for j in range(0, NYbins):
            val = input.GetBinContent(i + 1, j + 1)
            if val > max:
                max = val
                max_x = i
                max_y = j

    canvas = ROOT.TCanvas("c2", "c2", 700, 500)
    canvas.Print("gaus.pdf[")

    like_slice = ROOT.TH1D("like_slice", "like_slice", NYbins, y_low, y_high)
    smear_contrib = ROOT.TH1D("smear_contrib", "smear_contrib", NYbins, y_low, y_high)
    bestfit_slice = ROOT.TH1D("like_slice", "like_slice", NYbins, y_low, y_high)
    boundary = ROOT.TH1D("boundary", "boundary", NYbins, y_low, y_high)
    left = yaxis.FindBin(0.00238)
    right = yaxis.FindBin(0.00265)
    for i in range(left - 7, left):
        boundary.SetBinContent(i, input.GetBinContent(max_x, i))
    for i in range(right + 1, right + 8):
        boundary.SetBinContent(i, input.GetBinContent(max_x, i))
    boundary.SetTitle("Likelihood around original 90% CL")
    boundary.Draw()
    canvas.Print("gaus.pdf")

    receive = []
    for i in range(0, NYbins):
        receive.append(ROOT.TH1D("receive" + str(i), "receive" + str(i), NYbins, y_low, y_high))
    for i in range(0, NXbins):
        x_smear = xaxis.GetBinCenter(i + 1)
        if i == max_x:
            for j in range(0, NYbins):
                bestfit_slice.SetBinContent(j + 1, input.GetBinContent(i + 1, j + 1))
        elif 0.4 < x_smear < 0.65:
            for j in range(0, NYbins):
                like_slice.SetBinContent(j + 1, input.GetBinContent(i + 1, j + 1))
            like_slice.SetTitle("Likelihood through sin^{2}#theta_{23} = " + str(x_smear) + " (" + str(i) + ")")
            like_slice.Draw()
            canvas.Print("gaus.pdf")

        for j in range(0, NYbins):
            y_smear = yaxis.GetBinCenter(j + 1)
            likelihood = input.GetBinContent(i + 1, j + 1)
            gauss_bin = ROOT.TF1("gauss_bin", "gaus(0)", y_low, y_high)
            gauss_bin.SetParameters(norm, y_smear, error)
            for k in range(0, NYbins):
                if j != k:
                    y_temp = yaxis.GetBinCenter(k + 1)
                    integral = gauss_bin.Integral(y_temp - 0.5 * bin_width, y_temp + 0.5 * bin_width)
                    output.Fill(x_smear, y_temp, integral * likelihood)
                    if i == max_x:
                        smear_contrib.SetBinContent(k, integral * likelihood)
                        receive[k - 1].SetBinContent(j, integral * likelihood)
                else:
                    if i == max_x:
                        smear_contrib.SetBinContent(k, 0)
            if i == max_x:
                smear_contrib.SetTitle("Spread of " + str(x_smear) + " " + str(y_smear))
                smear_contrib.Draw()
                canvas.Print("gaus.pdf")
    for i in range(1, NYbins + 1):
        receive[i - 1].SetTitle("Inputs into " + str(xaxis.GetBinCenter(max_x)) + " " + str(yaxis.GetBinCenter(i)))
        receive[i - 1].Draw()
        canvas.Print("gaus.pdf")
    bestfit_slice.SetTitle("Likelihood through sin^{2}#theta_{23} best fit")
    bestfit_slice.Draw()
    opt = "SAME"
    for i in (left, right, right + 1):
        receive[i - 1].Scale(10)
        receive[i - 1].SetLineColor(ROOT.kRed)
        receive[i - 1].Draw(opt)
        receive[i - 1].SetTitle("Gaussian at 90CL")

    canvas.Print("gaus.pdf")

    canvas.Print("gaus.pdf]")
    return output

def smear_1d(input, error=3.556768e-5):
    # Smear the input TH1D
    NXbins = input.GetNbinsX()
    xaxis = input.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)
    bin_width = xaxis.GetBinWidth(1)

    norm = 1.0 / (error * np.sqrt(2 * np.pi))
    output = input.Clone()
    output.SetDirectory(0)

    for i in range(0, NXbins):
        x_smear = xaxis.GetBinCenter(i + 1)
        likelihood = input.GetBinContent(i + 1)
        gauss_bin = ROOT.TF1("gauss_bin", "gaus(0)", x_low, x_high)
        gauss_bin.SetParameters(norm, x_smear, error)
        for k in range(0, NXbins):
            if i != k:
                x_temp = xaxis.GetBinCenter(k + 1)
                integral = gauss_bin.Integral(x_temp - 0.5 * bin_width, x_temp + 0.5 * bin_width)
                output.Fill(x_temp, integral * likelihood)
    return output


def marginalise(input, axis="x"):
    # Smear the input TH2D - assumes smear variable is Y
    NXbins = input.GetNbinsX()
    xaxis = input.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)

    NYbins = input.GetNbinsY()
    yaxis = input.GetYaxis()
    y_low = yaxis.GetBinLowEdge(1)
    y_high = yaxis.GetBinUpEdge(NYbins)

    if axis == "x":
        output = ROOT.TH1D("dmsq23", "dmsq23", NYbins, yaxis.GetXmin(), yaxis.GetXmax())
        for r in range(1, NYbins + 1):
            val = 0.0
            for c in range(1, NXbins + 1):
                val += input.GetBinContent(c, r)
            val /= NXbins
            output.SetBinContent(r, val)
    elif axis == "y":
        output = ROOT.TH1D("mixing23", "mixing23", NXbins, xaxis.GetXmin(), xaxis.GetXmax())
        for c in range(1, NXbins + 1):
            val = 0.0
            for r in range(1, NYbins + 1):
                val += input.GetBinContent(c, r)
            val /= NYbins
            output.SetBinContent(c, val)

    return output

def chisq_to_likelihood(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output.SetBinContent(i, j, np.exp(-0.5 * input.GetBinContent(i, j)))
    return output

def likelihood_to_chisq(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output.SetBinContent(i, j, -2.0 * np.log(input.GetBinContent(i, j)))
    return output

def make_delta_chisq(input):
    output = input.Clone()
    output.SetDirectory(0)
    minval = input.GetMinimum()
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output.SetBinContent(i, j, input.GetBinContent(i, j) - minval)
    return output

def chisq_to_likelihood_1d(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    for i in range(1, NXbins + 1):
        output.SetBinContent(i, np.exp(-0.5 * input.GetBinContent(i)))
    return output

def likelihood_to_chisq_1d(input):
    output = input.Clone()
    output.SetDirectory(0)
    NXbins = input.GetNbinsX()
    for i in range(1, NXbins + 1):
        output.SetBinContent(i, -2.0 * np.log(input.GetBinContent(i)))
    return output

def make_delta_chisq_1d(input):
    output = input.Clone()
    output.SetDirectory(0)
    minval = input.GetMinimum()
    NXbins = input.GetNbinsX()
    for i in range(1, NXbins + 1):
        output.SetBinContent(i, input.GetBinContent(i) - minval)
    return output

def print_matrix(input):
    for r in range(input.shape[1]):
        for c in range(input.shape[0]):
            print "%.4f " % input[r][c],
        print ""

def fill_histo(input):
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            input.SetBinContent(i, j, ((j - 1) * NYbins + i) / 25.0)

def histo_to_nparray(input):
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    output = np.zeros((NYbins, NXbins))
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            output[j - 1][i - 1] = input.GetBinContent(i, j)
    return output

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

def histo_to_nparray_1d(input):
    NXbins = input.GetNbinsX()
    output = np.zeros(NXbins)
    for i in range(1, NXbins + 1):
        output[i - 1] = input.GetBinContent(i)
    return output

def copy_nparray_to_histo_1d(input_array, output_histo):
    NXbins = output_histo.GetNbinsX()
    if NXbins != len(input_array):
        print "Error: Shape mismatch between array and histogram:"
        print "Array %d Histogram %d" % (len(input_array), NXbins)
        sys.exit()
    for i in range(1, NXbins + 1):
        output_histo.SetBinContent(i, input_array[i - 1])
    return output_histo

def print_histo(input):
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for j in range(1, NYbins + 1):
        for i in range(1, NXbins + 1):
            print "%.4f" % input.GetBinContent(i, j),
        print ""

def process_2d(original, root_filename, pdf_filename):
    NXbins = original.GetNbinsX()
    xaxis = original.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)

    NYbins = original.GetNbinsY()
    yaxis = original.GetYaxis()
    y_low = yaxis.GetBinLowEdge(1)
    y_high = yaxis.GetBinUpEdge(NYbins)
    bin_width = yaxis.GetBinWidth(1)

    temp = chisq_to_likelihood(original)
    original_array = histo_to_nparray(original)
    scale = 10.0
    resampled_array = resample(original_array, scale)
    resampled_histogram = ROOT.TH2D("atmos", "atmos", int(scale * original.GetNbinsX()), original.GetXaxis().GetXmin(), original.GetXaxis().GetXmax(), int(scale * original.GetNbinsY()), original.GetYaxis().GetXmin(), original.GetYaxis().GetXmax())
    copy_nparray_to_histo(resampled_array, resampled_histogram)
    likelihood_histogram = chisq_to_likelihood(resampled_histogram)
    likelihood_histogram.SetTitle("Resampled likelihood surface")
    likelihood_histogram.Draw("COLZ")
    canvas.Print(pdf_filename)

    smeared_histogram = smear(likelihood_histogram)
    smeared_histogram.SetTitle("Smeared likelihood surface")
    smeared_histogram.Draw("COLZ")
    canvas.Print(pdf_filename)

    chisq_histogram = likelihood_to_chisq(smeared_histogram)
    chisqmin = chisq_histogram.GetMinimum()
    chisq_histogram.SetContour(2)
    chisq_histogram.SetContourLevel(0, chisqmin + 2.3)
    chisq_histogram.SetContourLevel(1, chisqmin + 4.61)

    root_out_file = ROOT.TFile(root_filename, "RECREATE")
    root_out_file.cd()
    xaxis_label = xaxis.GetTitle()
    yaxis_label = yaxis.GetTitle()
    likelihood_histogram.GetXaxis().SetTitle(xaxis_label)
    likelihood_histogram.GetYaxis().SetTitle(yaxis_label)
    likelihood_histogram.Write("likelihood_2D")
    chisq_histogram.GetXaxis().SetTitle(xaxis_label)
    chisq_histogram.GetYaxis().SetTitle(yaxis_label)
    chisq_histogram.Write("chisq_map_2D")
    dchisq_histogram = make_delta_chisq(chisq_histogram)
    dchisq_histogram.GetXaxis().SetTitle(xaxis_label)
    dchisq_histogram.GetYaxis().SetTitle(yaxis_label)
    dchisq_histogram.Write("dchisq_map_2D")

    minbin = chisq_histogram.GetMinimumBin()
    binx = ROOT.Long(0)
    biny = ROOT.Long(0)
    dummy = ROOT.Long(0)
    chisq_histogram.GetBinXYZ(minbin, binx, biny, dummy);
  
    xmin = chisq_histogram.GetXaxis().GetBinCenter(binx);
    ymin = chisq_histogram.GetYaxis().GetBinCenter(biny);

    marker = ROOT.TMarker(xmin, ymin, ROOT.kFullStar);
    marker.SetMarkerColor(ROOT.kBlack);
    marker.SetMarkerSize(1.5);
    marker.Write("best_fit_point_marker")

    chisq_histogram.Draw("CONT LIST")
    canvas.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    list68 = contours.At(0)
    iter = ROOT.TIter(list68)
    opt = ""
    for i, g in enumerate(iter):
        name = "g68_valor_" + str(i)
        g.SetLineWidth(2)
        g.SetLineStyle(2)
        g.SetLineColor(ROOT.kRed)
        root_out_file.cd()
        g.Draw(opt + "CONT3")
        g.Write(name)
        opt = "SAME"

    list90 = contours.At(1)
    iter = ROOT.TIter(list90)
    opt = "SAME"
    for i, g in enumerate(iter):
        name = "g90_valor_" + str(i)
        g.SetLineWidth(2)
        g.SetLineStyle(1)
        g.SetLineColor(ROOT.kRed)
        root_out_file.cd()
        g.Draw(opt + "CONT3")
        g.Write(name)
        opt = "SAME"

    canvas.Print(pdf_filename + "]")
    root_out_file.Close()

def process_1d(original, root_filename, pdf_filename):
    NXbins = original.GetNbinsX()
    xaxis = original.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)

    original_dchisq_histogram = make_delta_chisq_1d(original)
    original_dchisq_histogram.SetTitle("Original #Delta#chi^2 surface")
    original_dchisq_histogram.Draw("")
    canvas.Print(pdf_filename)

    temp = chisq_to_likelihood_1d(original)
    original_array = histo_to_nparray_1d(original)
    scale = 10.0
    resampled_array = resample(original_array, scale)
    resampled_histogram = ROOT.TH1D("atmos", "atmos", int(scale * original.GetNbinsX()), original.GetXaxis().GetXmin(), original.GetXaxis().GetXmax())
    copy_nparray_to_histo_1d(resampled_array, resampled_histogram)
    likelihood_histogram = chisq_to_likelihood_1d(resampled_histogram)

    smeared_histogram = smear_1d(likelihood_histogram)

    chisq_histogram = likelihood_to_chisq_1d(smeared_histogram)
    chisqmin = chisq_histogram.GetMinimum()
    chisq_histogram.SetContour(2)
    chisq_histogram.SetContourLevel(0, chisqmin + 2.3)
    chisq_histogram.SetContourLevel(1, chisqmin + 4.61)

    root_out_file = ROOT.TFile(root_filename, "RECREATE")
    root_out_file.cd()
    xaxis_label = xaxis.GetTitle()
    likelihood_histogram.GetXaxis().SetTitle(xaxis_label)
    likelihood_histogram.Write("likelihood")
    chisq_histogram.GetXaxis().SetTitle(xaxis_label)
    chisq_histogram.Write("h_chisq")
    dchisq_histogram = make_delta_chisq_1d(chisq_histogram)
    dchisq_histogram.GetXaxis().SetTitle(xaxis_label)
    dchisq_histogram.Write("h_deltachisq")

    dchisq_histogram.SetTitle("Smeared #Delta#chi^2 surface")
    dchisq_histogram.Draw("")
    canvas.Print(pdf_filename)

    canvas.Print(pdf_filename + "]")
    root_out_file.Close()

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input contour file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output contour file.")
    parser.add_argument("-r", "--root-filename", dest="root_filename", default="", help="Output root file.")
    args = parser.parse_args()
    canvas = ROOT.TCanvas("c", "c", 700, 500)
    canvas.Print(args.out_filename + "[")

    in_file = ROOT.TFile(args.in_filename, "READ")
    if in_file.GetListOfKeys().Contains("chisq_map_2D"):
        original = in_file.Get("chisq_map_2D")
        process_2d(original, args.root_filename, args.out_filename)
    else:
        original = in_file.Get("h_chisq")
        process_1d(original, args.root_filename, args.out_filename)

    in_file.Close()

