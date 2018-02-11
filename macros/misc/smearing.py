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

def convolve(sub_sample, filter):
    return np.sum(sub_sample * filter)

def resample(input, scale):
    """Apply Catmull-Rom cubic interpolation"""
    if len(input.shape) == 1:
        # 1D resampling
        None
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

def smear_new(input, error=3.556768e-5):
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
    
    for i in range(1, NXbins + 1):
        x_smear = xaxis.GetBinCenter(i)
        for j in range(1, NYbins + 1):
            y_smear = yaxis.GetBinCenter(j)
            for k in range(1, NYbins + 1):
                likelihood = input.GetBinContent(i + 1, j + 1)
                gauss_bin = ROOT.TF1("gauss_bin", "gaus(0)", y_low, y_high)
                gauss_bin.SetParameters(norm, y_smear, error)
                if j != k:
                    y_temp = yaxis.GetBinCenter(k)
                    integral = gauss_bin.Integral(y_temp - 0.5 * bin_width, y_temp + 0.5 * bin_width)
                    likelihood = input.GetBinContent(i, k)
                    input.Fill(x_smear, y_smear, integral * likelihood)
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

def print_histo(input):
    NXbins = input.GetNbinsX()
    NYbins = input.GetNbinsY()
    for j in range(1, NYbins + 1):
        for i in range(1, NXbins + 1):
            print "%.4f" % input.GetBinContent(i, j),
        print ""


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input contour file.")
    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output contour file.")
    parser.add_argument("-r", "--root-filename", dest="root_filename", default="", help="Output root file.")
    args = parser.parse_args()
    canvas = ROOT.TCanvas("c", "c", 700, 500)
    canvas.Print(args.out_filename + "[")

    in_file = ROOT.TFile(args.in_filename, "READ")
    original = in_file.Get("chisq_map_2D")

    NXbins = original.GetNbinsX()
    xaxis = original.GetXaxis()
    x_low = xaxis.GetBinLowEdge(1)
    x_high = xaxis.GetBinUpEdge(NXbins)

    NYbins = original.GetNbinsY()
    yaxis = original.GetYaxis()
    y_low = yaxis.GetBinLowEdge(1)
    y_high = yaxis.GetBinUpEdge(NYbins)
    bin_width = yaxis.GetBinWidth(1)

    min = float(np.inf)
    for i in range(1, NXbins + 1):
        for j in range(1, NYbins + 1):
            val = original.GetBinContent(i, j)
            if val < min:
                min = val
                min_x = i
                min_y = j

    left_bin = 0
    right_bin = 0
    i = min_x
    for j in range(1, NYbins + 1):
        if left_bin == 0:
            if (original.GetBinContent(i, j) - min) < 4.61:
                left_bin = j
        elif right_bin == 0:
            if (original.GetBinContent(i, j) - min) > 4.61:
                right_bin = j
    print "2 sigma interval [", original.GetYaxis().GetBinCenter(left_bin), original.GetYaxis().GetBinCenter(right_bin), "]"

    temp = chisq_to_likelihood(original)
    temp.SetTitle("Original likelihood surface")
    temp.Draw("COLZ")
    canvas.Print(args.out_filename)
    original_array = histo_to_nparray(original)
    scale = 10.0
    resampled_array = resample(original_array, scale)
    resampled_histogram = ROOT.TH2D("atmos", "atmos", int(scale * original.GetNbinsX()), original.GetXaxis().GetXmin(), original.GetXaxis().GetXmax(), int(scale * original.GetNbinsY()), original.GetYaxis().GetXmin(), original.GetYaxis().GetXmax())
    copy_nparray_to_histo(resampled_array, resampled_histogram)
    likelihood_histogram = chisq_to_likelihood(resampled_histogram)
    likelihood_histogram.SetTitle("Resampled likelihood surface")
    likelihood_histogram.Draw("COLZ")
    canvas.Print(args.out_filename)

    smeared_histogram = smear(likelihood_histogram)
    smeared_histogram.SetTitle("Smeared likelihood surface")
    smeared_histogram.Draw("COLZ")
    canvas.Print(args.out_filename)

    chisq_histogram = likelihood_to_chisq(smeared_histogram)
    chisqmin = chisq_histogram.GetMinimum()
    chisq_histogram.SetContour(2)
    chisq_histogram.SetContourLevel(0, chisqmin + 2.3)
    chisq_histogram.SetContourLevel(1, chisqmin + 4.61)

    dmsq_original = marginalise(temp)
    dmsq_original.SetTitle("Original likelihood")
    dmsq_original.Draw()
    canvas.Print(args.out_filename)
    dmsq_smeared = marginalise(smeared_histogram)
    dmsq_smeared.SetTitle("Smeared likelihood")
    dmsq_smeared.Draw()
    canvas.Print(args.out_filename)

    theta23_original = marginalise(temp, "y")
    theta23_original.SetTitle("Original likelihood")
    theta23_original.Draw()
    canvas.Print(args.out_filename)
    theta23_smeared = marginalise(smeared_histogram, "y")
    theta23_smeared.SetTitle("Smeared likelihood")
    theta23_smeared.Draw()
    canvas.Print(args.out_filename)

    root_out_file = ROOT.TFile(args.root_filename, "RECREATE")
    root_out_file.cd()
    likelihood_histogram.Write("likelihood_2D")
    chisq_histogram.Write("chisq_map_2D")
    dchisq_histogram = make_delta_chisq(chisq_histogram)
    dchisq_histogram.Write("dchisq_map_2D")
    contours = array.array('d', [ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 ])
    dchisq_histogram.SetContour(len(contours), contours)
    dchisq_histogram.Draw("COLZ")

    dmsq_original = likelihood_to_chisq(dmsq_original)
    dmsq_original = make_delta_chisq(dmsq_original)
    dmsq_original.SetTitle("Original #Delta#chi^{2} surface")
    dmsq_original.Draw()
    canvas.Print(args.out_filename)

    theta23_original = likelihood_to_chisq(theta23_original)
    theta23_original = make_delta_chisq(theta23_original)
    theta23_original.SetTitle("Original #Delta#chi^{2} surface")
    theta23_original.Draw()
    canvas.Print(args.out_filename)

    dmsq_smeared = likelihood_to_chisq(dmsq_smeared)
    dmsq_smeared = make_delta_chisq(dmsq_smeared)
    dmsq_smeared.SetTitle("Smeared #Delta#chi^{2} surface")
    dmsq_smeared.Draw()
    canvas.Print(args.out_filename)

    theta23_smeared = likelihood_to_chisq(theta23_smeared)
    theta23_smeared = make_delta_chisq(theta23_smeared)
    theta23_smeared.SetTitle("Smeared #Delta#chi^{2} surface")
    theta23_smeared.Draw()
    canvas.Print(args.out_filename)

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

    canvas.Print(args.out_filename + "]")
    root_out_file.Close()
    in_file.Close()


#    print "Original:"
#    print_histo(original)
#    print "Resampled:"
#    print_histo(resampled_histogram)
#    print "Likelihood:"
#    print_histo(likelihood_histogram)
#    print "Smeared:"
#    print_histo(smeared_histogram)
#    print "Chisq:"
#    print_histo(chisq_histogram)

