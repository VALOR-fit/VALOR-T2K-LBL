#!/usr/bin/env python

import argparse
import ROOT
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
        elif y >= input.shape[1]:
            y = input.shape[1] - 1
        for dx in range(-1, 3):
            x = col + dx
            if x < 0:
                x = 0
            elif x >= input.shape[0]:
                x = input.shape[0] - 1
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

def print_matrix(input):
    for r in range(input.shape[1]):
        for c in range(input.shape[0]):
            print "%.4f " % input[r][c],
        print ""


if __name__=="__main__":
#    parser = argparse.ArgumentParser(description="Index splines")
#    parser.add_argument("-f", "--filename", dest="in_filename", default="", help="Input p-theta fake data file.")
#    parser.add_argument("-o", "--output-filename", dest="out_filename", default="", help="Output Valor fake data file.")
#    args = parser.parse_args()
    print ""
    image = np.reshape(np.array([range(1,26)]), (5, 5))
#    image = np.ones((5, 5))
    print_matrix(image)
    print ""
    output = resample(image, 3.0)
    print_matrix(output)
