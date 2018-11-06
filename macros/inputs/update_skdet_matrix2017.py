#!/usr/bin/env python

import ROOT
import math
import array
import numpy as np
import sys
import os
import argparse

env_var = "T2KOSCANA"
  
ROOT.gROOT.SetBatch(1)

def add_escale(matrix, escale=0.0):
    n = matrix.GetNrows()
    new_matrix = ROOT.TMatrixD(n + 1, n + 1)
    for r in range(n):
        for c in range(n):
            new_matrix[r][c] = matrix(r,c)
    print escale
    new_matrix[n][n] = escale**2
    return new_matrix

def load_sk_total_matrix(filename):
    """ Loads an SK matrix with detector and FSI+SI+PN (or some subset) components """
    sk_file = ROOT.TFile(filename)
    sk_only_matrix = ROOT.TMatrixD(sk_file.Get("SKJointError_Erec_Only"))
    fsi_only_matrix = ROOT.TMatrixD(sk_file.Get("FSI_SI_PN_Total"))
    summed_matrix = ROOT.TMatrixD(sk_file.Get("SKJointError_Erec_Total"))
    sk_file.Close()

    return sk_only_matrix, fsi_only_matrix, summed_matrix

def load_sk_detector_matrix(filename):
    """ Loads an SK matrix with only detector components """
    sk_file = ROOT.TFile(filename)
    sk_only_matrix = ROOT.TMatrixD(sk_file.Get("Covariance_Matrix_Erec"))
    sk_file.Close()

    return sk_only_matrix

def make_sk_total_matrix(filename, det_matrix, fsi_matrix):
    n = det_matrix.GetNrows()
    sk_matrix = ROOT.TMatrixD(n, n)
    print "Det:"
    print_matrix(det_matrix)
    print "FSI:"
    print_matrix(fsi_matrix)
    for r in range(n):
        for c in range(n):
            sk_matrix[r][c] = det_matrix(r,c) + fsi_matrix(r,c)
    print "SK:"
    print_matrix(sk_matrix)
    sk_file = ROOT.TFile(filename, "RECREATE")
    sk_matrix.Write("SKJointError_Erec_Total")
    det_matrix.Write("SKJointError_Erec_Only")
    fsi_matrix.Write("FSI_SI_PN_Total")
    sk_file.Close()

def print_matrix(matrix):
    n = matrix.GetNrows()
    for i in range(n):
        print "%02d: det: %.8f" % (i, matrix(i, i))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Index splines")
    parser.add_argument("-d", "--sk-detector-file", dest="sk_det_filename", default="", help="Input file containing SK detector matrix.")
    parser.add_argument("-t", "--sk-total-file", dest="sk_tot_filename", default="", help="Input file containing SK detector matrix with FSI+SI+PN components.")
    parser.add_argument("-f", "--filename", dest="out_filename", default="", help="Output filename for the full SK matrix.")

    args = parser.parse_args()
    _, fsi_matrix, _ = load_sk_total_matrix(args.sk_tot_filename)
    det_matrix = load_sk_detector_matrix(args.sk_det_filename)
    # Old matrix already had escale in it, detector only doesn't, so just add a zero component
    det_matrix = add_escale(det_matrix)
    make_sk_total_matrix(args.out_filename, det_matrix, fsi_matrix)
  
