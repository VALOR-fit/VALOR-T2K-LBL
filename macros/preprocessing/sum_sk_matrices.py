#!/usr/bin/env python

import ROOT
import math
import array
import numpy as np
import sys
import os

env_var = "T2KOSCANA"
  
ROOT.gROOT.SetBatch(1)

def main(skdetfsi_filename, skdet_filename, output_file="skdetfsi.root",escale=0.024):
  
  skdetfsi_file  = ROOT.TFile(skdetfsi_filename)
  skdet_file = ROOT.TFile(skdet_filename)

  fsi_matrix = skdetfsi_file.Get("FSI_SI_PN_Total")
  skdet_matrix = skdet_file.Get("Covariance_Matrix_Erec")

  n = fsi_matrix.GetNrows()
  skdet_new_matrix = ROOT.TMatrixD(n, n)
  summed_matrix = ROOT.TMatrixD(n, n)

  for r in range(n - 1):
    for c in range(n - 1):
      skdet_new_matrix[r][c] = skdet_matrix[r][c]

  for r in range(n):
    for c in range(n):
      summed_matrix[r][c] = skdet_new_matrix[r][c] + fsi_matrix[r][c]
  
  out_tf = ROOT.TFile(output_file,"RECREATE")

  summed_matrix.Write("SKJointError_Erec_Total")
  skdet_new_matrix.Write("SKJointError_Erec_Only")
  fsi_matrix.Write("FSI_SI_PN_Total")

  out_tf.Close()
    
if __name__=="__main__":
  main(sys.argv[1], sys.argv[2])
  
