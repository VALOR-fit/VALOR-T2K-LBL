#!/usr/bin/env python

import ROOT
import math
import array
import numpy as np
import sys
import os

env_var = "T2KOSCANA"
  
ROOT.gROOT.SetBatch(1)

n_bins_mu = 6
n_bins_e  = 12
first_bin_mu_in = n_bins_e
first_bin_mu_out = 0
first_bin_e_in = 0
first_bin_e_out = n_bins_mu

skip_row_col = [36,39,42,45]

nbins_per_beam_mode = n_bins_mu + n_bins_e
nbins_total = nbins_per_beam_mode*2

currents = dict( [ ("fhc","250kA" ) ,
                   ("rhc","m250kA") ])
#______________________________________________________________________________________________________________________________

def main(skdet_file, output_file="skdetfsi.root",escale=0.024):
  print "Reading", skdet_file 
  skdet_tf  = ROOT.TFile(skdet_file)

  fsi_only_matrix = skdet_tf.Get("FSI_SI_PN_Total")
  sk_only_matrix = skdet_tf.Get("SKJointError_Erec_Only")
  summed_matrix = skdet_tf.Get("SKJointError_Erec_Total")

  n = summed_matrix.GetNrows()
  print "Inputs"
  print "    FSI       SK    SKDETFSI"
  for i in range(0,n):
    print "%02d: %.8f   %.8f   %.8f" % (i, fsi_only_matrix(i,i), sk_only_matrix(i,i), summed_matrix(i,i))

  offset = 66 #Start of SKDet parameters
  out_tf = ROOT.TFile(output_file,"RECREATE")

  reordered_range = range(12,18)+range(0,12)+range(30,36)+range(18,30)+range(36,44)+range(44,45)
  summed_reordered = reorder_tmatrix(summed_matrix, reordered_range)
  sk_reordered = reorder_tmatrix(sk_only_matrix, reordered_range)
  fsi_reordered = reorder_tmatrix(fsi_only_matrix, reordered_range)

  escale_index = summed_reordered.GetNrows() - 1
  summed_reordered[escale_index][escale_index] = escale * escale
  summed_reordered.Write("skdetfsi")
  sk_reordered.Write("skdet")
  fsi_reordered.Write("fsi")

  print ""
  print "Outputs"
  print "    FSI       SK    SKDETFSI"
  for i in range(0,n):
    print "%02d: %.8f   %.8f   %.8f" % (i, fsi_reordered[i][i], sk_reordered[i][i], summed_reordered[i][i])


#  skdet_matrix_strip = skip_row_col_tmatrix(ko_matrix_reordered,skip_row_col)
#  skdet_matrix_strip.Write("covmat_skdetfsi_nozerobins_with_escale")

  
  out_tf.Close()

  skdet_tf  = ROOT.TFile(output_file)
  fsi_only_matrix = skdet_tf.Get("fsi")
  print ""
  print "Outputs after saving"
  print "    FSI       SK    SKDETFSI"
  for i in range(0,n):
    print "%02d: %.8f" % (i, fsi_only_matrix[i][i])

    
def reorder_tmatrix(matrix,new_order):
  n = matrix.GetNrows()
  assert(n==matrix.GetNcols())
  assert(n==len(new_order))
  new_matrix = ROOT.TMatrixD(n,n)
  for inew,iold  in enumerate(new_order):
    for jnew,jold in enumerate(new_order):
      new_matrix[inew][jnew] = matrix[iold][jold]
  return new_matrix

def truncate_tmatrix(matrix,offset):
  n_old = matrix.GetNrows()
  n = n_old - offset
  assert(n_old==matrix.GetNcols())
  new_matrix = ROOT.TMatrixD(n,n)
  for i in range(0,n):
    for j in range(0,n):
      new_matrix[i][j] = matrix[i+offset][j+offset]
  return new_matrix

def skip_row_col_tmatrix(matrix,skip):
  n_old = matrix.GetNrows()
  n = n_old - len(skip) 
  assert(n_old==matrix.GetNcols())
  new_matrix = ROOT.TMatrixD(n,n)
  i = 0
  inew = 0
  while i < n_old:
    j = 0
    jnew = 0
    if i in skip:
      i = i + 1
      continue
    while j < n_old:
      if j in skip:
        j = j + 1
        continue
      new_matrix[inew][jnew] = matrix[i][j]
      j = j + 1
      jnew = jnew + 1
    i = i + 1
    inew = inew + 1
  return new_matrix

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] skdet.root fsi_numu.root fsi_nue.root",description="Add up SKdet and FSI matrices")
  
  parser.add_option('-o', '--output_file', dest='output_file', default="skdetfsi.root", help="Current spline version.")
  parser.add_option('-e', '--escale', dest='escale', default=0.024, type=float, help="EScale systematic")
  
  o,other = parser.parse_args()
  
  main(other[0],o.output_file,o.escale)
  
  if other[3:]:
    print "Warning - ignored arguments: ",other[3:]
