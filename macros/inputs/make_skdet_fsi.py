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

nbins_per_beam_mode = n_bins_mu + n_bins_e
nbins_total = nbins_per_beam_mode*2

currents = dict( [ ("fhc","250kA" ) ,
                   ("rhc","m250kA") ])
#______________________________________________________________________________________________________________________________

def main(skdet_file, fsi_numu_file, fsi_nue_file, output_file="skdetfsi.root",escale=0.024):
  
  fsi_mu_tf = ROOT.TFile(fsi_numu_file)
  fsi_e_tf  = ROOT.TFile(fsi_nue_file)
  
  assert(fsi_mu_tf )
  assert(fsi_e_tf  )
  
  skdet_tf  = ROOT.TFile(skdet_file)
  skdet_matrix = skdet_tf.Get("postfit_cov") 
  #skdet_matrix = skdet_tf.Get("SKJointErrorMatrix_Out") 
  #skdet_matrix = skdet_tf.Get("Covariance_Matrix_Erec")
  assert(skdet_matrix.GetNrows()==nbins_total)
  
  out_tf = ROOT.TFile(output_file,"RECREATE")
  
  new_skdet_order_fhc = range(first_bin_mu_in, first_bin_mu_in+n_bins_mu)+range(first_bin_e_in, first_bin_e_in+n_bins_e)
  new_skdet_order_all = new_skdet_order_fhc + [n+nbins_per_beam_mode for n in new_skdet_order_fhc]
  
  big_skdet_matrix  = reorder_tmatrix(skdet_matrix, new_skdet_order_all)
  
  big_fsi_matrix    = ROOT.TMatrixD(nbins_total,nbins_total)
  big_pn_matrix     = ROOT.TMatrixD(nbins_total,nbins_total)
  big_fsi_matrix.Zero()
  big_pn_matrix .Zero()
  
  e_indices_out  = []
  mu_indices_out = []
  
  fsi_mu = dict()
  fsi_e  = dict()
  pn_mu  = dict()
  pn_e   = dict()
  
  # Put together the different submatrices
  for i_mode,beam_mode in enumerate(["fhc","rhc"]):
    
    fsi_mu[beam_mode] = make_tmatrix_from_th2(fsi_mu_tf.Get("FSIandSICovMat_"+currents[beam_mode]))
    fsi_e [beam_mode] = make_tmatrix_from_th2(fsi_e_tf .Get("FSIandSICovMat_"+currents[beam_mode]))
    
    pn_e [beam_mode] = make_tmatrix_from_th2(fsi_e_tf .Get("PNCovMat_"+currents[beam_mode]))
    
    pn_mu[beam_mode] = ROOT.TMatrixD(n_bins_mu,n_bins_mu)
    pn_mu[beam_mode].Zero()
    
    assert(fsi_mu[beam_mode].GetNrows()==n_bins_mu)
    assert(fsi_mu[beam_mode].GetNcols()==n_bins_mu)
    assert(fsi_e [beam_mode].GetNrows()==n_bins_e )
    assert(fsi_e [beam_mode].GetNcols()==n_bins_e )
    assert(pn_mu [beam_mode].GetNrows()==n_bins_mu)
    assert(pn_mu [beam_mode].GetNcols()==n_bins_mu)
    assert(pn_e  [beam_mode].GetNrows()==n_bins_e )
    assert(pn_e  [beam_mode].GetNcols()==n_bins_e )
    
    first_mu_row  = i_mode * nbins_per_beam_mode
    first_e_row   = first_mu_row + n_bins_mu
    last_mu_row   = first_e_row - 1
    last_e_row    = last_mu_row + n_bins_e
    
    big_fsi_matrix.SetSub(first_mu_row, first_mu_row, fsi_mu[beam_mode])
    big_fsi_matrix.SetSub(first_e_row , first_e_row , fsi_e [beam_mode])
    big_pn_matrix .SetSub(first_mu_row, first_mu_row, pn_mu [beam_mode])
    big_pn_matrix .SetSub(first_e_row , first_e_row , pn_e  [beam_mode])
    
    
  big_fsi_pn_matrix = big_fsi_matrix + big_pn_matrix
  big_skdet_fsi_matrix = big_fsi_pn_matrix + big_skdet_matrix
  
  mu_indices = [] 
  e_indices  = [] 
  
  # Write single beam mode fit matrices
  for i_mode,beam_mode in enumerate(["fhc","rhc"]):
    
    # Joint
    first_row = i_mode*nbins_per_beam_mode
    last_row  = first_row+nbins_per_beam_mode-1
    big_fsi_matrix      .GetSub(first_row,last_row,first_row,last_row).Write("covmat_fsi_only_"   +beam_mode)
    big_fsi_pn_matrix   .GetSub(first_row,last_row,first_row,last_row).Write("covmat_fsi_pn_only_"+beam_mode)
    big_pn_matrix       .GetSub(first_row,last_row,first_row,last_row).Write("covmat_pn_only_"    +beam_mode)
    big_skdet_matrix    .GetSub(first_row,last_row,first_row,last_row).Write("covmat_skdet_only_" +beam_mode)
    big_skdet_fsi_matrix.GetSub(first_row,last_row,first_row,last_row).Write("covmat_skdetfsi_"   +beam_mode) 
    
    # Numu
    first_row_mu = first_row+first_bin_mu_out
    last_row_mu  = first_row_mu+n_bins_mu-1
    big_fsi_matrix      .GetSub(first_row_mu,last_row_mu,first_row_mu,last_row_mu).Write("covmat_fsi_only_"   +beam_mode+"_numu")
    big_fsi_pn_matrix   .GetSub(first_row_mu,last_row_mu,first_row_mu,last_row_mu).Write("covmat_fsi_pn_only_"+beam_mode+"_numu")
    big_pn_matrix       .GetSub(first_row_mu,last_row_mu,first_row_mu,last_row_mu).Write("covmat_pn_only_"    +beam_mode+"_numu")
    big_skdet_matrix    .GetSub(first_row_mu,last_row_mu,first_row_mu,last_row_mu).Write("covmat_skdet_only_" +beam_mode+"_numu")
    big_skdet_fsi_matrix.GetSub(first_row_mu,last_row_mu,first_row_mu,last_row_mu).Write("covmat_skdetfsi_"   +beam_mode+"_numu")
    
    # NuE
    first_row_e = first_row+first_bin_e_out
    last_row_e  = first_row_e+n_bins_e-1
    big_fsi_matrix      .GetSub(first_row_e,last_row_e,first_row_e,last_row_e).Write("covmat_fsi_only_"   +beam_mode+"_nue")
    big_fsi_pn_matrix   .GetSub(first_row_e,last_row_e,first_row_e,last_row_e).Write("covmat_fsi_pn_only_"+beam_mode+"_nue")
    big_pn_matrix       .GetSub(first_row_e,last_row_e,first_row_e,last_row_e).Write("covmat_pn_only_"    +beam_mode+"_nue")
    big_skdet_matrix    .GetSub(first_row_e,last_row_e,first_row_e,last_row_e).Write("covmat_skdet_only_" +beam_mode+"_nue")
    big_skdet_fsi_matrix.GetSub(first_row_e,last_row_e,first_row_e,last_row_e).Write("covmat_skdetfsi_"   +beam_mode+"_nue") 
    
    # store the indices so we can retrieve the right submatrices for both mode numu or both mode nue
    mu_indices.extend(range(first_row_mu,last_row_mu +1))
    e_indices .extend(range(first_row_e ,last_row_e  +1))
    
  
  # And the RHC/FHC matrices for each set of analyses
  fsi_bothmodes_mu       = ROOT.TMatrixD(2*n_bins_mu,2*n_bins_mu)
  fsi_pn_bothmodes_mu    = ROOT.TMatrixD(2*n_bins_mu,2*n_bins_mu)
  pn_bothmodes_mu        = ROOT.TMatrixD(2*n_bins_mu,2*n_bins_mu)
  skdet_bothmodes_mu     = ROOT.TMatrixD(2*n_bins_mu,2*n_bins_mu)
  skdet_fsi_bothmodes_mu = ROOT.TMatrixD(2*n_bins_mu,2*n_bins_mu)
  
  fsi_bothmodes_e       = ROOT.TMatrixD(2*n_bins_e,2*n_bins_e)
  fsi_pn_bothmodes_e    = ROOT.TMatrixD(2*n_bins_e,2*n_bins_e)
  pn_bothmodes_e        = ROOT.TMatrixD(2*n_bins_e,2*n_bins_e)
  skdet_bothmodes_e     = ROOT.TMatrixD(2*n_bins_e,2*n_bins_e)
  skdet_fsi_bothmodes_e = ROOT.TMatrixD(2*n_bins_e,2*n_bins_e)
  
  for newx,oldx in enumerate(mu_indices):
    for newy,oldy in enumerate(mu_indices):
      fsi_bothmodes_mu       [newx][newy] = big_fsi_matrix      [oldx][oldy]
      fsi_pn_bothmodes_mu    [newx][newy] = big_fsi_pn_matrix   [oldx][oldy]
      pn_bothmodes_mu        [newx][newy] = big_pn_matrix       [oldx][oldy]
      skdet_bothmodes_mu     [newx][newy] = big_skdet_matrix    [oldx][oldy]
      skdet_fsi_bothmodes_mu [newx][newy] = big_skdet_fsi_matrix[oldx][oldy]
      
  for newx,oldx in enumerate(e_indices):
    for newy,oldy in enumerate(e_indices):
      fsi_bothmodes_e       [newx][newy] = big_fsi_matrix      [oldx][oldy]
      fsi_pn_bothmodes_e    [newx][newy] = big_fsi_pn_matrix   [oldx][oldy]
      pn_bothmodes_e        [newx][newy] = big_pn_matrix       [oldx][oldy]
      skdet_bothmodes_e     [newx][newy] = big_skdet_matrix    [oldx][oldy]
      skdet_fsi_bothmodes_e [newx][newy] = big_skdet_fsi_matrix[oldx][oldy]
  
  # Write full-size matrices
  big_fsi_matrix      .Write("covmat_fsi_only_full")
  big_fsi_pn_matrix   .Write("covmat_fsi_pn_only_full")
  big_pn_matrix       .Write("covmat_pn_only_full")
  big_skdet_matrix    .Write("covmat_skdet_only_full")
  big_skdet_fsi_matrix.Write("covmat_skdetfsi_full")
  
  ## Write broken down matrices
  fsi_bothmodes_mu.Write("covmat_fsi_only_numu")
  fsi_bothmodes_e .Write("covmat_fsi_only_nue")
  
  fsi_pn_bothmodes_mu.Write("covmat_fsi_pn_only_numu")
  fsi_pn_bothmodes_e .Write("covmat_fsi_pn_only_nue")
  
  pn_bothmodes_mu.Write("covmat_pn_only_numu")
  pn_bothmodes_e .Write("covmat_pn_only_nue")
  
  skdet_bothmodes_mu.Write("covmat_skdet_only_numu")
  skdet_bothmodes_e .Write("covmat_skdet_only_nue")
  
  skdet_fsi_bothmodes_mu.Write("covmat_skdetfsi_numu")
  skdet_fsi_bothmodes_e .Write("covmat_skdetfsi_nue")
  
  # Make big matrices with extra row for EScale
  escale_elem = escale*escale
  
  add_diag_row_to_tmatrix(big_skdet_matrix      , escale_elem).Write("covmat_skdet_only_full"    +"_with_escale")
  add_diag_row_to_tmatrix(big_skdet_fsi_matrix  , escale_elem).Write("covmat_skdetfsi_full"      +"_with_escale")
  add_diag_row_to_tmatrix(skdet_bothmodes_mu    , escale_elem).Write("covmat_skdet_only_numu"    +"_with_escale")
  add_diag_row_to_tmatrix(skdet_bothmodes_e     , escale_elem).Write("covmat_skdet_only_nue"     +"_with_escale")
  add_diag_row_to_tmatrix(skdet_fsi_bothmodes_mu, escale_elem).Write("covmat_skdetfsi_numu"      +"_with_escale")
  add_diag_row_to_tmatrix(skdet_fsi_bothmodes_e , escale_elem).Write("covmat_skdetfsi_nue"       +"_with_escale")
  
  out_tf.Close()

def make_tmatrix_from_th2(hist):
  nrows = hist.GetXaxis().GetNbins()
  ncols = hist.GetYaxis().GetNbins()
  matr = ROOT.TMatrixD(nrows,ncols)
  
  for i in xrange(nrows):
    for j in xrange(ncols):
      matr[i][j] = hist.GetBinContent(i+1,j+1)
  
  return matr
  
def reorder_tmatrix(matrix,new_order):
  n = matrix.GetNrows()
  assert(n==matrix.GetNcols())
  assert(n==len(new_order))
  new_matrix = ROOT.TMatrixD(n,n)
  for inew,iold  in enumerate(new_order):
    for jnew,jold in enumerate(new_order):
      new_matrix[inew][jnew] = matrix[iold][jold]
  return new_matrix

def add_diag_row_to_tmatrix(matrix,new_element):
  new_matrix = ROOT.TMatrixD(matrix.GetNrows()+1,matrix.GetNcols()+1).Zero()
  new_matrix.SetSub(0,0,matrix)
  new_matrix[matrix.GetNrows()][matrix.GetNcols()] = new_element
  return new_matrix
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] skdet.root fsi_numu.root fsi_nue.root",description="Add up SKdet and FSI matrices")
  
  parser.add_option('-o', '--output_file', dest='output_file', default="skdetfsi.root", help="Current spline version.")
  parser.add_option('-e', '--escale', dest='escale', default=0.024, type=float, help="EScale systematic")
  
  o,other = parser.parse_args()
  
  main(other[0],other[1],other[2],o.output_file,o.escale)
  
  if other[3:]:
    print "Warning - ignored arguments: ",other[3:]
