#!/usr/bin/env python

import ROOT as R
import math
import array
import sys
import os

R.PyConfig.IgnoreCommandLineOptions = True

base_dir = os.environ["VALOR_EXPT"]
try:
  R.valor
except:
  R.gROOT.ProcessLine(".x %s/root/valor_expt.C"%base_dir);
  
R.gROOT.SetBatch(1)
  
valor = R.valor
util  = valor.util
conv  = valor.conv
opt   = valor.Opt.Instance()
env   = valor.Env.Instance()
rnd   = valor.Rnd.Instance()

arridx = util.ArrIdx.Instance()

treename = "valor_dune_nd_ntp"
stepsname = "mcsteps"

#______________________________________________________________________________________________________________________________

def main(input_fpl_file, fit_number, output_file, mcmc_files, burnin = 1000, make_plots = False):
  infile = R.TFile(input_fpl_file)
  assert(infile.IsOpen())
  intree = infile.Get(treename)
  assert(intree)
  
  intree.GetEntry(fit_number)
  
  fparm = intree.best_fit_params
  
  outfile = R.TFile(output_file,"RECREATE")
  out_tree = intree.CloneTree(0)
  
  npar = fparm.Size()
  
  n_entries_so_far = 0
  
  pca = R.TPrincipal(npar)
  
  var_mins = []
  var_maxs = []
  for ip in xrange(npar):
    var_mins.append( 9E9)
    var_maxs.append(-9E9)
  
  init_print_done = False
  
  for f in mcmc_files:
    mcfile = R.TFile(f)
    steptree = mcfile.Get(stepsname)
    n_entries = steptree.GetEntries()
    
    if (not init_print_done):
      print "Reading ",len(mcmc_files)," - this file has ",n_entries," entries."
      print "Minus burnin (",burnin,"), expecting ",len(mcmc_files)*(n_entries-burnin), "total steps."
      #~ init_print_done = True
    
    for iev in xrange(burnin, steptree.GetEntries()):
      steptree.GetEntry(iev)
      par = steptree.params
      
      n_entries_so_far += 1
      
      if n_entries_so_far % 1000 == 0:
        print "Reading step ",n_entries_so_far
      
      pca.AddRow(par)
      
      for ip in xrange(npar):
        var_mins[ip] = min(var_mins[ip],par[ip])
        var_maxs[ip] = max(var_maxs[ip],par[ip])
    #end event loop
    del steptree
    mcfile.Close()
    del mcfile
          
    print "Finished reading",f,"for",n_entries_so_far," steps (ignored",burnin,"for burn-in)"
  
  pca_mat = util.LowerTriangularMatrixToSymMatrix(pca.GetCovarianceMatrix())
  pca_vec = pca.GetMeanValues()
  
  # Cholesky decompose
  chol = util.CholeskyDecomposition(pca_mat)
  
  cholT = R.TMatrixD(R.TMatrix.kTransposed, chol)
  output_mat = util.MatrixToSymMatrix(chol * cholT)
  
  fparm.SetCovMatrix(output_mat)
  out_tree.Fill()
  outfile.cd()
  out_tree.Write()
  
  if make_plots:
    
    hists = []
    for ip in xrange(npar):
      parname = fparm.NameByIndex(ip)
      bf      = fparm.BestFitValueByIndex(ip)
      # keep best-fit centered
      max_diff = max(abs(var_mins[ip]-bf),abs(var_maxs[ip]-bf))
      h = R.TH1D(parname,parname,200,bf-max_diff,bf+max_diff)
      hists.append(h)
    
    print "\nReading them all again (ugh) for plotting"
    n_entries_so_far = 0
    init_print_done = False
    for f in mcmc_files:
      mcfile = R.TFile(f)
      steptree = mcfile.Get(stepsname)
      n_entries = steptree.GetEntries()
      
      if (not init_print_done):
        print "Reading ",len(mcmc_files)," - this file has ",n_entries," entries."
        print "Minus burnin (",burnin,"), expecting ",len(mcmc_files)*(n_entries-burnin), "total steps."
        #~ init_print_done = True
      
      for iev in xrange(burnin, steptree.GetEntries()):
        steptree.GetEntry(iev)
        par = steptree.params
      
        n_entries_so_far += 1
        
        if n_entries_so_far % 1000 == 0:
          print "Reading step for plots",n_entries_so_far
          
        par = steptree.params
        for ip in xrange(npar):
          hists[ip].Fill(par[ip])
      # End tree loop
      del steptree
      mcfile.Close()
      del mcfile
    # End file loop
        
    canv = R.TCanvas("c","c",1000,1000)
    canv.SetTopMargin  (0.1)
    canv.SetLeftMargin (0.1)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.05)
    
    plot_file = R.TFile("mcmc_hists.root","RECREATE")
    
    for ip,h in enumerate(hists):
      
      hmin,hmax = h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()
      gauss = R.TF1("mygaus_%d"%ip, "gaus(0)",hmin,hmax)
      gauss.SetParameter(0,h.GetMaximum())
      gauss.SetParameter(1,h.GetMean   ())
      gauss.SetParameter(2,h.GetRMS    ())
      gauss.SetParLimits(2,conv.kASmallNum,1)
      h.Fit(gauss,"RM0")
      gauss.SetLineColor(R.kRed)
      
      title_pave = R.TPaveText(0.3,0.9,0.7,1.0,"NDC")
      title_pave.AddText(h.GetTitle())
      
      fitted_pave = R.TPaveText(0.12,0.8,0.25,0.95,"NDC")
      fitted_pave.AddText("Fitted: ")
      fitted_pave.AddText("#mu"   +" = %.3lf"%gauss.GetParameter(1))
      fitted_pave.AddText("#sigma"+" = %.3lg"%gauss.GetParameter(2))
      fitted_pave.SetTextColor(R.kRed)
      cov_pave = R.TPaveText(0.75,0.8,0.93,0.95,"NDC")
      cov_pave.AddText("Covariance: ")
      cov_pave.AddText("#mu"   +" = %.3lf"%pca_vec[ip])
      cov_pave.AddText("#sigma"+" = %.3lg"%math.sqrt(output_mat[ip][ip]))
      format_pave(fitted_pave)
      format_pave(cov_pave)
      format_pave(title_pave)
      
      h.SetMinimum(0)
      h.Draw("hist")
      #~ func_hist.Draw("same")
      gauss.Draw("same")
      fitted_pave.Draw()
      cov_pave.Draw()
      title_pave.Draw()
      canv.Update()
      canv.Print("mcmc_steps_"+h.GetName()+".png")
      h.Write()
      canv.Write(h.GetName()+"_canv")
  # end if make plots
    
  pca_vec.Print()
  
  print "Done!"
  
def format_pave(pave):
  pave.SetTextFont(22)  
  pave.SetShadowColor(0)
  pave.SetFillColor(0)  
  #~ pave.SetFillStyle(0)  
  pave.SetLineColor(0)  
  pave.SetLineStyle(1)
  
#~ def tf1_to_th1d(func,hfit):
  #~ func.SetRange(hfit.GetXaxis().GetXmin(),hfit.GetXaxis().GetXmax())
  #~ func_hist = R.TH1D(hfit)
      #~ 
  #~ bin_width = func_hist.GetXaxis().GetBinWidth(1)
  #~ for i in xrange(func_hist.GetXaxis().GetNbins()):
    #~ bc = func_hist.GetXaxis().GetBinCenter(i+1)
    #~ val = func.Eval(bc)
    #~ func_hist.SetBinContent(i+1, val*bin_width)
  #~ func_hist.SetLineColor(R.kRed)
  #~ func_hist.SetLineWidth(2)
  #~ return func_hist

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-o', '--output', dest='output_file', default="combined_mcmc.root", help="Name for output FPL file.")
  parser.add_option('-i', '--input' , dest='input_file', default="fparm.root" , help="Input file result file")
  parser.add_option('-n', '--fit-number', dest='fit_number', default=0, type=int, help="Input tree entry")
  parser.add_option('-b', '--burn-in' , dest='burnin', default=1000, type=int, help="Number of steps to skip per chain")
  parser.add_option('-p', '--plots' , dest="make_plots", default=False, action="store_true", help="Make summary plots")
  
  o,other = parser.parse_args()
  
  main(o.input_file, o.fit_number, o.output_file, other, o.burnin, o.make_plots)
  
