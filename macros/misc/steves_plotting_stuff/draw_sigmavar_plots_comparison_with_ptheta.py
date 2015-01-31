#!/usr/bin/python
import sys
import ROOT 
import os

env_var = "T2KOSCANA"

pthet_filename = '/home/steve/analysis/systematic_variations/ptheta_var.root'

valor_filename = '/home/steve/analysis/systematic_variations/systematic_variations_RHC1Rmu.root'

pthet_pattern     = "[PARAM]_[VAR]"
valor_pattern     = "Ratio_RHC1Rmu_[PARAM]_[VAR]"

# paired things have valor, then pthet
variations = [ ("p1","3"),
               ("p3","4"),
               ("m1","2"),
               ("m3","1") ]

parameters = [( 'f_banffmaqe_rhc'          ,  'Par_25'  ) ,
              ( 'f_banffpfermi_rhc'        ,  'Par_26'  ) ,
              ( 'f_banffmec_rhc'           ,  'Par_27'  ) ,
              ( 'f_banffbindE_rhc'         ,  'Par_28'  ) ,
              ( 'f_banffca5_rhc'           ,  'Par_29'  ) ,
              ( 'f_banffmares_rhc'         ,  'Par_30'  ) ,
              ( 'f_banffbgres_rhc'         ,  'Par_31'  ) ,
              ( 'f_banffnorm_nuetonumu_rhc',  'Par_32' ) ,
              ( 'f_banffshapeCCoth_rhc'    ,  'Par_33' ) ,
              ( 'f_banffnorm_cccoh_rhc'    ,  'Par_34' ) ,
              ( 'f_banffnorm_nccoh_rhc'    ,  'Par_35' ) ,
              ( 'f_banffnorm_ncoth_rhc'    ,  'Par_36' ) ,
              ( 'f_banff00_rhc'            ,  'Par_0' ) , 
              ( 'f_banff01_rhc'            ,  'Par_1' ) , 
              ( 'f_banff02_rhc'            ,  'Par_2' ) , 
              ( 'f_banff03_rhc'            ,  'Par_3' ) , 
              ( 'f_banff04_rhc'            ,  'Par_4' ) , 
              ( 'f_banff05_rhc'            ,  'Par_5' ) , 
              ( 'f_banff06_rhc'            ,  'Par_6' ) , 
              ( 'f_banff07_rhc'            ,  'Par_7' ) , 
              ( 'f_banff08_rhc'            ,  'Par_8' ) , 
              ( 'f_banff09_rhc'            ,  'Par_9' ) , 
              ( 'f_banff10_rhc'            ,  'Par_10' ) , 
              ( 'f_banff11_rhc'            ,  'Par_11' ) , 
              ( 'f_banff12_rhc'            ,  'Par_12' ) , 
              ( 'f_banff13_rhc'            ,  'Par_13' ) , 
              ( 'f_banff14_rhc'            ,  'Par_14' ) , 
              ( 'f_banff15_rhc'            ,  'Par_15' ) , 
              ( 'f_banff16_rhc'            ,  'Par_16' ) , 
              ( 'f_banff17_rhc'            ,  'Par_17' ) , 
              ( 'f_banff18_rhc'            ,  'Par_18' ) , 
              ( 'f_banff19_rhc'            ,  'Par_19' ) , 
              ( 'f_banff20_rhc'            ,  'Par_20' ) , 
              ( 'f_banff21_rhc'            ,  'Par_21' ) , 
              ( 'f_banff22_rhc'            ,  'Par_22' ) , 
              ( 'f_banff23_rhc'            ,  'Par_23' ) , 
              ( 'f_banff24_rhc'            ,  'Par_24' ) , 
              ( 'f_skdetfsi000_rhc'        ,  'Par_50') , 
              ( 'f_skdetfsi001_rhc'        ,  'Par_51') , 
              ( 'f_skdetfsi002_rhc'        ,  'Par_52') , 
              ( 'f_skdetfsi003_rhc'        ,  'Par_53') , 
              ( 'f_skdetfsi004_rhc'        ,  'Par_54') , 
              ( 'f_skdetfsi005_rhc'        ,  'Par_55') ,
              #~ ( 'f_sk_e_scale_rhc' ,  'skdet_par_36')
              ]
              
      
line_colors = dict()
line_colors["m1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["m3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def main():
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)
  
  pthet_file = ROOT.TFile(pthet_filename)
  valor_file = ROOT.TFile(valor_filename)
  
  assert(pthet_file)
  assert(valor_file)
    
  valor_hists             = dict() # only valor names will be used as keys
  pthet_hists             = dict() # only valor names will be used as keys
  valor_hists_binno       = dict() # only valor names will be used as keys
  pthet_hists_binno       = dict() # only valor names will be used as keys
  canvases                = dict()
  
  for iparam,(valor_param,pthet_param) in enumerate(parameters):
    
    valor_hists             [valor_param] = dict()
    pthet_hists             [valor_param] = dict()
    valor_hists_binno       [valor_param] = dict()
    pthet_hists_binno [valor_param] = dict()
    
    for valor_var,pthet_var in variations:
      
      valor_histname = valor_pattern
      for old,new in [("[PARAM]", valor_param),
                      ("[VAR]"  , valor_var  )] :
        valor_histname = valor_histname.replace(old,new)
                                              
      pthet_histname = pthet_pattern
      for old,new in [("[PARAM]"    , pthet_param    ),
                      ("[VAR]"      , pthet_var      ) ]:
        pthet_histname = pthet_histname.replace(old,new)
                                             
      valor_hists[valor_param][valor_var] = valor_file.Get(valor_histname)
      
      if not valor_hists[valor_param][valor_var]:
        print "Failed to get VALOR histogram ",valor_histname
        
      pthet_hists[valor_param][valor_var] = pthet_file.Get(pthet_histname)
      if not pthet_hists[valor_param][valor_var]:
        print "Failed to get MaCh3 histogram ",pthet_histname
        
      valor_hists[valor_param][valor_var].SetLineWidth(1)
      valor_hists[valor_param][valor_var].SetLineColor(line_colors[valor_var])
        
      pthet_hists[valor_param][valor_var].SetName(valor_hists[valor_param][valor_var].GetName()+"_pthet")
      pthet_hists[valor_param][valor_var].SetLineStyle(2)
      pthet_hists[valor_param][valor_var].SetLineWidth(2)
      pthet_hists[valor_param][valor_var].SetLineColor(valor_hists[valor_param][valor_var].GetLineColor()+2)
      
      valor_hists_binno[valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_param][valor_var])
      pthet_hists_binno[valor_param][valor_var] = hist_to_bin_number(pthet_hists[valor_param][valor_var])
      
    hist_limits = []
    for valor_var,pthet_var in variations:
      for bin_content in (hist_to_list(valor_hists[valor_param][valor_var])+hist_to_list(pthet_hists[valor_param][valor_var])):
        if abs(bin_content) > 1E-5:
          hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    pdf_printing = [""]*len(parameters)
    pdf_printing[ 0] = "("
    pdf_printing[-1] = ")"
    
    canvases[valor_param] = ROOT.TCanvas(valor_param,valor_param,1280,720)
    option=""
    for valor_var,pthet_var in variations:
      valor_hists[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+pthet_param)
      valor_hists[valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
      valor_hists[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_hists[valor_param][valor_var].Draw(option)
    
    canvases[valor_param].Update()
    canvases[valor_param].Print("syst_var_valor_vs_ptheta_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[valor_param+"_binno"] = ROOT.TCanvas(valor_param+"_binno",valor_param+"_binno",1280,720)
    option=""
    for valor_var,pthet_var in variations:
      valor_hists_binno[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+pthet_param)
      valor_hists_binno[valor_param][valor_var].GetXaxis().SetTitle("Bin number")
      valor_hists_binno[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists_binno[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists_binno[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists_binno[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_hists_binno[valor_param][valor_var].Draw(option)
    
    canvases[valor_param+"_binno"].Update()
    canvases[valor_param+"_binno"].Print("syst_var_valor_vs_ptheta_binno"+".pdf"+pdf_printing[iparam])
    

def hist_to_bin_number(h):
  nbins = h.GetXaxis().GetNbins()
  hnew = ROOT.TH1D(h.GetName()+"_binno",h.GetTitle()+" by bins",nbins,0,nbins)
  
  hnew.SetLineColor(h.GetLineColor())
  hnew.SetLineStyle(h.GetLineStyle())
  hnew.SetLineWidth(h.GetLineWidth())
  
  for ibin in xrange(1,nbins+1):
    hnew.Fill(ibin-0.5, h.GetBinContent(ibin))
    
  hnew.GetXaxis().SetTitle("Bin number")
  
  return hnew
  
def ratio_hist(hvar,hnom):
  nbins = hvar.GetXaxis().GetNbins()
  hratio = ROOT.TH1D(hvar)
  hratio.SetName (hvar.GetName ()+"_ratio")
  hratio.SetTitle(hvar.GetTitle()+" Ratio")
  hratio.GetYaxis().SetTitle("Ratio")
  hratio.Divide(hnom)
  return hratio
  

def hist_to_list(h):
  nbins = h.GetXaxis().GetNbins()
  return [h.GetBinContent(ibin) for ibin in xrange(1,nbins+1)]
  
if __name__=="__main__":
  main()

