#!/usr/bin/env python
import sys
import ROOT 
import os

env_var = "T2KOSCANA"

mach3_filename = '/home/steve/analysis/systematic_variations/sigmavar2014SK.root'
pthet_filename = '/home/steve/analysis/systematic_variations/systvariation.root'
valor_filename = '/home/steve/analysis/systematic_variations/systematic_variations_RHC1Rmu.root'

mach3_nom_pattern = "[CURRENT]h[SELECTION]_nom"
mach3_pattern     = "[CURRENT][SELECTION]_[PARAM]_sig_[VAR]"
valor_pattern     = "Ratio_RHC1Rmu_[PARAM]_[VAR]"

# paired things have valor, then mach3
variations = [ ("p1","p1"),
               ("p3","p3"),
               ("m1","n1"),
               ("m3","n3") ]
               
mach3_current = "anu_"
mach3_selection = "sknumubar"

parameters = [( 'f_banffmaqe_rhc'          ,  'xsec_par_0'  ) ,
            #(                             ,  'xsec_par_1'  ) ,       # pF carbon
            #(                             ,  'xsec_par_2'  ) ,       # MEC carbon
            #(                             ,  'xsec_par_3'  ) ,       # Eb carbon
              ( 'f_banffpfermi_rhc'        ,  'xsec_par_4'  ) ,
              ( 'f_banffmec_rhc'           ,  'xsec_par_5'  ) ,
              ( 'f_banffbindE_rhc'         ,  'xsec_par_6'  ) ,
              ( 'f_banffca5_rhc'           ,  'xsec_par_7'  ) ,
              ( 'f_banffmares_rhc'         ,  'xsec_par_8'  ) ,
              ( 'f_banffbgres_rhc'         ,  'xsec_par_9'  ) ,
              ( 'f_banffnorm_nuetonumu_rhc',  'xsec_par_10' ) ,
              ( 'f_banffshapeCCoth_rhc'    ,  'xsec_par_11' ) ,
             #(                            ,  'xsec_par_12' ) ,       # CCcoh carbon
              ( 'f_banffnorm_cccoh_rhc'    ,  'xsec_par_13' ) ,
              ( 'f_banffnorm_nccoh_rhc'    ,  'xsec_par_14' ) ,
              ( 'f_banffnorm_ncoth_rhc'    ,  'xsec_par_15' ) ,
             #(                            ,  'xsec_par_16' ) ,
             #(                            ,  'xsec_par_17' ) ,
             #(                            ,  'xsec_par_18' ) ,
             #(                            ,  'xsec_par_19' ) , 
             #(                            ,  'xsec_par_20' ) , 
             #(                            ,  'xsec_par_21' ) , 
             #(                            ,  'xsec_par_22' ) , 
             #(                            ,  'xsec_par_23' ) , 
             # ( ,  'flux_par_50' ) , 
             # ( ,  'flux_par_51' ) , 
             # ( ,  'flux_par_52' ) , 
             # ( ,  'flux_par_53' ) , 
             # ( ,  'flux_par_54' ) , 
             # ( ,  'flux_par_55' ) , 
             # ( ,  'flux_par_56' ) , 
             # ( ,  'flux_par_57' ) , 
             # ( ,  'flux_par_58' ) , 
             # ( ,  'flux_par_59' ) , 
             # ( ,  'flux_par_60' ) , 
             # ( ,  'flux_par_61' ) , 
             # ( ,  'flux_par_62' ) , 
             # ( ,  'flux_par_63' ) , 
             # ( ,  'flux_par_64' ) , 
             # ( ,  'flux_par_65' ) , 
             # ( ,  'flux_par_66' ) , 
             # ( ,  'flux_par_67' ) , 
             # ( ,  'flux_par_68' ) , 
             # ( ,  'flux_par_69' ) , 
             # ( ,  'flux_par_70' ) , 
             # ( ,  'flux_par_71' ) , 
             # ( ,  'flux_par_72' ) , 
             # ( ,  'flux_par_73' ) , 
             # ( ,  'flux_par_74' ) , 
              ( 'f_banff00_rhc',  'flux_par_75' ) , 
              ( 'f_banff01_rhc',  'flux_par_76' ) , 
              ( 'f_banff02_rhc',  'flux_par_77' ) , 
              ( 'f_banff03_rhc',  'flux_par_78' ) , 
              ( 'f_banff04_rhc',  'flux_par_79' ) , 
              ( 'f_banff05_rhc',  'flux_par_80' ) , 
              ( 'f_banff06_rhc',  'flux_par_81' ) , 
              ( 'f_banff07_rhc',  'flux_par_82' ) , 
              ( 'f_banff08_rhc',  'flux_par_83' ) , 
              ( 'f_banff09_rhc',  'flux_par_84' ) , 
              ( 'f_banff10_rhc',  'flux_par_85' ) , 
              ( 'f_banff11_rhc',  'flux_par_86' ) , 
              ( 'f_banff12_rhc',  'flux_par_87' ) , 
              ( 'f_banff13_rhc',  'flux_par_88' ) , 
              ( 'f_banff14_rhc',  'flux_par_89' ) , 
              ( 'f_banff15_rhc',  'flux_par_90' ) , 
              ( 'f_banff16_rhc',  'flux_par_91' ) , 
              ( 'f_banff17_rhc',  'flux_par_92' ) , 
              ( 'f_banff18_rhc',  'flux_par_93' ) , 
              ( 'f_banff19_rhc',  'flux_par_94' ) , 
              ( 'f_banff20_rhc',  'flux_par_95' ) , 
              ( 'f_banff21_rhc',  'flux_par_96' ) , 
              ( 'f_banff22_rhc',  'flux_par_97' ) , 
              ( 'f_banff23_rhc',  'flux_par_98' ) , 
              ( 'f_banff24_rhc',  'flux_par_99' ) , 
             #( ,  'skdet_par_0' ) , 
             #( ,  'skdet_par_1' ) , 
             #( ,  'skdet_par_2' ) , 
             #( ,  'skdet_par_3' ) , 
             #( ,  'skdet_par_4' ) , 
             #( ,  'skdet_par_5' ) , 
             #( ,  'skdet_par_6' ) , 
             #( ,  'skdet_par_7' ) , 
             #( ,  'skdet_par_8' ) , 
             #( ,  'skdet_par_9' ) , 
             #( ,  'skdet_par_10') , 
             #( ,  'skdet_par_11') , 
             
             # ( ,  'skdet_par_12') , 
             # ( ,  'skdet_par_13') , 
             # ( ,  'skdet_par_14') , 
             # ( ,  'skdet_par_15') , 
             # ( ,  'skdet_par_16') , 
             # ( ,  'skdet_par_17') , 
             
             # ( ,  'skdet_par_18') , 
             # ( ,  'skdet_par_19') , 
             # ( ,  'skdet_par_20') , 
             # ( ,  'skdet_par_21') , 
             # ( ,  'skdet_par_22') , 
             # ( ,  'skdet_par_23') , 
             # ( ,  'skdet_par_24') , 
             # ( ,  'skdet_par_25') , 
             # ( ,  'skdet_par_26') , 
             # ( ,  'skdet_par_27') , 
             # ( ,  'skdet_par_28') , 
             # ( ,  'skdet_par_29') , 
              
              ( 'f_skdetfsi000_rhc',  'skdet_par_30') , 
              ( 'f_skdetfsi001_rhc',  'skdet_par_31') , 
              ( 'f_skdetfsi002_rhc',  'skdet_par_32') , 
              ( 'f_skdetfsi003_rhc',  'skdet_par_33') , 
              ( 'f_skdetfsi004_rhc',  'skdet_par_34') , 
              ( 'f_skdetfsi005_rhc',  'skdet_par_35') ,
              
              ( 'f_sk_e_scale_rhc' ,  'skdet_par_36') ]
              
      
line_colors = dict()
line_colors["m1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["m3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def main():
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)
  
  pthet_file = ROOT.TFile(pthet_filename)
  mach3_file = ROOT.TFile(mach3_filename)
  valor_file = ROOT.TFile(valor_filename)
  
  assert(mach3_file)
  assert(valor_file)
  
  mach3_nom_histname = mach3_nom_pattern
  for old,new in [("[CURRENT]"  , mach3_current  ),
                  ("[SELECTION]", mach3_selection) ]:
    mach3_nom_histname = mach3_nom_histname.replace(old,new)
                                                    
  mach3_nom_hist = mach3_file.Get(mach3_nom_histname)
  if not mach3_nom_hist:
    print "Failed to get MaCh3 nominal histogram ",mach3_nom_histname
    
  pthet_nom_hist = pthet_file.Get(mach3_nom_histname)
  if not pthet_nom_hist:
    print "Failed to get P-theta nominal histogram ",mach3_nom_histname
  
  mach3_nom_hist_binno = hist_to_bin_number(mach3_nom_hist)
  pthet_nom_hist_binno = hist_to_bin_number(pthet_nom_hist)
    
  valor_hists             = dict() # only valor names will be used as keys
  valor_hists_binno       = dict() # only valor names will be used as keys
  mach3_hists             = dict() # only valor names will be used as keys
  mach3_ratio_hists       = dict() # only valor names will be used as keys
  mach3_hists_binno       = dict() # only valor names will be used as keys
  mach3_ratio_hists_binno = dict() # only valor names will be used as keys
  pthet_hists             = dict() # only valor names will be used as keys
  pthet_ratio_hists       = dict() # only valor names will be used as keys
  pthet_hists_binno       = dict() # only valor names will be used as keys
  pthet_ratio_hists_binno = dict() # only valor names will be used as keys
  canvases                = dict()
  
  for iparam,(valor_param,mach3_param) in enumerate(parameters):
    
    valor_hists             [valor_param] = dict()
    valor_hists_binno       [valor_param] = dict()
    mach3_hists             [valor_param] = dict()
    mach3_ratio_hists       [valor_param] = dict()
    mach3_hists_binno       [valor_param] = dict()
    mach3_ratio_hists_binno [valor_param] = dict()
    pthet_hists             [valor_param] = dict()
    pthet_ratio_hists       [valor_param] = dict()
    pthet_hists_binno       [valor_param] = dict()
    pthet_ratio_hists_binno [valor_param] = dict()
    
    for valor_var,mach3_var in variations:
      
      valor_histname = valor_pattern
      for old,new in [("[PARAM]", valor_param),
                      ("[VAR]"  , valor_var  )] :
        valor_histname = valor_histname.replace(old,new)
                                              
      mach3_histname = mach3_pattern
      for old,new in [("[PARAM]"    , mach3_param    ),
                      ("[CURRENT]"  , mach3_current  ),
                      ("[SELECTION]", mach3_selection),
                      ("[VAR]"      , mach3_var      ) ]:
        mach3_histname = mach3_histname.replace(old,new)
                                             
      valor_hists[valor_param][valor_var] = valor_file.Get(valor_histname)
      
      if not valor_hists[valor_param][valor_var]:
        print "Failed to get VALOR histogram ",valor_histname
        
      mach3_hists[valor_param][valor_var] = mach3_file.Get(mach3_histname)
      
      if not mach3_hists[valor_param][valor_var]:
        print "Failed to get MaCh3 histogram ",mach3_histname
        
      pthet_hists[valor_param][valor_var] = pthet_file.Get(mach3_histname)
      
      if not pthet_hists[valor_param][valor_var]:
        print "Failed to get P-theta histogram ",mach3_histname
        
      valor_hists[valor_param][valor_var].SetLineWidth(1)
      valor_hists[valor_param][valor_var].SetLineColor(line_colors[valor_var])
        
      mach3_hists[valor_param][valor_var].SetName(valor_hists[valor_param][valor_var].GetName()+"_mach3")
      mach3_hists[valor_param][valor_var].SetLineStyle(2)
      mach3_hists[valor_param][valor_var].SetLineWidth(2)
      mach3_hists[valor_param][valor_var].SetLineColor(valor_hists[valor_param][valor_var].GetLineColor()+2)
      
      pthet_hists[valor_param][valor_var].SetName(valor_hists[valor_param][valor_var].GetName()+"_pthet")
      pthet_hists[valor_param][valor_var].SetLineStyle(3)
      pthet_hists[valor_param][valor_var].SetLineWidth(2)
      pthet_hists[valor_param][valor_var].SetLineColor(valor_hists[valor_param][valor_var].GetLineColor()+4)
      
      valor_hists_binno[valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_param][valor_var])
      mach3_hists_binno[valor_param][valor_var] = hist_to_bin_number(mach3_hists[valor_param][valor_var])
      pthet_hists_binno[valor_param][valor_var] = hist_to_bin_number(pthet_hists[valor_param][valor_var])
      
      mach3_ratio_hists      [valor_param][valor_var] = ratio_hist(mach3_hists      [valor_param][valor_var],mach3_nom_hist)
      mach3_ratio_hists_binno[valor_param][valor_var] = ratio_hist(mach3_hists_binno[valor_param][valor_var],mach3_nom_hist_binno)
      
      pthet_ratio_hists      [valor_param][valor_var] = ratio_hist(pthet_hists      [valor_param][valor_var],pthet_nom_hist)
      pthet_ratio_hists_binno[valor_param][valor_var] = ratio_hist(pthet_hists_binno[valor_param][valor_var],pthet_nom_hist_binno)
      
      # Make sure that if there's zero nominal events, we get a ratio of 1.
      for ibin in xrange(mach3_nom_hist.GetXaxis().GetNbins()):
        if abs(mach3_nom_hist.GetBinContent(ibin+1)) < 1E-5:
          mach3_ratio_hists      [valor_param][valor_var].SetBinContent(ibin+1,1)
          mach3_ratio_hists_binno[valor_param][valor_var].SetBinContent(ibin+1,1)
          
      for ibin in xrange(pthet_nom_hist.GetXaxis().GetNbins()):
        if abs(pthet_nom_hist.GetBinContent(ibin+1)) < 1E-5:
          pthet_ratio_hists      [valor_param][valor_var].SetBinContent(ibin+1,1)
          pthet_ratio_hists_binno[valor_param][valor_var].SetBinContent(ibin+1,1)
          
    # Start making those plots.
    pdf_printing = [""]*len(parameters)
    pdf_printing[ 0] = "("
    pdf_printing[-1] = ")"
    
    # Plot mach3 vs valor
    hist_limits = []
    for valor_var,mach3_var in variations:
      for bin_content in (hist_to_list(valor_hists      [valor_param][valor_var])+
                          hist_to_list(mach3_ratio_hists[valor_param][valor_var])):
        hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    canv_key = valor_param+"_valor_vs_mach3"
    canv_key_binno = canv_key+"_binno"
    
    canvases[canv_key] = ROOT.TCanvas(canv_key,canv_key,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists[valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
      valor_hists[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists[valor_param][valor_var].Draw(option)
      option = "same"
      mach3_ratio_hists[valor_param][valor_var].Draw(option)
    
    canvases[canv_key].Update()
    canvases[canv_key].Print("syst_var_valor_vs_mach3_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[canv_key_binno] = ROOT.TCanvas(canv_key_binno,canv_key_binno,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists_binno[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists_binno[valor_param][valor_var].GetXaxis().SetTitle("Bin number")
      valor_hists_binno[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists_binno[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists_binno[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists_binno[valor_param][valor_var].Draw(option)
      option = "same"
      mach3_ratio_hists_binno[valor_param][valor_var].Draw(option)
    
    canvases[canv_key_binno].Update()
    canvases[canv_key_binno].Print("syst_var_valor_vs_mach3_binno"+".pdf"+pdf_printing[iparam])
    
    # Plot ptheta vs valor
    hist_limits = []
    for valor_var,mach3_var in variations:
      for bin_content in (hist_to_list(valor_hists[valor_param][valor_var])+hist_to_list(pthet_ratio_hists[valor_param][valor_var])):
        hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    canv_key = valor_param+"_valor_vs_pthet"
    canv_key_binno = canv_key+"_binno"
    
    canvases[canv_key] = ROOT.TCanvas(canv_key,canv_key,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists[valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
      valor_hists[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_ratio_hists[valor_param][valor_var].Draw(option)
    
    canvases[canv_key].Update()
    canvases[canv_key].Print("syst_var_valor_vs_ptheta_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[canv_key_binno] = ROOT.TCanvas(canv_key_binno,canv_key_binno,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists_binno[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists_binno[valor_param][valor_var].GetXaxis().SetTitle("Bin number")
      valor_hists_binno[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists_binno[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists_binno[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists_binno[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_ratio_hists_binno[valor_param][valor_var].Draw(option)
    
    canvases[canv_key_binno].Update()
    canvases[canv_key_binno].Print("syst_var_valor_vs_ptheta_binno"+".pdf"+pdf_printing[iparam])
    
    # Plot pthet vs mach3
    hist_limits = []
    for valor_var,mach3_var in variations:
      for bin_content in (hist_to_list(mach3_ratio_hists[valor_param][valor_var])+
                          hist_to_list(pthet_ratio_hists[valor_param][valor_var])):
        hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    canv_key = valor_param+"_mach3_vs_pthet"
    canv_key_binno = canv_key+"_binno"
    
    canvases[canv_key] = ROOT.TCanvas(canv_key,canv_key,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      mach3_ratio_hists[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+mach3_param+" vs "+mach3_param)
      mach3_ratio_hists[valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
      mach3_ratio_hists[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      mach3_ratio_hists[valor_param][valor_var].SetMaximum(hist_max)
      mach3_ratio_hists[valor_param][valor_var].SetMinimum(hist_min)
      mach3_ratio_hists[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_ratio_hists[valor_param][valor_var].Draw(option)
    
    canvases[canv_key].Update()
    canvases[canv_key].Print("syst_var_mach3_vs_ptheta_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[canv_key_binno] = ROOT.TCanvas(canv_key_binno,canv_key_binno,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      mach3_ratio_hists_binno[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+mach3_param+" vs "+mach3_param)
      mach3_ratio_hists_binno[valor_param][valor_var].GetXaxis().SetTitle("Bin number")
      mach3_ratio_hists_binno[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      mach3_ratio_hists_binno[valor_param][valor_var].SetMaximum(hist_max)
      mach3_ratio_hists_binno[valor_param][valor_var].SetMinimum(hist_min)
      mach3_ratio_hists_binno[valor_param][valor_var].Draw(option)
      option = "same"
      pthet_ratio_hists_binno[valor_param][valor_var].Draw(option)
    
    canvases[canv_key_binno].Update()
    canvases[canv_key_binno].Print("syst_var_mach3_vs_ptheta_binno"+".pdf"+pdf_printing[iparam])
    
    # Plot all
    hist_limits = []
    for valor_var,mach3_var in variations:
      for bin_content in (hist_to_list(valor_hists      [valor_param][valor_var]) +
                          hist_to_list(mach3_ratio_hists[valor_param][valor_var]) +
                          hist_to_list(pthet_ratio_hists[valor_param][valor_var])):
          hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    canv_key = valor_param+"_valor_vs_mach3_vs_pthet"
    canv_key_binno = canv_key+"_binno"
    
    canvases[canv_key] = ROOT.TCanvas(canv_key,canv_key,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists[valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
      valor_hists[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists[valor_param][valor_var].Draw(option)
      option = "same"
      mach3_ratio_hists[valor_param][valor_var].Draw(option)
      pthet_ratio_hists[valor_param][valor_var].Draw(option)
    
    canvases[canv_key].Update()
    canvases[canv_key].Print("syst_var_valor_vs_mach3_vs_ptheta_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[canv_key_binno] = ROOT.TCanvas(canv_key_binno,canv_key_binno,1280,720)
    option=""
    for valor_var,mach3_var in variations:
      valor_hists_binno[valor_param][valor_var].SetTitle("Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
      valor_hists_binno[valor_param][valor_var].GetXaxis().SetTitle("Bin number")
      valor_hists_binno[valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
      valor_hists_binno[valor_param][valor_var].SetMaximum(hist_max)
      valor_hists_binno[valor_param][valor_var].SetMinimum(hist_min)
      valor_hists_binno[valor_param][valor_var].Draw(option)
      option = "same"
      mach3_ratio_hists_binno[valor_param][valor_var].Draw(option)
      pthet_ratio_hists_binno[valor_param][valor_var].Draw(option)
    
    canvases[canv_key_binno].Update()
    canvases[canv_key_binno].Print("syst_var_valor_vs_mach3_vs_ptheta_binno"+".pdf"+pdf_printing[iparam])
    
    
  return

#_______________________________________________________________________

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

