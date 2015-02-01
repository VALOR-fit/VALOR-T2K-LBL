#!/usr/bin/python
import sys
import ROOT 
import os

env_var = "T2KOSCANA"

pthet_filename = '/home/steve/analysis/systematic_variations/systvariation.root'

mach3_filename = '/home/steve/analysis/systematic_variations/sigmavar2014SK.root'

pthet_nom_pattern = "[CURRENT]h[SELECTION]_nom"
pthet_pattern     = "[CURRENT][SELECTION]_[PARAM]_sig_[VAR]"
mach3_nom_pattern = "[CURRENT]h[SELECTION]_nom"
mach3_pattern     = "[CURRENT][SELECTION]_[PARAM]_sig_[VAR]"

# paired things have valor, then pthet
variations = [ ("p1","p1"),
               ("p3","p3"),
               ("n1","n1"),
               ("n3","n3") ]
               
pthet_current = "anu_"
pthet_selection = "sknumubar"
mach3_current = "anu_"
mach3_selection = "sknumubar"

parameters = [( 'xsec_par_0'   ,  'xsec_par_0'  ) ,
              ( 'xsec_par_4'   ,  'xsec_par_4'  ) ,
              ( 'xsec_par_5'   ,  'xsec_par_5'  ) ,
              ( 'xsec_par_6'   ,  'xsec_par_6'  ) ,
              ( 'xsec_par_7'   ,  'xsec_par_7'  ) ,
              ( 'xsec_par_8'   ,  'xsec_par_8'  ) ,
              ( 'xsec_par_9'   ,  'xsec_par_9'  ) ,
              ( 'xsec_par_10'  ,  'xsec_par_10' ) ,
              ( 'xsec_par_11'  ,  'xsec_par_11' ) ,
              ( 'xsec_par_13'  ,  'xsec_par_13' ) ,
              ( 'xsec_par_14'  ,  'xsec_par_14' ) ,
              ( 'xsec_par_15'  ,  'xsec_par_15' ) ,
              ( 'flux_par_75'  ,  'flux_par_75' ) , 
              ( 'flux_par_76'  ,  'flux_par_76' ) , 
              ( 'flux_par_77'  ,  'flux_par_77' ) , 
              ( 'flux_par_78'  ,  'flux_par_78' ) , 
              ( 'flux_par_79'  ,  'flux_par_79' ) , 
              ( 'flux_par_80'  ,  'flux_par_80' ) , 
              ( 'flux_par_81'  ,  'flux_par_81' ) , 
              ( 'flux_par_82'  ,  'flux_par_82' ) , 
              ( 'flux_par_83'  ,  'flux_par_83' ) , 
              ( 'flux_par_84'  ,  'flux_par_84' ) , 
              ( 'flux_par_85'  ,  'flux_par_85' ) , 
              ( 'flux_par_86'  ,  'flux_par_86' ) , 
              ( 'flux_par_87'  ,  'flux_par_87' ) , 
              ( 'flux_par_88'  ,  'flux_par_88' ) , 
              ( 'flux_par_89'  ,  'flux_par_89' ) , 
              ( 'flux_par_90'  ,  'flux_par_90' ) , 
              ( 'flux_par_91'  ,  'flux_par_91' ) , 
              ( 'flux_par_92'  ,  'flux_par_92' ) , 
              ( 'flux_par_93'  ,  'flux_par_93' ) , 
              ( 'flux_par_94'  ,  'flux_par_94' ) , 
              ( 'flux_par_95'  ,  'flux_par_95' ) , 
              ( 'flux_par_96'  ,  'flux_par_96' ) , 
              ( 'flux_par_97'  ,  'flux_par_97' ) , 
              ( 'flux_par_98'  ,  'flux_par_98' ) , 
              ( 'flux_par_99'  ,  'flux_par_99' ) , 
              ( 'skdet_par_30' ,  'skdet_par_30') , 
              ( 'skdet_par_31' ,  'skdet_par_31') , 
              ( 'skdet_par_32' ,  'skdet_par_32') , 
              ( 'skdet_par_33' ,  'skdet_par_33') , 
              ( 'skdet_par_34' ,  'skdet_par_34') , 
              ( 'skdet_par_35' ,  'skdet_par_35') ,
              ( 'skdet_par_36' ,  'skdet_par_36') ]
              
      
line_colors = dict()
line_colors["n1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["n3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def main():
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)
  
  pthet_file = ROOT.TFile(pthet_filename)
  mach3_file = ROOT.TFile(mach3_filename)
  
  assert(pthet_file)
  assert(mach3_file)
  
  pthet_nom_histname = pthet_nom_pattern
  for old,new in [("[CURRENT]"  , pthet_current  ),
                  ("[SELECTION]", pthet_selection) ]:
    pthet_nom_histname = pthet_nom_histname.replace(old,new)
                                                    
  pthet_nom_hist = pthet_file.Get(pthet_nom_histname)
  if not pthet_nom_hist:
    print "Failed to get P-theta nominal histogram ",pthet_nom_histname
    
  mach3_nom_histname = mach3_nom_pattern
  for old,new in [("[CURRENT]"  , mach3_current  ),
                  ("[SELECTION]", mach3_selection) ]:
    mach3_nom_histname = mach3_nom_histname.replace(old,new)
                                                    
  mach3_nom_hist = mach3_file.Get(mach3_nom_histname)
  if not mach3_nom_hist:
    print "Failed to get mach3 nominal histogram ",mach3_nom_histname
  
  pthet_nom_hist_binno = hist_to_bin_number(pthet_nom_hist)
  mach3_nom_hist_binno = hist_to_bin_number(mach3_nom_hist)
    
  mach3_hists             = dict() # only mach3 names will be used as keys
  pthet_hists             = dict() # only mach3 names will be used as keys
  pthet_ratio_hists       = dict() # only mach3 names will be used as keys
  mach3_ratio_hists       = dict() # only mach3 names will be used as keys
  mach3_hists_binno       = dict() # only mach3 names will be used as keys
  pthet_hists_binno       = dict() # only mach3 names will be used as keys
  pthet_ratio_hists_binno = dict() # only mach3 names will be used as keys
  mach3_ratio_hists_binno = dict() # only mach3 names will be used as keys
  canvases                = dict()
  
  for iparam,(mach3_param,pthet_param) in enumerate(parameters):
    
    mach3_hists             [mach3_param] = dict()
    pthet_hists             [mach3_param] = dict()
    mach3_ratio_hists       [mach3_param] = dict()
    pthet_ratio_hists       [mach3_param] = dict()
    mach3_hists_binno       [mach3_param] = dict()
    pthet_hists_binno       [mach3_param] = dict()
    mach3_ratio_hists_binno [mach3_param] = dict()
    pthet_ratio_hists_binno [mach3_param] = dict()
    
    for mach3_var,pthet_var in variations:
      
      print mach3_var
                                              
      mach3_histname = mach3_pattern
      for old,new in [("[PARAM]"    , mach3_param    ),
                      ("[CURRENT]"  , mach3_current  ),
                      ("[SELECTION]", mach3_selection),
                      ("[VAR]"      , mach3_var      ) ]:
        mach3_histname = mach3_histname.replace(old,new)
        
      pthet_histname = pthet_pattern
      for old,new in [("[PARAM]"    , pthet_param    ),
                      ("[CURRENT]"  , pthet_current  ),
                      ("[SELECTION]", pthet_selection),
                      ("[VAR]"      , pthet_var      ) ]:
        pthet_histname = pthet_histname.replace(old,new)
                                             
      mach3_hists[mach3_param][mach3_var] = mach3_file.Get(mach3_histname)
      
      if not mach3_hists[mach3_param][mach3_var]:
        print "Failed to get VALOR histogram ",mach3_histname
        
      pthet_hists[mach3_param][mach3_var] = pthet_file.Get(pthet_histname)
      if not pthet_hists[mach3_param][mach3_var]:
        print "Failed to get MaCh3 histogram ",pthet_histname
        
      mach3_hists[mach3_param][mach3_var].SetLineWidth(1)
      mach3_hists[mach3_param][mach3_var].SetLineStyle(2)
      mach3_hists[mach3_param][mach3_var].SetLineColor(line_colors[mach3_var])
        
      pthet_hists[mach3_param][mach3_var].SetName(mach3_hists[mach3_param][mach3_var].GetName()+"_pthet")
      pthet_hists[mach3_param][mach3_var].SetLineStyle(3)
      pthet_hists[mach3_param][mach3_var].SetLineWidth(2)
      pthet_hists[mach3_param][mach3_var].SetLineColor(mach3_hists[mach3_param][mach3_var].GetLineColor()+2)
      
      pthet_hists_binno[mach3_param][mach3_var] = hist_to_bin_number(pthet_hists[mach3_param][mach3_var])
      pthet_ratio_hists      [mach3_param][mach3_var] = ratio_hist(pthet_hists      [mach3_param][mach3_var],pthet_nom_hist)
      pthet_ratio_hists_binno[mach3_param][mach3_var] = ratio_hist(pthet_hists_binno[mach3_param][mach3_var],pthet_nom_hist_binno)
      
      mach3_hists_binno[mach3_param][mach3_var] = hist_to_bin_number(mach3_hists[mach3_param][mach3_var])
      mach3_ratio_hists      [mach3_param][mach3_var] = ratio_hist(mach3_hists      [mach3_param][mach3_var],mach3_nom_hist)
      mach3_ratio_hists_binno[mach3_param][mach3_var] = ratio_hist(mach3_hists_binno[mach3_param][mach3_var],mach3_nom_hist_binno)
    
    hist_limits = []
    for mach3_var,pthet_var in variations:
      for bin_content in (hist_to_list(mach3_ratio_hists[mach3_param][mach3_var])+hist_to_list(pthet_ratio_hists[mach3_param][mach3_var])):
        if abs(bin_content) > 1E-5:
          hist_limits.append(abs(bin_content-1))
    
    hist_var = 1.1*max(hist_limits)
    hist_max = hist_var + 1
    hist_min = 1-hist_var
    
    pdf_printing = [""]*len(parameters)
    pdf_printing[ 0] = "("
    pdf_printing[-1] = ")"
    
    canvases[mach3_param] = ROOT.TCanvas(mach3_param,mach3_param,1280,720)
    option=""
    for mach3_var,pthet_var in variations:
      mach3_ratio_hists[mach3_param][mach3_var].SetTitle("Systematic tweak ratio: "+mach3_param+" vs "+pthet_param)
      mach3_ratio_hists[mach3_param][mach3_var].GetXaxis().SetTitle("E_{reco}")
      mach3_ratio_hists[mach3_param][mach3_var].GetYaxis().SetTitle("Ratio to nominal")
      mach3_ratio_hists[mach3_param][mach3_var].SetMaximum(hist_max)
      mach3_ratio_hists[mach3_param][mach3_var].SetMinimum(hist_min)
      mach3_ratio_hists[mach3_param][mach3_var].Draw(option)
      option = "same"
      pthet_ratio_hists[mach3_param][mach3_var].Draw(option)
    
    canvases[mach3_param].Update()
    canvases[mach3_param].Print("syst_var_mach3_vs_ptheta_ereco"+".pdf"+pdf_printing[iparam])
    
    canvases[mach3_param+"_binno"] = ROOT.TCanvas(mach3_param+"_binno",mach3_param+"_binno",1280,720)
    option=""
    for mach3_var,pthet_var in variations:
      mach3_ratio_hists_binno[mach3_param][mach3_var].SetTitle("Systematic tweak ratio: "+mach3_param+" vs "+pthet_param)
      mach3_ratio_hists_binno[mach3_param][mach3_var].GetXaxis().SetTitle("Bin number")
      mach3_ratio_hists_binno[mach3_param][mach3_var].GetYaxis().SetTitle("Ratio to nominal")
      mach3_ratio_hists_binno[mach3_param][mach3_var].SetMaximum(hist_max)
      mach3_ratio_hists_binno[mach3_param][mach3_var].SetMinimum(hist_min)
      mach3_ratio_hists_binno[mach3_param][mach3_var].Draw(option)
      option = "same"
      pthet_ratio_hists_binno[mach3_param][mach3_var].Draw(option)
    
    canvases[mach3_param+"_binno"].Update()
    canvases[mach3_param+"_binno"].Print("syst_var_mach3_vs_ptheta_binno"+".pdf"+pdf_printing[iparam])
    

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

