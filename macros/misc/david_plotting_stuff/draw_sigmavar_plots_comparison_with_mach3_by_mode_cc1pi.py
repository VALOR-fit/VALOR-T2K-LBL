#!/usr/bin/python
import sys
import ROOT
import os

env_var = "T2KOSCANA"

parameters = [

  # VALOR parameter name T2KPhysParamId , MaCh3 Parameter names NEUT/NIWG 2015c
  
  ( 'f_banffmaqe_fhc_rhc'          ,  'xsec_par_0'  ) ,    #   0: MaQE                               
  #(                                ,  'xsec_par_1'  ) ,    #   1: pF (C)                               
  #(                                ,  'xsec_par_2'  ) ,    #   2: MEC (C)                            
  #(                                ,  'xsec_par_3'  ) ,    #   3: Eb (C)                                                                                                         
  ( 'f_banffpfermi_fhc_rhc'        ,  'xsec_par_4'  ) ,    #   4: pF (O)                             
  ( 'f_banffmec_fhc_rhc'           ,  'xsec_par_5'  ) ,    #   5: MEC (O)                            
  ( 'f_banffbindE_fhc_rhc'         ,  'xsec_par_6'  ) ,    #   6: Eb (O)
  ( 'f_banffmecbar_fhc_rhc'        ,  'xsec_par_7'  ) ,    #   7: $\bar{\nu}$ MEC (O)                
  ( 'f_banffca5_fhc_rhc'           ,  'xsec_par_8'  ) ,    #   8: CA5                                   
  ( 'f_banffmares_fhc_rhc'         ,  'xsec_par_9'  ) ,    #   9: MaRES                              
  ( 'f_banffbgres_fhc_rhc'         ,  'xsec_par_10' ) ,    #  10: Isospin Bkg.                       
  ( 'f_banffnorm_nuetonumu_fhc_rhc',  'xsec_par_11' ) ,    #  11: $\nu_{e}/\nu_{\mu}$                
  ( 'f_banffnorm_nuebartonumubar_fhc_rhc'  ,  'xsec_par_12' ) ,    #  12: $\bar{\nu_{e}}/\bar{\nu_{\mu}}$                                                                              
  ( 'f_banffshapeCCoth_fhc_rhc'    ,  'xsec_par_13' ) ,    #  13: CC Oth Shape                                                                                                                                                                           
  #(     ,  'xsec_par_14' ) ,    #  14: CC Coh (C)
  ( 'f_banffnorm_cccoh_fhc_rhc'    ,  'xsec_par_15' ) ,    #  15: CC Coh (O)                                            
  ( 'f_banffnorm_nccoh_fhc_rhc'    ,  'xsec_par_16' ) ,    #  16: NC Coh                                                
  ( 'f_banffnorm_nc1gamma_fhc_rhc' ,  'xsec_par_17' ) ,    #  17: NC 1 Gamma
  #(                                ,  'xsec_par_18' ) ,    #  18: NC Oth Near                                           
  ( 'f_banffnorm_ncoth_fhc_rhc'     ,  'xsec_par_19' ) ,    #  19: NC Oth Far                                            
  #(                                ,  'xsec_par_20' ) ,    #  20: FSI Inel. Low E                                       
  #(                                ,  'xsec_par_21' ) ,    #  21: FSI Inel. High E                                      
  #(                                ,  'xsec_par_22' ) ,    #  22: FSI Pion Prod.                                        
  #(                                ,  'xsec_par_23' ) ,    #  23: FSI Pion Abs.                                         
  #(                                ,  'xsec_par_24' ) ,    #  24: FSI Ch. Exch. Low E                                   
  #(                                ,  'xsec_par_25' ) ,    #  25: FSI Ch. Exch. High E                                  
  
  ( 'f_banff00' ,  'flux_par_50' ) ,
  ( 'f_banff01' ,  'flux_par_51' ) ,
  ( 'f_banff02' ,  'flux_par_52' ) ,
  ( 'f_banff03' ,  'flux_par_53' ) ,
  ( 'f_banff04' ,  'flux_par_54' ) ,
  ( 'f_banff05' ,  'flux_par_55' ) ,
  ( 'f_banff06' ,  'flux_par_56' ) ,
  ( 'f_banff07' ,  'flux_par_57' ) ,
  ( 'f_banff08' ,  'flux_par_58' ) ,
  ( 'f_banff09' ,  'flux_par_59' ) ,
  ( 'f_banff10' ,  'flux_par_60' ) ,
  ( 'f_banff11' ,  'flux_par_61' ) ,
  ( 'f_banff12' ,  'flux_par_62' ) ,
  ( 'f_banff13' ,  'flux_par_63' ) ,
  ( 'f_banff14' ,  'flux_par_64' ) ,
  ( 'f_banff15' ,  'flux_par_65' ) ,
  ( 'f_banff16' ,  'flux_par_66' ) ,
  ( 'f_banff17' ,  'flux_par_67' ) ,
  ( 'f_banff18' ,  'flux_par_68' ) ,
  ( 'f_banff19' ,  'flux_par_69' ) ,
  ( 'f_banff20' ,  'flux_par_70' ) ,
  ( 'f_banff21' ,  'flux_par_71' ) ,
  ( 'f_banff22' ,  'flux_par_72' ) ,
  ( 'f_banff23' ,  'flux_par_73' ) ,
  ( 'f_banff24' ,  'flux_par_74' ) ,

  # Have no entries in lowest energy bins see SKFSISystConv.h
  # kListOfMultiRingSKDetFSIParams
  
  #( 'f_skdetfsi006_multiring',  '') , 
  ( 'f_skdetfsi007_multiring',  'skdet_par_36') ,
  ( 'f_skdetfsi008_multiring',  'skdet_par_37') ,
  #( 'f_skdetfsi009_multiring',  '') ,
  ( 'f_skdetfsi010_multiring',  'skdet_par_38') ,
  ( 'f_skdetfsi011_multiring',  'skdet_par_39') ,
  #( 'f_skdetfsi012_multiring',  '') ,
  ( 'f_skdetfsi013_multiring',  'skdet_par_40') ,
  ( 'f_skdetfsi014_multiring',  'skdet_par_41') ,
  #( 'f_skdetfsi015_multiring',  '') ,
  ( 'f_skdetfsi016_multiring',  'skdet_par_42') ,
  ( 'f_skdetfsi017_multiring',  'skdet_par_43') ,

  ( 'f_sk_e_scale' ,  'skdet_par_44'),

]

if single :
  modes = [
           ( ''                 , ''                  )
           ]
else :
  modes = [
           ( ''                , ''                  ) ,
		  ( '_numuccqe'        , "_sm_numu_ccqe"      ) ,  
		  ( '_numucc1pi'       , "_sm_numu_cc1pi"     ) ,
		  ( '_numucccoh'       , "_sm_numu_cccoh"     ) ,
		  ( '_numuccmec'       , "_sm_numu_mec"       ) ,
		  ( '_numuccoth'       , "_sm_numu_ccoth"     ) ,
		  ( '_numunc1piC'      , "_sm_numu_ncpipm"    ) ,
		  ( '_numunc1pi0'      , "_sm_numu_ncpi0"     ) ,
		  ( '_numunccoh'       , "_sm_numu_nccoh"     ) ,
		  ( '_numuncoth'       , "_sm_numu_ncoth"     ) ,
		  ( '_numunc1gamma'    , "_sm_numu_nc1gamma"  ) , # Adding missing NC 1gamma

		  ( '_numubarccqe'     , "_sm_numub_ccqe"     ) ,
		  ( '_numubarcc1pi'    , "_sm_numub_cc1pi"    ) ,
		  ( '_numubarcccoh'    , "_sm_numub_cccoh"    ) ,
		  ( '_numubarccmec'    , "_sm_numub_mec"      ) ,
		  ( '_numubarccoth'    , "_sm_numub_ccoth"    ) ,
		  ( '_numubarnc1piC'   , "_sm_numub_ncpipm"   ) ,
		  ( '_numubarnc1pi0'   , "_sm_numub_ncpi0"    ) ,
		  ( '_numubarnccoh'    , "_sm_numub_nccoh"    ) ,
		  ( '_numubarncoth'    , "_sm_numub_ncoth"    ) ,
		  ( '_numubarnc1gamma' , "_sm_numub_nc1gamma" ) , # Adding missing NC 1gamma  

		  ( '_nueccqe'         , "_sm_nue_ccqe"       ) ,
		  ( '_nuecc1pi'        , "_sm_nue_cc1pi"      ) ,
		  ( '_nuecccoh'        , "_sm_nue_cccoh"      ) ,
		  ( '_nueccmec'        , "_sm_nue_mec"        ) ,
		  ( '_nueccoth'        , "_sm_nue_ccoth"      ) ,
		  ( '_nuenc1piC'       , "_sm_nue_ncpipm"     ) ,
		  ( '_nuenc1pi0'       , "_sm_nue_ncpi0"      ) ,
		  ( '_nuenccoh'        , "_sm_nue_nccoh"      ) ,
		  ( '_nuencoth'        , "_sm_nue_ncoth"      ) ,
		  ( '_nuenc1gamma'     , "_sm_nue_nc1gamma"   ) , # Adding missing NC 1gamma  

		  ( '_oscnueccqe'      , "_sm_signue_ccqe"    ) ,
		  ( '_oscnuecc1pi'     , "_sm_signue_cc1pi"   ) ,
		  ( '_oscnuecccoh'     , "_sm_signue_cccoh"   ) ,
		  ( '_oscnueccmec'     , "_sm_signue_mec"     ) ,
		  ( '_oscnueccoth'     , "_sm_signue_ccoth"   ) ,

		  ( '_nuebarccqe'      , "_sm_nueb_ccqe"      ) ,
		  ( '_nuebarcc1pi'     , "_sm_nueb_cc1pi"     ) ,
		  ( '_nuebarcccoh'     , "_sm_nueb_cccoh"     ) ,
		  ( '_nuebarccmec'     , "_sm_nueb_mec"       ) ,
		  ( '_nuebarccoth'     , "_sm_nueb_ccoth"     ) ,
		  ( '_nuebarnc1piC'    , "_sm_nueb_ncpipm"    ) ,
		  ( '_nuebarnc1pi0'    , "_sm_nueb_ncpi0"     ) ,
		  ( '_nuebarnccoh'     , "_sm_nueb_nccoh"     ) ,
		  ( '_nuebarncoth'     , "_sm_nueb_ncoth"     ) ,
		  ( '_nuebarnc1gamma'  , "_sm_nueb_nc1gamma"  ) , # Adding missing NC 1gamma  

		  ( '_oscnuebarccqe'   , "_sm_signueb_ccqe"   ) ,
		  ( '_oscnuebarcc1pi'  , "_sm_signueb_cc1pi"  ) ,
		  ( '_oscnuebarcccoh'  , "_sm_signueb_cccoh"  ) ,
		  ( '_oscnuebarccmec'  , "_sm_signueb_mec"    ) ,
		  ( '_oscnuebarccoth'  , "_sm_signueb_ccoth"  ) ]

line_colors = dict()
line_colors["m1"] = ROOT.kViolet
line_colors["p1"] = ROOT.kCyan
line_colors["m3"] = ROOT.kRed
line_colors["p3"] = ROOT.kBlue

def main(mach3_filename,valor_filename,mach3_selection,tag,sample):

  mach3_nom_pattern = "[CURRENT][SELECTION][MODE]_nom"
  mach3_pattern     = "[CURRENT][SELECTION]_[PARAM][MODE]_sig_[VAR]"
#  valor_nom_pattern = "Spectrum_"+sample+"_nominal"
#  valor_nom_pattern     = "Spectrum_"+sample+"[MODE]_[PARAM]_[VAR]" 
  valor_pattern     = "Ratio_"+sample+"[MODE]_[PARAM]_[VAR]" 
#  valor_pattern     = valor_nom_pattern
  
  # paired objects have valor, then mach3
  variations = [ ("p1","p1"),
                 ("p3","p3"),
                 ("m1","n1"),
                 ("m3","n3") ]

  mach3_current = "" #"anu_" # only looking at the FHC current samples, and are not labelled in the current round
  mach3_selection_1 = mach3_selection # nue CC 1pi sample "sknue1pi" "sknumu"
  
  ROOT.gROOT.SetBatch(1)
  ROOT.gStyle.SetOptStat(0)

  mach3_file = ROOT.TFile(mach3_filename)
  valor_file = ROOT.TFile(valor_filename)

  assert(mach3_file)
  assert(valor_file)

  # valor_nom_histname = valor_nom_pattern
  
  # valor_nom_hist = valor_file.Get(valor_nom_histname)
  # if not valor_nom_hist:
  #   print "Failed to get VALOR nominal histogram ", valor_nom_histname

  # valor_nom_hist_binno = hist_to_bin_number(valor_nom_hist)
  
  mach3_nom_hists         = dict() # only valor names will be used as keys
  mach3_nom_hists_binno   = dict() # only valor names will be used as keys
  
  valor_hists             = dict() # only valor names will be used as keys

  mach3_hists             = dict() # only valor names will be used as keys
  mach3_ratio_hists       = dict() # only valor names will be used as keys

  valor_hists_binno       = dict() # only valor names will be used as keys
  mach3_hists_binno       = dict() # only valor names will be used as keys
  mach3_ratio_hists_binno = dict() # only valor names will be used as keys

  canvases                = dict()

  pdf_printing = [""]*(len(parameters)*len(modes))
  pdf_printing[ 0] = "("
  pdf_printing[-1] = ")"
  pdf_printing_iter_ereco = iter(pdf_printing)
  pdf_printing_iter_binno = iter(pdf_printing)

  for imode,(valor_mode,mach3_mode) in enumerate(modes):
    valor_hists             [valor_mode] = dict()
    mach3_hists             [valor_mode] = dict()
    mach3_ratio_hists       [valor_mode] = dict()
    valor_hists_binno       [valor_mode] = dict()
    mach3_hists_binno       [valor_mode] = dict()
    mach3_ratio_hists_binno [valor_mode] = dict()

    if mach3_mode == "" :
      mach3_selection = "h" + mach3_selection_1
    else :
      mach3_selection = mach3_selection_1

    mach3_nom_histname = mach3_nom_pattern

    for old,new in [("[CURRENT]"  , mach3_current  ),
                    ("[MODE]"     , mach3_mode     ),
                    ("[SELECTION]", mach3_selection) ]:
      mach3_nom_histname = mach3_nom_histname.replace(old,new)

    mach3_nom_hists[valor_mode] = mach3_file.Get(mach3_nom_histname)
    if not mach3_nom_hists[valor_mode]:
      print "Failed to get MaCh3 nominal histogram ",mach3_nom_histname

    print mach3_nom_hists[valor_mode]

    mach3_nom_hists_binno[valor_mode] = hist_to_bin_number(mach3_nom_hists[valor_mode])


    for iparam,(valor_param,mach3_param) in enumerate(parameters):

      valor_hists             [valor_mode][valor_param] = dict()
      mach3_hists             [valor_mode][valor_param] = dict()
      mach3_ratio_hists       [valor_mode][valor_param] = dict()
      valor_hists_binno       [valor_mode][valor_param] = dict()
      mach3_hists_binno       [valor_mode][valor_param] = dict()
      mach3_ratio_hists_binno [valor_mode][valor_param] = dict()

      for valor_var,mach3_var in variations:

        valor_histname = valor_pattern
        for old,new in [("[PARAM]", valor_param),
                        ("[MODE]" , valor_mode ),
                        ("[VAR]"  , valor_var  )] :
          valor_histname = valor_histname.replace(old,new)

        mach3_histname = mach3_pattern
        for old,new in [("[PARAM]"    , mach3_param    ),
                        ("[CURRENT]"  , mach3_current  ),
                        ("[MODE]"     , mach3_mode     ),
                        ("[SELECTION]", mach3_selection_1),
                        ("[VAR]"      , mach3_var      ) ]:
          mach3_histname = mach3_histname.replace(old,new)

        valor_hists[valor_mode][valor_param][valor_var] = valor_file.Get(valor_histname)

        if not valor_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get VALOR histogram ",valor_histname

        mach3_hists[valor_mode][valor_param][valor_var] = mach3_file.Get(mach3_histname)
        if not mach3_hists[valor_mode][valor_param][valor_var]:
          print "Failed to get MaCh3 histogram ",mach3_histname

        valor_hists[valor_mode][valor_param][valor_var].SetLineWidth(1)
        valor_hists[valor_mode][valor_param][valor_var].SetLineColor(line_colors[valor_var])

        mach3_hists[valor_mode][valor_param][valor_var].SetName(valor_hists[valor_mode][valor_param][valor_var].GetName()+"_mach3")
        mach3_hists[valor_mode][valor_param][valor_var].SetLineStyle(2)
        mach3_hists[valor_mode][valor_param][valor_var].SetLineWidth(2)
        mach3_hists[valor_mode][valor_param][valor_var].SetLineColor(valor_hists[valor_mode][valor_param][valor_var].GetLineColor()+2)

        valor_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(valor_hists[valor_mode][valor_param][valor_var])
        mach3_hists_binno[valor_mode][valor_param][valor_var] = hist_to_bin_number(mach3_hists[valor_mode][valor_param][valor_var])

        mach3_ratio_hists      [valor_mode][valor_param][valor_var] = ratio_hist(mach3_hists      [valor_mode][valor_param][valor_var],mach3_nom_hists[valor_mode])
        mach3_ratio_hists_binno[valor_mode][valor_param][valor_var] = ratio_hist(mach3_hists_binno[valor_mode][valor_param][valor_var],mach3_nom_hists_binno[valor_mode])

        set_zero_to_one(mach3_ratio_hists      [valor_mode][valor_param][valor_var])
        set_zero_to_one(mach3_ratio_hists_binno[valor_mode][valor_param][valor_var])

      hist_limits = []
      for valor_var,mach3_var in variations:
        for bin_content in (hist_to_list(valor_hists[valor_mode][valor_param][valor_var])+hist_to_list(mach3_ratio_hists[valor_mode][valor_param][valor_var])):
          if abs(bin_content) > 1E-5:
            hist_limits.append(abs(bin_content-1))

      hist_var = 1.1*max(hist_limits)
      hist_max = hist_var + 1
      hist_min = 1-hist_var

      canv_name = valor_mode+"_"+valor_param

      canvases[canv_name] = ROOT.TCanvas(canv_name,canv_name,1280,720)
      option=""
      for valor_var,mach3_var in variations:
        valor_hists[valor_mode][valor_param][valor_var].SetTitle(mach3_mode+" Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
        valor_hists[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("E_{reco}")
        valor_hists[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        mach3_ratio_hists[valor_mode][valor_param][valor_var].Draw(option)
#        mach3_hists[valor_mode][valor_param][valor_var].Draw(option)

      canvases[canv_name].SetLeftMargin(0.125)   
      canvases[canv_name].SetBottomMargin(0.125)          
      canvases[canv_name].Update()
      canvases[canv_name].Print("syst_var_mode_ereco"+tag+sample+".pdf"+next(pdf_printing_iter_ereco))

      canv_name_binno = canv_name+"_binno"

      canvases[canv_name_binno] = ROOT.TCanvas(canv_name_binno,canv_name_binno,1280,720)
      option=""
      for valor_var,mach3_var in variations:
        valor_hists_binno[valor_mode][valor_param][valor_var].SetTitle(mach3_mode+" Systematic tweak ratio: "+valor_param+" vs "+mach3_param)
        valor_hists_binno[valor_mode][valor_param][valor_var].GetXaxis().SetTitle("Bin number")
        valor_hists_binno[valor_mode][valor_param][valor_var].GetYaxis().SetTitle("Ratio to nominal")
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMaximum(hist_max)
        valor_hists_binno[valor_mode][valor_param][valor_var].SetMinimum(hist_min)
        valor_hists_binno[valor_mode][valor_param][valor_var].Draw(option)
        option = "same"
        mach3_ratio_hists_binno[valor_mode][valor_param][valor_var].Draw(option)

      canvases[canv_name_binno].Update()
      canvases[canv_name_binno].Print("syst_var_mode_binno"+tag+sample+".pdf"+next(pdf_printing_iter_binno))

  print "Done!"
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

def set_zero_to_one(h):
  nbins = h.GetXaxis().GetNbins()
  for ibin in xrange(1,nbins+1):
    if abs(h.GetBinContent(ibin) < 1E-10):
      h.SetBinContent(ibin,1)
  return

def hist_to_list(h):
  nbins = h.GetXaxis().GetNbins()
  return [h.GetBinContent(ibin) for ibin in xrange(1,nbins+1)]


if __name__=="__main__":

  import argparse
  parser = argparse.ArgumentParser(description="Get input files from the command line")
  parser.add_argument('-m', '--mach3-file', dest='mach3filename', help="MaCh3 file containing ratio of variations")
  parser.add_argument('-v', '--valor-file', dest='valorfilename', help="VALOR file containing ratio of variations")
  parser.add_argument('-s', '--mach3-selection', dest='mach3selection', help="MaCh3 selection")

  parser.add_argument('-t', '--tag', dest='tag', help='oscillated or unoscillated')
  parser.add_argument('-r', '--valor-sample', dest='sample', help='valor sample')
  parser.add_argument('-f', '--full-spectrum', dest='single', help='only draw full spectrum comparison')
  
  args = parser.parse_args()

  main(args.mach3filename,args.valorfilename,args.mach3selection,args.tag,args.sample)


