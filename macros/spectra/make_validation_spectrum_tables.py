#!/usr/bin/env python

import ROOT
import math
import array
import sys
import os
import subprocess

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

species   = ["numu","nue","numubar","nuebar","signue", "signuebar"]
#processes      = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother"]
#processes_used = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother"]
processes      = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother", "NC1gamma"]
processes_used = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother", "NC1gamma"]


#samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re"]
#samples   = ["MultiRe"]

samples   = []
if t2kopt.fFit_1Rmu     : samples.append("1Rmu")
if t2kopt.fFit_1Re      : samples.append("1Re")
if t2kopt.fFit_RHC_1Rmu : samples.append("RHC1Rmu")
if t2kopt.fFit_RHC_1Re  : samples.append("RHC1Re")
if t2kopt.fFit_MultiRmu : samples.append("MultiRmu")
if t2kopt.fFit_MultiRe  : samples.append("MultiRe")

col1_width = 12
coln_width = 20
to3dp = "{0:.3f}"
to5dp = "{0:.5f}"
to7dp = "{0:.5f}"

modes_to_species_and_process = dict()

modes_to_species_and_process[conv.kSKmNuMuCCQE]         = ("numu", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuCC1pi]        = ("numu", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuCCcoh]        = ("numu", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuCC2p2h]       = ("numu", "2p2h")
modes_to_species_and_process[conv.kSKmNuMuCCoth]        = ("numu", "CCother")
modes_to_species_and_process[conv.kSKmNuMuNC1piC]       = ("numu", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuNC1pi0]       = ("numu", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuNCcoh]        = ("numu", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuNCoth]        = ("numu", "NCother")  # unified defini,.qtion
modes_to_species_and_process[conv.kSKmNuMuNC1Gamma]     = ("numu", "NC1gamma")  # unified definition
#~ modes_to_species_and_process[conv.kSKmNuMuNCoth_mu]     = ("numu", "NCother")  #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuMuNCoth_e]      = ("numu", "NCother")  #nu_e appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCQE]      = ("oscillated numu CCQE
#~ modes_to_species_and_process[conv.kSKmOscNuMuCC1pi]     = ("oscillated numu CC 1pi
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCcoh]     = ("oscillated numu coh
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCoth]     = ("oscillated numu CC other")
modes_to_species_and_process[conv.kSKmNuMuBarCCQE]      = ("numubar", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuBarCC1pi]     = ("numubar", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuBarCCcoh]     = ("numubar", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarCC2p2h]     = ("numubar", "2p2h")
modes_to_species_and_process[conv.kSKmNuMuBarCCoth]     = ("numubar", "CCother")
modes_to_species_and_process[conv.kSKmNuMuBarNC1piC]    = ("numubar", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuBarNC1pi0]    = ("numubar", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuBarNCcoh]     = ("numubar", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarNCoth]     = ("numubar", "NCother") #unified definition
modes_to_species_and_process[conv.kSKmNuMuBarNC1Gamma]  = ("numubar", "NC1gamma") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuMuBarNCoth_mu]  = ("numubar", "NCother")  #nu_mu appearance analysis definition"
#~ modes_to_species_and_process[conv.kSKmNuMuBarNCoth_e]   = ("numubar", "NCother")  #nu_e appearance analysis definition"
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCQE]   = ("oscillated numubar CCQE")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCC1pi]  = ("oscillated numubar CC 1pi")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCcoh]  = ("oscillated numubar CC coh")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCoth]  = ("oscillated numubar CC other")
modes_to_species_and_process[conv.kSKmNuECCQE]          = ("nue","CCQE")
modes_to_species_and_process[conv.kSKmNuECC1pi]         = ("nue","CC1pi")
modes_to_species_and_process[conv.kSKmNuECCcoh]         = ("nue","CCcoherent")
modes_to_species_and_process[conv.kSKmNuECC2p2h]         = ("nue","2p2h")
modes_to_species_and_process[conv.kSKmNuECCoth]         = ("nue","CCother")
modes_to_species_and_process[conv.kSKmNuENC1piC]        = ("nue","NC1piPM")
modes_to_species_and_process[conv.kSKmNuENC1pi0]        = ("nue","NC1pi0")
modes_to_species_and_process[conv.kSKmNuENCcoh]         = ("nue","NCcoherent")
modes_to_species_and_process[conv.kSKmNuENCoth]         = ("nue","NCother") #unified definition
modes_to_species_and_process[conv.kSKmNuENC1Gamma]      = ("nue","NC1gamma") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuENCoth_mu]      = ("nue","NCother") # nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuENCoth_e]       = ("nue","NCother") #nu_e appearance analysis definition"
modes_to_species_and_process[conv.kSKmOscNuEBarCCQE]    = ("signuebar","CCQE")
modes_to_species_and_process[conv.kSKmOscNuEBarCC1pi]   = ("signuebar","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuEBarCCcoh]   = ("signuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuEBarCC2p2h]  = ("signuebar","2p2h")
modes_to_species_and_process[conv.kSKmOscNuEBarCCoth]   = ("signuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarCCQE]       = ("nuebar","CCQE")
modes_to_species_and_process[conv.kSKmNuEBarCC1pi]      = ("nuebar","CC1pi")
modes_to_species_and_process[conv.kSKmNuEBarCCcoh]      = ("nuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarCC2p2h]     = ("nuebar","2p2h")
modes_to_species_and_process[conv.kSKmNuEBarCCoth]      = ("nuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarNC1piC]     = ("nuebar","NC1piPM")
modes_to_species_and_process[conv.kSKmNuEBarNC1pi0]     = ("nuebar","NC1pi0")
modes_to_species_and_process[conv.kSKmNuEBarNCcoh]      = ("nuebar","NCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarNCoth]      = ("nuebar","NCother") #unified definition
modes_to_species_and_process[conv.kSKmNuEBarNC1Gamma]   = ("nuebar","NC1gamma") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuEBarNCoth_mu]   = ("nuebar","NCother") #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuEBarNCoth_e]    = ("nuebar","NCother") #nu_e appearance analysis definition
modes_to_species_and_process[conv.kSKmOscNuECCQE]       = ("signue","CCQE")
modes_to_species_and_process[conv.kSKmOscNuECC1pi]      = ("signue","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuECCcoh]      = ("signue","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuECC2p2h]     = ("signue","2p2h")
modes_to_species_and_process[conv.kSKmOscNuECCoth]      = ("signue","CCother")

fancy_species = dict()
fancy_species ["numu"      ] = "$\\nu_\\mu$"
fancy_species ["nue"       ] = "$\\nu_e$"
fancy_species ["numubar"   ] = "$\\bar{\\nu}_{\\mu}$"
fancy_species ["nuebar"    ] = "$\\bar{\\nu}_{e}$"
fancy_species ["signue"    ] = "Osc. $\\nu_e$"
fancy_species ["signuebar" ] = "Osc. $\\bar{\\nu}_{e}$"
fancy_species ["total"     ] = "Total"

signal_for_sample = dict()
signal_for_sample["1Rmu"   ] = "numu"
signal_for_sample["1Re"    ] = "signue"
signal_for_sample["RHC1Rmu"] = "numubar"
signal_for_sample["RHC1Re" ] = "signuebar"
signal_for_sample["MultiRe"] = "signue"

#______________________________________________________________________________________________________________________________

def main(config_file,mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp,hierarchy,angle_scheme,octant,dmsqscheme,output_file,nobanff,tweakedbanff,label,
         useBar = False, useNormBar = True,  mixing12Bar = 0, mixing13Bar = 0, mixing23Bar = 0,dm12sqBar = 0, dm23sqBar = 0, deltacpBar = 0 , betaBar = 0):
  
  opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
  opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
  opt.gOscOctantChoice = util.OctantFromString(octant)
  opt.gOscDmsqScheme = util.DmsqSchemeFromString(dmsqscheme)
  
  opt.SetDefaultParameterValues()
  
  if useBar:
    if useNormBar:
      phyp = valor.PhysHypothesis("3f_oscbar_normbar")
    else:
      phyp = valor.PhysHypothesis("3f_oscbar")
  else:
    if useNormBar:
      phyp = valor.PhysHypothesis("3f_osc_normbar")
    else:
      phyp = valor.PhysHypothesis("3f_osc")

  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  skmc = skexp.GetSuperkMC()
  
  if mixing12 < 0.:
    mixing12 = opt.GetGlobalBestFit(conv.kPhypMixing12)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing12),": ",mixing12
  if mixing13 < 0.:
    #~ mixing13 = opt.GetGlobalBestFit(conv.kPhypMixing13)
    if hierarchy=="normal": mixing13 = conv.kT2KBestFitTheta13_NH
    else:                   mixing13 = conv.kT2KBestFitTheta13_IH
    print "Using default ("+hierarchy+") value for ",util.ParamAsString(conv.kPhypMixing13),": ",mixing13
  if mixing23 < 0.:
    mixing23 = opt.GetGlobalBestFit(conv.kPhypMixing23)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing23),": ",mixing23
  if dm12sq < 0.:
    dm12sq = opt.GetGlobalBestFit(conv.kPhypDm12sq)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm12sq),": ",dm12sq
  if dm23sq < 0.:
    dm23sq = opt.GetGlobalBestFit(conv.kPhypDm23sq)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm23sq),": ",dm23sq
  if deltacp <= -9:
    if hierarchy=="normal": deltacp = conv.kT2KBestFitDeltaCP_NH
    else:                   deltacp = conv.kT2KBestFitDeltaCP_IH
    print "Using default ("+hierarchy+") value for ",util.ParamAsString(conv.kPhypCPV),": ",deltacp
    
  # barred parameters
  if mixing12Bar < 0.:
    mixing12Bar = opt.GetGlobalBestFit(conv.kPhypMixing12Bar)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing12Bar),": ",mixing12Bar
  if mixing13Bar < 0.:
    #~ mixing13Bar = opt.GetGlobalBestFit(conv.kPhypMixing13Bar)
    if hierarchy=="normal": mixing13Bar = conv.kT2KBestFitTheta13_NH
    else:                   mixing13Bar = conv.kT2KBestFitTheta13_IH
    print "Using default ("+hierarchy+") value for ",util.ParamAsString(conv.kPhypMixing13Bar),": ",mixing13Bar
  if mixing23Bar < 0.:
    mixing23Bar = opt.GetGlobalBestFit(conv.kPhypMixing23Bar)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing23Bar),": ",mixing23Bar
  if dm12sqBar < 0.:
    dm12sqBar = opt.GetGlobalBestFit(conv.kPhypDm12sqBar)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm12sqBar),": ",dm12sqBar
  if dm23sqBar < 0.:
    dm23sqBar = opt.GetGlobalBestFit(conv.kPhypDm23sqBar)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm23sqBar),": ",dm23sqBar
  if deltacpBar <= -9:
    if hierarchy=="normal": deltacpBar = conv.kT2KBestFitDeltaCP_NH
    else:                   deltacpBar = conv.kT2KBestFitDeltaCP_IH
    print "Using default ("+hierarchy+") value for ",util.ParamAsString(conv.kPhypCPVBar),": ",deltacp

  if betaBar < 0.:
    betaBar = opt.GetGlobalBestFit(conv.kPhypBetaBar)
    print "Using default value for ",util.ParamAsString(conv.kPhypBetaBar),": ",betaBar
   

  phyp.SetCurrValue(conv.kPhypMixing12,  mixing12);
  phyp.SetCurrValue(conv.kPhypMixing13,  mixing13);
  phyp.SetCurrValue(conv.kPhypMixing23,  mixing23);
  phyp.SetCurrValue(conv.kPhypCPV     ,   deltacp);
  phyp.SetCurrValue(conv.kPhypDm12sq  ,    dm12sq);
  phyp.SetCurrValue(conv.kPhypDm23sq  ,    dm23sq);
  
  print util.ParamAsString(conv.kPhypMixing12).ljust(25), mixing12
  print util.ParamAsString(conv.kPhypMixing13).ljust(25), mixing13
  print util.ParamAsString(conv.kPhypMixing23).ljust(25), mixing23
  print util.ParamAsString(conv.kPhypCPV     ).ljust(25), deltacp
  print util.ParamAsString(conv.kPhypDm12sq  ).ljust(25), dm12sq 
  print util.ParamAsString(conv.kPhypDm23sq  ).ljust(25), dm23sq 
  
  if useBar:
    phyp.SetCurrValue(conv.kPhypMixing12Bar ,  mixing12Bar);
    phyp.SetCurrValue(conv.kPhypMixing13Bar ,  mixing13Bar);
    phyp.SetCurrValue(conv.kPhypMixing23Bar ,  mixing23Bar);
    phyp.SetCurrValue(conv.kPhypCPVBar      ,   deltacpBar);
    phyp.SetCurrValue(conv.kPhypDm12sqBar   ,    dm12sqBar);
    phyp.SetCurrValue(conv.kPhypDm23sqBar   ,    dm23sqBar);
  
    print
    print util.ParamAsString(conv.kPhypMixing12Bar).ljust(25), mixing12Bar
    print util.ParamAsString(conv.kPhypMixing13Bar).ljust(25), mixing13Bar
    print util.ParamAsString(conv.kPhypMixing23Bar).ljust(25), mixing23Bar
    print util.ParamAsString(conv.kPhypCPVBar     ).ljust(25), deltacpBar
    print util.ParamAsString(conv.kPhypDm12sqBar  ).ljust(25), dm12sqBar 
    print util.ParamAsString(conv.kPhypDm23sqBar  ).ljust(25), dm23sqBar 
  
  if useNormBar:
    phyp.SetCurrValue(conv.kPhypBetaBar     ,      betaBar);
      
    print util.ParamAsString(conv.kPhypBetaBar    ).ljust(25), betaBar

  print "\n"
  
  # Disable BANFF tuning if requested
  if(nobanff):
    print "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..."
    skexp.GetSuperkMC().FillTemplateArrays(False)
    # New interface to Correlated Systematics -- see BANFF.h
    skexp.GetBANFF().Tweak(0,0,"untuned")
     
    skexp.GetBANFF().PrintSystematicVariations()
    
    #~ skexp.GetNuXSec().PrintSystematicVariations()

  # Special tweak to compare with prefit values
  # BANFF XSec parameters for convenience
  # 50 (default: 1)  BANFF; Norm; 2p2h                  kNuisBANFF2p2h                
  # 51 (default: 1)  BANFF; CA5 RES                     kNuisBANFFCA5                
  # 52 (default: 1)  BANFF; BgRES Isospin 1/2           kNuisBANFFBgRES              
  # 53 (default: 1)  BANFF; Ma QE                       kNuisBANFFMAQE               
  # 54 (default: 1)  BANFF; Ma RES                      kNuisBANFFMARES              
  # 55 (default: 1)  BANFF; Second current class axial  kNuisBANFFSCCA               
  # 56 (default: 1)  BANFF; Second current class vector kNuisBANFFSCCV               
  # 57 (default: 1)  BANFF; Fermi Momentum              kNuisBANFFFermiMom           
  # 58 (default: 0)  BANFF; Shape; CC Oth               kNuisBANFFShapeCCoth         
  # 59 (default: 1)  BANFF; Binding energy              kNuisBANFFBindE              
  # 60 (default: 1)  BANFF; Norm; CC Coh                kNuisBANFFNormCCcoh          
  # 61 (default: 1)  BANFF; Norm; NC Coh                kNuisBANFFNormNCcoh          
  # 62 (default: 1)  BANFF; Norm; NC Oth                kNuisBANFFNormNCoth          
  # 63 (default: 1)  BANFF; Norm; NuE To NuMu           kNuisBANFFNormNuEToNuMu      
  # 64 (default: 1)  BANFF; Norm; NC 1 Gamma            kNuisBANFFNormNC1Gamma       
  # 65 (default: 1)  BANFF; Norm; NuEBar To NuMuBar     kNuisBANFFNormNuEBarToNuMuBar
  # 66 (default: 1)  BANFF; Norm; 2p2hBar               kNuisBANFF2p2hBar             

  if(tweakedbanff):
    print "- Tweaking BANFF tuning to match ptheta analysis"
#    skexp.GetSuperkMC().FillTemplateArrays(False)

    v = skexp.GetBANFF().CurrentValues()
    v_untuned = skexp.GetBANFF().UntunedValues()
    v=v_untuned
    # for i in range(50,67):
    #   #v[i] = v_untuned[i];
    #   print "Parameter: ", util.ParamAsString(skexp.GetBANFF().Idx2Param(i))
    #   if ( skexp.GetBANFF().Idx2Param(i) == conv.kNuisBANFFMAQE ):
    #     print "Param at :", i, util.ParamAsString(skexp.GetBANFF().Idx2Param(i)), "Postfit value set ", v[i]
    #     continue;
    #   elif (util.IsNormXSecParam(skexp.GetBANFF().Idx2Param(i))==False):
    #     print "\t", util.ParamAsString(skexp.GetBANFF().Idx2Param(i)), " set to untuned values"
    #     v[i] = v_untuned[i];

    #v[50] =     # BANFF; Norm; 2p2h                        1     | .  #50 MA QE, type:3                       
    v[51] = 0.950495049505  # BANFF; CA5 RES                          0.78  | .  #51 pF, type:3                          
    v[52] = 0.738461538462  # BANFF; BgRES Isospin 1/2                1     |   #52 2p2h O norm, type:2                  
    v[53] = 1.0  # BANFF; Ma QE                            0.93  | .  #53 Eb, type:3                          
    v[54] = 1.12631578947  # BANFF; Ma RES                           0.88  |   #54 2p2h nubar, type:2                   
    v[55] = 1.0  # BANFF; Second current class axial       1     | .  #55 CA5_RES, type:3                     
    v[56] = 1.0  # BANFF; Second current class vector      1     | .  #56 MANFF_RES, type:3                   
    v[57] = 1.0  # BANFF; Fermi Momentum                   1.04  | .  #57 BG_RES, type:3                      
    # v[58] = 0.0  # BANFF; Shape; CC Oth                    -0.03 |  .  #58 nue/numu xsec, type:4               
    # v[59] =    # BANFF; Norm; CC Coh                     1     | .  #60 CCother Shape/DIS Multi-PI, type:3  
    # v[60] =    # BANFF; Norm; NC Coh                     1     |   #61 CC coherent norm, type:2            
    # v[61] =    # BANFF; Norm; NC Oth                     1     |   #62 NC coherent norm, type:2            
    # v[62] =    # BANFF; Norm; NuE To NuMu                1     |   #63 NC gamma norm, type:2               
    # v[63] =    # BANFF; Norm; NC 1 Gamma                 1     |   #64 NC other norm, type:2               
    # v[64] =    # BANFF; Norm; NuEBar To NuMuBar          1     |   #65 RPA, type:3                         
    # v[65] =    # BANFF; Norm; 2p2hBar                     0.57  |
    v[66] =  0.59   # BANFF; BeRPA A
    v[67] =  1.05   # BANFF; BeRPA B
    v[68] =  1.13   # BANFF; BeRPA D
    v[69] =  0.88   # BANFF; BeRPA E
    v[70] =  1.2    # BANFF; BeRPA U
    # v[72] =    # BANFF; Shape; 2p2h nu
    # v[73] =    # BANFF; Shape; 2p2h nubar

    skexp.GetBANFF().Tweak(1,v,"")
    #skexp.GetBANFF().PrintSystematicVariations()

    print "Special tweak for comparison with prefit systematic variations "
  skexp.GetBANFF().PrintSystematicVariations()

  # Make the dataset
  print "Making the data set"
  skexp.Generate(phyp)
  
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt,a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\begin{document}
""")

  # Write osc params
  params = [conv.kPhypMixing12, conv.kPhypMixing13, conv.kPhypMixing23, conv.kPhypDm12sq , conv.kPhypDm23sq, conv.kPhypCPV]
  
  bar_params = [conv.kPhypMixing12Bar, conv.kPhypMixing13Bar, conv.kPhypMixing23Bar, conv.kPhypDm12sqBar , conv.kPhypDm23sqBar, conv.kPhypCPVBar, conv.kPhypBetaBar]
  
  normbar_params = [conv.kPhypBetaBar]
  
  tex_output.write("\n\\section*{Job details}\n\n")
  
  tex_output.write("\n\\begin{minipage}[t]{0.5\\textwidth}")
  tex_output.write("\n\\begin{itemize}\n")
  for param in params:
    tex_output.write("\\item[] "+util.ParamAsStringTex(param)+" = "+str(phyp.CurrValue(param))+"\n")
  
  tex_output.write("\\end{itemize}\n\\end{minipage}")
  
  if useBar:
    tex_output.write("\n\\begin{minipage}[t]{0.5\\textwidth}")
    tex_output.write("\n\\begin{itemize}\n")
    for param in bar_params:
      tex_output.write("\\item[] "+util.ParamAsStringTex(param)+" = "+str(phyp.CurrValue(param))+"\n")
    
    tex_output.write("\\end{itemize}\n\\end{minipage}")

  if useNormBar:
    tex_output.write("\n\\begin{minipage}[t]{0.5\\textwidth}")
    tex_output.write("\n\\begin{itemize}\n")
    for param in normbar_params:
      tex_output.write("\\item[] "+util.ParamAsStringTex(param)+" = "+str(phyp.CurrValue(param))+"\n")
    
    tex_output.write("\\end{itemize}\n\\end{minipage}")
    
  tex_output.write("\n\\\\\n\n")
  
  banff_status = ["Enabled","Disabled"][int(nobanff)]

  #banff_tweaked = ["Enabled","Disabled"][int(tweakedbanff)]
  
  tex_output.write("\n\nBANFF tunings: "+ banff_status +"\n\n")
  # tex_output.write("\n\nSpecial preBANFF tunings: "+ banff_tweaked +"\n\n")
  tex_output.write("Hierarchy: "+hierarchy.capitalize()+"\n\n")
  if label: tex_output.write(label+"\n\n")
  
  for isample,sample_name in enumerate(samples):
    
    nev_tot = 0.
    
    event_counts_per_species = dict()
    total_event_counts_per_species = dict()
    total_event_counts_per_process = dict()
    
    for proc in processes:
      event_counts_per_species[proc] = dict()
      total_event_counts_per_process[proc] = 0.
      
    for spec in species:
      total_event_counts_per_species[spec] = 0.
  
    s = util.SKSampleFromString(sample_name)
    signal_events     = 0.
    background_events = 0.
    
    # Count events
    for mode,(spec,proc) in modes_to_species_and_process.iteritems():
      
      #~ if util.Is1RESample (s) and (mode in [conv.kSKmNuMuNCoth_mu, conv.kSKmNuMuBarNCoth_mu, conv.kSKmNuENCoth_mu,  conv.kSKmNuEBarNCoth_mu ]): continue  
      #~ if util.Is1RMuSample(s) and (mode in [conv.kSKmNuMuNCoth_e , conv.kSKmNuMuBarNCoth_e , conv.kSKmNuENCoth_e ,  conv.kSKmNuEBarNCoth_e  ]): continue  
      
      #spectrum = skexp.CurrPredictedErecoSpectrum(s,mode)
      spectrum = skexp.CurrPredictedErecoVsQ2Spectrum(s,mode)
      if not spectrum:
        #~ print "Bad spectrum: ","\t",spec,"\t",proc
        event_counts_per_species[proc][spec] = "Unused"
        continue
      
      nev = spectrum.Integral()
      nev_tot += nev
      total_event_counts_per_species[spec] += nev
      total_event_counts_per_process[proc] += nev
      event_counts_per_species[proc][spec] = to5dp.format(nev)
      
      if (signal_for_sample[sample_name]==spec): signal_events+=nev
      else: background_events+=nev
      
      #~ print spec,"\t",proc,"\t",event_counts_per_species[proc][spec]
      
    
    # Start writing output
    tex_output.write("\n\\section*{"+sample_name+"}\n")
    tex_output.write("\\begin{tabular} { | l ||"+("r | "*len(species))+"| r | }\n\\hline\n")
    tex_output.write(" "+"".ljust(col1_width)+" & "+ " & ".join([fancy_species[x].ljust(coln_width) for x in species+["total"]]) + " \\\\ \\hline\\hline\n")
    
    for proc in processes_used:
      ev = event_counts_per_species[proc]
      if not total_event_counts_per_process[proc]: continue
      tex_output.write(" "+proc.ljust(col1_width)+" & ")
      ev_counts = []
      for spec in species:
        try:
          ev_counts.append(ev[spec].ljust(coln_width))
        except KeyError:
          ev_counts.append("N/A".ljust(coln_width))
      
      ev_counts.append(to5dp.format(total_event_counts_per_process[proc]).ljust(coln_width))
          
      tex_output.write(" & ".join(ev_counts) + " \\\\ \\hline\n")
    
    # Write totals
    tex_output.write("\\hline\n")
    tex_output.write(" "+"Total".ljust(col1_width)+" & ")
    string_totals = [to5dp.format(total_event_counts_per_species[spec]).ljust(coln_width) for spec in species]
    string_totals.append(to7dp.format(nev_tot).ljust(coln_width))
    tex_output.write(" & ".join(string_totals) + " \\\\ \\hline\n")
    
    tex_output.write("\n\\end{tabular}\n")
    #total_spectrum  = skexp.CurrPredictedErecoSpectrum(s)
    total_spectrum  = skexp.CurrPredictedErecoVsQ2Spectrum(s)
    total_events = total_spectrum.Integral()
    
    print util.SampleAsString(s),total_events
    
    tex_output.write("\nTotal: "           +to5dp.format(total_events)     +"    ("+to5dp.format(nev_tot)+")\n\n")
    tex_output.write("\nTotal signal: "    +to5dp.format(signal_events)    +"\n\n")
    tex_output.write("\nTotal background: "+to5dp.format(background_events)+"\n\n")
    
    tex_output.write("POT: "+to5dp.format(skmc.AnalyzedSamplePOT(s)*1E-20)+"E20\n\n")
    
    if (isample%2): tex_output.write("\n\n\\clearpage\n\n")
    
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file    

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")

  # # T2K+reactors (sinsqtheta) - Normal Hierarchy
  parser.add_option('--mixing12', dest='mixing12',default=0.   ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  parser.add_option('--mixing13', dest='mixing13',default=0.  ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  parser.add_option('--mixing23', dest='mixing23',default=0.   ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.53E-5 ,type=float, help='Value of small mass splitting')
  parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.509E-3 ,type=float, help='Value of large mass splitting')
  parser.add_option('--deltacp' , dest='deltacp' ,default=0.    ,type=float, help='Value of delta CP')
  parser.add_option('--mixing12bar', dest='mixing12Bar',default=0.   ,type=float, help='Barred value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  parser.add_option('--mixing13bar', dest='mixing13Bar',default=0.  ,type=float, help='Barred value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  parser.add_option('--mixing23bar', dest='mixing23Bar',default=0.   ,type=float, help='Barred value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  parser.add_option('--dm12sqbar'  , dest='dm12sqBar'  ,default=7.53E-5 ,type=float, help='Barred value of small mass splitting')
  parser.add_option('--dm23sqbar'  , dest='dm23sqBar'  ,default=2.509E-3 ,type=float, help='Barred value of large mass splitting')
  parser.add_option('--deltacpbar' , dest='deltacpBar' ,default=0.   ,type=float, help='Barred value of delta CP')
  parser.add_option('--betabar'    , dest='betaBar'    ,default=1.      ,type=float, help='Barred value of nuebar appearance normalization ')
  parser.add_option('--hierarchy'  , dest="hierarchy"  , default="normal", help='Mass hierarchy: "normal" or "inverted"')

  # T2K+reactors best-fit (sinsqtheta) - Asimov data set (A)
  # parser.add_option('--mixing12', dest='mixing12',default=0.304   ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13', dest='mixing13',default=0.0217  ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23', dest='mixing23',default=0.528   ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.53E-5 ,type=float, help='Value of small mass splitting')
  # parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.509E-3 ,type=float, help='Value of large mass splitting')
  # parser.add_option('--deltacp' , dest='deltacp' ,default=-1.601    ,type=float, help='Value of delta CP')
  # parser.add_option('--mixing12bar', dest='mixing12Bar',default=0.304   ,type=float, help='Barred value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13bar', dest='mixing13Bar',default=0.0217  ,type=float, help='Barred value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23bar', dest='mixing23Bar',default=0.528   ,type=float, help='Barred value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sqbar'  , dest='dm12sqBar'  ,default=7.53E-5 ,type=float, help='Barred value of small mass splitting')
  # parser.add_option('--dm23sqbar'  , dest='dm23sqBar'  ,default=2.509E-3 ,type=float, help='Barred value of large mass splitting')
  # parser.add_option('--deltacpbar' , dest='deltacpBar' ,default=-1.601   ,type=float, help='Barred value of delta CP')
  # parser.add_option('--betabar'    , dest='betaBar'    ,default=1.      ,type=float, help='Barred value of nuebar appearance normalization ')
  # parser.add_option('--hierarchy'  , dest="hierarchy"  , default="normal", help='Mass hierarchy: "normal" or "inverted"')

#  # T2K+reactors best-fit (sinsqtheta) - Asimov data set (B)
#  parser.add_option('--mixing12', dest='mixing12',default=0.304   ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
#  parser.add_option('--mixing13', dest='mixing13',default=0.0217  ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
#  parser.add_option('--mixing23', dest='mixing23',default=0.45   ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
#  parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.53E-5 ,type=float, help='Value of small mass splitting')
#  parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.509E-3 ,type=float, help='Value of large mass splitting')
#  parser.add_option('--deltacp' , dest='deltacp' ,default=0.    ,type=float, help='Value of delta CP')
#  parser.add_option('--mixing12bar', dest='mixing12Bar',default=0.304   ,type=float, help='Barred value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
#  parser.add_option('--mixing13bar', dest='mixing13Bar',default=0.0217  ,type=float, help='Barred value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
#  parser.add_option('--mixing23bar', dest='mixing23Bar',default=0.45   ,type=float, help='Barred value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
#  parser.add_option('--dm12sqbar'  , dest='dm12sqBar'  ,default=7.53E-5 ,type=float, help='Barred value of small mass splitting')
#  parser.add_option('--dm23sqbar'  , dest='dm23sqBar'  ,default=2.509E-3 ,type=float, help='Barred value of large mass splitting')
#  parser.add_option('--deltacpbar' , dest='deltacpBar' ,default=0.   ,type=float, help='Barred value of delta CP')
#  parser.add_option('--betabar'    , dest='betaBar'    ,default=1.      ,type=float, help='Barred value of nuebar appearance normalization ')
#  parser.add_option('--hierarchy'  , dest="hierarchy"  , default="normal", help='Mass hierarchy: "normal" or "inverted"')

  # T2K+reactors best-fit (sinsqtheta) - Leila's weights in presentation
  # parser.add_option('--mixing12', dest='mixing12',default=0.304   ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13', dest='mixing13',default=0.0217  ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23', dest='mixing23',default=0.528   ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.53E-5 ,type=float, help='Value of small mass splitting')
  # parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.509E-3 ,type=float, help='Value of large mass splitting')
  # parser.add_option('--deltacp' , dest='deltacp' ,default=0.    ,type=float, help='Value of delta CP')
  # parser.add_option('--mixing12bar', dest='mixing12Bar',default=0.304   ,type=float, help='Barred value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13bar', dest='mixing13Bar',default=0.0217  ,type=float, help='Barred value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23bar', dest='mixing23Bar',default=0.528   ,type=float, help='Barred value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sqbar'  , dest='dm12sqBar'  ,default=7.53E-5 ,type=float, help='Barred value of small mass splitting')
  # parser.add_option('--dm23sqbar'  , dest='dm23sqBar'  ,default=2.509E-3 ,type=float, help='Barred value of large mass splitting')
  # parser.add_option('--deltacpbar' , dest='deltacpBar' ,default=0.   ,type=float, help='Barred value of delta CP')
  # parser.add_option('--betabar'    , dest='betaBar'    ,default=1.      ,type=float, help='Barred value of nuebar appearance normalization ')
  # parser.add_option('--hierarchy'  , dest="hierarchy"  , default="normal", help='Mass hierarchy: "normal" or "inverted"')

   # T2K+reactors best-fit (sinsqtheta) - Konosuke's weights in presentation 
  # parser.add_option('--mixing12', dest='mixing12',default=math.pow(math.sin(0.5*math.asin(math.sqrt(0.8495))),2)   ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13', dest='mixing13',default=math.pow(math.sin(0.5*math.asin(math.sqrt(0.1))),2)  ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23', dest='mixing23',default= math.pow(math.sin(0.5*math.asin(math.sqrt(1))),2)   ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.6E-5 ,type=float, help='Value of small mass splitting')
  # parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.4E-3 ,type=float, help='Value of large mass splitting')
  # parser.add_option('--deltacp' , dest='deltacp' ,default=0    ,type=float, help='Value of delta CP')
  # parser.add_option('--mixing12bar', dest='mixing12Bar',default=0.846   ,type=float, help='Barred value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  # parser.add_option('--mixing13bar', dest='mixing13Bar',default=0.085  ,type=float, help='Barred value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  # parser.add_option('--mixing23bar', dest='mixing23Bar',default=1.   ,type=float, help='Barred value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  # parser.add_option('--dm12sqbar'  , dest='dm12sqBar'  ,default=7.53E-5 ,type=float, help='Barred value of small mass splitting')
  # parser.add_option('--dm23sqbar'  , dest='dm23sqBar'  ,default=2.509E-3 ,type=float, help='Barred value of large mass splitting')
  # parser.add_option('--deltacpbar' , dest='deltacpBar' ,default=-1.601   ,type=float, help='Barred value of delta CP')
  # parser.add_option('--betabar'    , dest='betaBar'    ,default=1.      ,type=float, help='Barred value of nuebar appearance normalization ')
  # parser.add_option('--hierarchy'  , dest="hierarchy"  , default="normal", help='Mass hierarchy: "normal" or "inverted"')

  
  parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta", help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
  parser.add_option('-t','--mixing-param-octant', dest="mixparamoctant", default="two"       , help='Mixing param octant: "one" or "two"')
  parser.add_option('-a','--dmsq-scheme', dest="dmsqscheme", default="auto"       , help='"abs23", "abs31", "auto" or "fogli"')
  
  parser.add_option('-f', '--file', dest='output_file', default="validation_no_osc_special.tex", help="Name for output TeX file.")
  
  
  #
  # NEW CONFIGURATION FILES (4/3/16 BANFF RELEASE)
  #
  parser.add_option('-c', '--config', dest='config', default="analysis_configuration_Run1-7cPOT_theta_banff_2p2h_berpa.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_FSIOnly.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_data_theta_1bin_RFGRPA_NIWG2015c_1p1h_SKDetOnly.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_data_theta_elike_RFGRPA_NIWG2015c_1p1h_ExpSummer2016.xml", help="Configuration file")
  
  #
  # OLD CONFIGURATION FILES
  #
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c_ExpSummer2016.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_1bin_RFGRPA_NIWG2015c.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_1bin_RFGRPA_NIWG2015c_Asimov.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_1bin_NIWG2015c.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_1bin_NIWG2015c_Asimov.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_theta_elike_RFGRPA_NIWG2015c.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_prefit_theta_elike_RFGRPA_NIWG2015c.xml", help="Configuration file")
  #parser.add_option('-c', '--config', dest='config', default="analysis_configuration_test_prefit_theta_1bin_NIWG2015c.xml", help="Configuration file")

  parser.add_option('-l', '--label' , dest='label' , default=None, help="Label to put in file.")
  
  parser.add_option('-d', '--disable-banff', dest='nobanff', action="store_true",default=False,help="Disable the BANFF tune. This returns to flux tuned only, no NIWG tune")

  parser.add_option('-s', '--special-banff-tweak', dest='tweakedbanff', action="store_true",default=False,help="Tweak BANFF prefit parameters to match what is done by ptheta")
  
  parser.add_option('-b', '--use-bar', dest='useBar', action="store_true",default=False,help="Use the 3f_osc_bar hypothesis. Without this, the bar mixing params will be ignored.")

  parser.add_option('-z', '--use-normbar', dest='useNormBar', action="store_true",default=False,help="Use the normbar hypothesis. Without this, the norm bar will be ignored.")

  v = optparse.Values()
  parser.parse_args(values=v) # v contains only the values that were specified, not defaults
  
  o,other = parser.parse_args()
  
  if not o.useBar: # Check we're not trying to specify a barred mixing parameter for without using the bar hypothesis
    bar_par = ["mixing12Bar","mixing13Bar","mixing23Bar","dm12sqBar","dm23sqBar","deltacpBar"]
    if any(i in bar_par for i in v.__dict__.keys()):
      print "If you want to specify barred mixing parameters, you need -b as well"
      exit(-1)

  if not o.useNormBar: # Check we're not trying to specify a barred mixing parameter for without using the bar hypothesis
    normbar_par = ["betaBar"]
    if any(i in normbar_par for i in v.__dict__.keys()):
      print "If you want to specify barred mixing parameters, you need -z as well"
      exit(-1)
  
  main(o.config,o.mixing12,o.mixing13,o.mixing23,o.dm12sq,o.dm23sq,o.deltacp,o.hierarchy,o.mixparamchoice,o.mixparamoctant,o.dmsqscheme,o.output_file,o.nobanff,o.tweakedbanff, o.label,
       o.useBar, o.useNormBar, o.mixing12Bar,o.mixing13Bar,o.mixing23Bar,o.dm12sqBar,o.dm23sqBar,o.deltacpBar,o.betaBar)
  
  if other:
    print "Warning - ignored arguments: ",other
