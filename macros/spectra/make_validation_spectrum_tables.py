#!/usr/bin/python

import ROOT
import math
import array
import numpy as np
import sys
import os

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re"]
species   = ["numu","nue","numubar","nuebar","signue"] #,"signuebar"]
processes = ["CCQE", "CC1pi", "CCcoherent", "CCother", "NC1pi0" , "NC1piPM", "NCcoherent", "NCother"]
processes_used = ["CCQE", "CC1pi", "CCcoherent", "CCother" , "NC1piPM", "NCother"]

col1_width = 12
coln_width = 20
to3dp = "{0:.3f}"

modes_to_species_and_process = dict()

modes_to_species_and_process[conv.kSKmNuMuCCQE]         = ("numu", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuCC1pi]        = ("numu", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuCCcoh]        = ("numu", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuCCoth]        = ("numu", "CCother")
modes_to_species_and_process[conv.kSKmNuMuNC1piC]       = ("numu", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuNC1pi0]       = ("numu", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuNCcoh]        = ("numu", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuNCoth_mu]     = ("numu", "NCother") #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuMuNCoth_e]      = ("numu","NCother") nu_e appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCQE]      = ("oscillated numu CCQE
#~ modes_to_species_and_process[conv.kSKmOscNuMuCC1pi]     = ("oscillated numu CC 1pi
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCcoh]     = ("oscillated numu coh
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCoth]     = ("oscillated numu CC other")
modes_to_species_and_process[conv.kSKmNuMuBarCCQE]      = ("numubar", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuBarCC1pi]     = ("numubar", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuBarCCcoh]     = ("numubar", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarCCoth]     = ("numubar", "CCother")
modes_to_species_and_process[conv.kSKmNuMuBarNC1piC]    = ("numubar", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuBarNC1pi0]    = ("numubar", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuBarNCcoh]     = ("numubar", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarNCoth_mu]  = ("numubar", "NCother") #nu_mu appearance analysis definition"
#~ modes_to_species_and_process[conv.kSKmNuMuBarNCoth_e]   = ("numubar", "NCother") nu_e appearance analysis definition")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCQE]   = ("oscillated numubar CCQE")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCC1pi]  = ("oscillated numubar CC 1pi")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCcoh]  = ("oscillated numubar CC coh")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCoth]  = ("oscillated numubar CC other")
modes_to_species_and_process[conv.kSKmNuECCQE]          = ("nue","CCQE")
modes_to_species_and_process[conv.kSKmNuECC1pi]         = ("nue","CC1pi")
modes_to_species_and_process[conv.kSKmNuECCcoh]         = ("nue","CCcoherent")
modes_to_species_and_process[conv.kSKmNuECCoth]         = ("nue","CCother")
modes_to_species_and_process[conv.kSKmNuENC1piC]        = ("nue","NC1piPM")
modes_to_species_and_process[conv.kSKmNuENC1pi0]        = ("nue","NC1pi0")
modes_to_species_and_process[conv.kSKmNuENCcoh]         = ("nue","NCcoherent")
modes_to_species_and_process[conv.kSKmNuENCoth_mu]      = ("nue","NCother") # nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuENCoth_e]       = ("nue","NCother") #nu_e appearance analysis definition")
modes_to_species_and_process[conv.kSKmOscNuEBarCCQE]    = ("nuebar","CCQE")
modes_to_species_and_process[conv.kSKmOscNuEBarCC1pi]   = ("nuebar","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuEBarCCcoh]   = ("nuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuEBarCCoth]   = ("nuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarCCQE]       = ("nuebar","CCQE")
modes_to_species_and_process[conv.kSKmNuEBarCC1pi]      = ("nuebar","CC1pi")
modes_to_species_and_process[conv.kSKmNuEBarCCcoh]      = ("nuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarCCoth]      = ("nuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarNC1piC]     = ("nuebar","NC1piPM")
modes_to_species_and_process[conv.kSKmNuEBarNC1pi0]     = ("nuebar","NC1pi0")
modes_to_species_and_process[conv.kSKmNuEBarNCcoh]      = ("nuebar","NCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarNCoth_mu]   = ("nuebar","NCother") #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuEBarNCoth_e]    = ("nuebar","NCother") #nu_e appearance analysis definition
modes_to_species_and_process[conv.kSKmOscNuECCQE]       = ("signue","CCQE")
modes_to_species_and_process[conv.kSKmOscNuECC1pi]      = ("signue","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuECCcoh]      = ("signue","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuECCoth]      = ("signue","CCother")

fancy_species = dict()
fancy_species ["numu"      ] = "$\\nu_\\mu$"
fancy_species ["nue"       ] = "$\\nu_e$"
fancy_species ["numubar"   ] = "$\\bar{\\nu_\\mu}$"
fancy_species ["nuebar"    ] = "$\\bar{\\nu_e}$"
fancy_species ["signue"    ] = "Osc. $\\nu_e$"
fancy_species ["signuebar" ] = "Osc. $\\bar{\\nu_e}$"

#______________________________________________________________________________________________________________________________

def main(config_file,mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp,hierarchy,angle_scheme,output_file):
  
  opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
  opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
  opt.SetDefaultParameterValues()
  
  phyp = valor.PhysHypothesis("3f_osc")
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  
  phyp.SetCurrValue(conv.kPhypMixing12,  mixing12);
  phyp.SetCurrValue(conv.kPhypMixing13,  mixing13);
  phyp.SetCurrValue(conv.kPhypMixing23,  mixing23);
  phyp.SetCurrValue(conv.kPhypCPV    ,   deltacp);
  phyp.SetCurrValue(conv.kPhypDm12sq,    dm12sq);
  phyp.SetCurrValue(conv.kPhypDm23sq,    dm23sq);
  
  #~ print mixing12
  #~ print mixing13
  #~ print mixing23
  #~ print deltacp 
  #~ print dm12sq  
  #~ print dm23sq  
  
  skexp.Generate(phyp)
  skmc = skexp.GetSuperkMC()
  
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt,a4paper]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\begin{document}
""")
  
  for sample_name in samples:
    
    event_counts_per_species = dict()
    for proc in processes:
      event_counts_per_species[proc] = dict()
  
    s = util.SKSampleFromString(sample_name)
    
    for mode,(spec,proc) in modes_to_species_and_process.iteritems():
      
      spectrum = skexp.CurrPredictedErecoSpectrum(s,mode)
      if not spectrum:
        #~ print "Bad spectrum: ","\t",spec,"\t",proc
        event_counts_per_species[proc][spec] = "Unused"
        continue
        
      event_counts_per_species[proc][spec] = to3dp.format(spectrum.Integral())
      
      #~ print spec,"\t",proc,"\t",event_counts_per_species[proc][spec]
      
    
    tex_output.write("\n\\section*{"+sample_name+"}\n")
    tex_output.write("\\begin{tabular} { | l ||"+("r | "*len(species))+"}\n\\hline\n")
    tex_output.write(" "+"".ljust(col1_width)+" & "+ " & ".join([fancy_species[x].ljust(coln_width) for x in species]) + " \\\\ \\hline\\hline\n")
    
    for proc in processes_used:
      ev = event_counts_per_species[proc]
      tex_output.write(" "+proc.ljust(col1_width)+" & ")
      ev_counts = []
      for spec in species:
        try:
          ev_counts.append(ev[spec].ljust(coln_width))
        except KeyError:
          ev_counts.append("N/A".ljust(coln_width))
          
      tex_output.write(" & ".join(ev_counts) + " \\\\ \\hline\n")
      
    tex_output.write("\n\\end{tabular}\n")
    total_spectrum  = skexp.CurrPredictedErecoSpectrum(s)
    total_events = total_spectrum.Integral()
    
    tex_output.write("\nTotal: "+to3dp.format(total_events)+"\n\n")
    tex_output.write("POT: "+to3dp.format(skmc.AnalyzedSamplePOT(s)*1E-20)+"E20\n\n")
    
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file

if __name__=="__main__":
  #~ main("analysis_configuration_run1to5.xml",0.,0.,0.,7.5E-5,2.4E-3,0,"normal","sinsqtheta","validation.tex")
  main("analysis_configuration_run1to5.xml", 0.846,0.093,0.999,7.53E-5,2.44E-3,0.,"normal","sinsq2theta","validation.tex")
  
