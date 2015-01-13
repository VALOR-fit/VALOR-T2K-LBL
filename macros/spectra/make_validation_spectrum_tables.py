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
  
valor   = ROOT.valor
util    = valor.util
conv    = valor.conv
opt     = valor.Opt.Instance()
t2kopt  = valor.T2KOpt.Instance()
env     = valor.Env.Instance()
rnd     = valor.Rnd.Instance()
arridx  = util.ArrIdx.Instance()

species   = ["numu","nue","numubar","nuebar","signue", "signuebar"]
processes = ["CCQE", "CC1pi", "CCcoherent","MEC", "CCother", "NC1pi0" , "NC1piPM", "NCcoherent", "NCother"]
processes_used = ["CCQE", "CC1pi", "CCcoherent", "MEC","CCother", "NC1pi0" , "NC1piPM", "NCcoherent", "NCother"]
#~ processes_used = ["CCQE", "CC1pi", "CCcoherent", "CCother" , "NC1piPM", "NCother"]

#~ samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re"]
samples   = []
if t2kopt.fFit_1Rmu    : samples.append("1Rmu")
if t2kopt.fFit_1Re     : samples.append("1Re")
if t2kopt.fFit_RHC_1Rmu: samples.append("RHC1Rmu")
if t2kopt.fFit_RHC_1Re : samples.append("RHC1Re")

col1_width = 12
coln_width = 20
to3dp = "{0:.3f}"

modes_to_species_and_process = dict()

modes_to_species_and_process[conv.kSKmNuMuCCQE]         = ("numu", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuCC1pi]        = ("numu", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuCCcoh]        = ("numu", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuCCMEC]        = ("numu", "MEC")
modes_to_species_and_process[conv.kSKmNuMuCCoth]        = ("numu", "CCother")
modes_to_species_and_process[conv.kSKmNuMuNC1piC]       = ("numu", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuNC1pi0]       = ("numu", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuNCcoh]        = ("numu", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuNCoth]        = ("numu", "NCother")  # unified definition
#~ modes_to_species_and_process[conv.kSKmNuMuNCoth_mu]     = ("numu", "NCother")  #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuMuNCoth_e]      = ("numu", "NCother")  #nu_e appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCQE]      = ("oscillated numu CCQE
#~ modes_to_species_and_process[conv.kSKmOscNuMuCC1pi]     = ("oscillated numu CC 1pi
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCcoh]     = ("oscillated numu coh
#~ modes_to_species_and_process[conv.kSKmOscNuMuCCoth]     = ("oscillated numu CC other")
modes_to_species_and_process[conv.kSKmNuMuBarCCQE]      = ("numubar", "CCQE")
modes_to_species_and_process[conv.kSKmNuMuBarCC1pi]     = ("numubar", "CC1pi")
modes_to_species_and_process[conv.kSKmNuMuBarCCcoh]     = ("numubar", "CCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarCCMEC]     = ("numubar", "MEC")
modes_to_species_and_process[conv.kSKmNuMuBarCCoth]     = ("numubar", "CCother")
modes_to_species_and_process[conv.kSKmNuMuBarNC1piC]    = ("numubar", "NC1piPM")
modes_to_species_and_process[conv.kSKmNuMuBarNC1pi0]    = ("numubar", "NC1pi0")
modes_to_species_and_process[conv.kSKmNuMuBarNCcoh]     = ("numubar", "NCcoherent")
modes_to_species_and_process[conv.kSKmNuMuBarNCoth]     = ("numubar", "NCother") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuMuBarNCoth_mu]  = ("numubar", "NCother")  #nu_mu appearance analysis definition"
#~ modes_to_species_and_process[conv.kSKmNuMuBarNCoth_e]   = ("numubar", "NCother")  #nu_e appearance analysis definition"
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCQE]   = ("oscillated numubar CCQE")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCC1pi]  = ("oscillated numubar CC 1pi")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCcoh]  = ("oscillated numubar CC coh")
#~ modes_to_species_and_process[conv.kSKmOscNuMuBarCCoth]  = ("oscillated numubar CC other")
modes_to_species_and_process[conv.kSKmNuECCQE]          = ("nue","CCQE")
modes_to_species_and_process[conv.kSKmNuECC1pi]         = ("nue","CC1pi")
modes_to_species_and_process[conv.kSKmNuECCcoh]         = ("nue","CCcoherent")
modes_to_species_and_process[conv.kSKmNuECCMEC]         = ("nue","MEC")
modes_to_species_and_process[conv.kSKmNuECCoth]         = ("nue","CCother")
modes_to_species_and_process[conv.kSKmNuENC1piC]        = ("nue","NC1piPM")
modes_to_species_and_process[conv.kSKmNuENC1pi0]        = ("nue","NC1pi0")
modes_to_species_and_process[conv.kSKmNuENCcoh]         = ("nue","NCcoherent")
modes_to_species_and_process[conv.kSKmNuENCoth]         = ("nue","NCother") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuENCoth_mu]      = ("nue","NCother") # nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuENCoth_e]       = ("nue","NCother") #nu_e appearance analysis definition"
modes_to_species_and_process[conv.kSKmOscNuEBarCCQE]    = ("signuebar","CCQE")
modes_to_species_and_process[conv.kSKmOscNuEBarCC1pi]   = ("signuebar","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuEBarCCcoh]   = ("signuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuEBarCCMEC]   = ("signuebar","MEC")
modes_to_species_and_process[conv.kSKmOscNuEBarCCoth]   = ("signuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarCCQE]       = ("nuebar","CCQE")
modes_to_species_and_process[conv.kSKmNuEBarCC1pi]      = ("nuebar","CC1pi")
modes_to_species_and_process[conv.kSKmNuEBarCCcoh]      = ("nuebar","CCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarCCMEC]      = ("nuebar","MEC")
modes_to_species_and_process[conv.kSKmNuEBarCCoth]      = ("nuebar","CCother")
modes_to_species_and_process[conv.kSKmNuEBarNC1piC]     = ("nuebar","NC1piPM")
modes_to_species_and_process[conv.kSKmNuEBarNC1pi0]     = ("nuebar","NC1pi0")
modes_to_species_and_process[conv.kSKmNuEBarNCcoh]      = ("nuebar","NCcoherent")
modes_to_species_and_process[conv.kSKmNuEBarNCoth]      = ("nuebar","NCother") #unified definition
#~ modes_to_species_and_process[conv.kSKmNuEBarNCoth_mu]   = ("nuebar","NCother") #nu_mu appearance analysis definition
#~ modes_to_species_and_process[conv.kSKmNuEBarNCoth_e]    = ("nuebar","NCother") #nu_e appearance analysis definition
modes_to_species_and_process[conv.kSKmOscNuECCQE]       = ("signue","CCQE")
modes_to_species_and_process[conv.kSKmOscNuECC1pi]      = ("signue","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuECCcoh]      = ("signue","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuECCMEC]      = ("signue","MEC")
modes_to_species_and_process[conv.kSKmOscNuECCoth]      = ("signue","CCother")

fancy_species = dict()
fancy_species ["numu"      ] = "$\\nu_\\mu$"
fancy_species ["nue"       ] = "$\\nu_e$"
fancy_species ["numubar"   ] = "$\\bar{\\nu_\\mu}$"
fancy_species ["nuebar"    ] = "$\\bar{\\nu_e}$"
fancy_species ["signue"    ] = "Osc. $\\nu_e$"
fancy_species ["signuebar" ] = "Osc. $\\bar{\\nu_e}$"
fancy_species ["total"     ] = "Total"

#______________________________________________________________________________________________________________________________

def main(config_file,mixing12,mixing13,mixing23,dm12sq,dm23sq,deltacp,hierarchy,angle_scheme,output_file,nobanff, label):
  
  opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
  opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
  opt.SetDefaultParameterValues()
  
  phyp = valor.PhysHypothesis("3f_osc")
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  skmc = skexp.GetSuperkMC()
  
  if mixing12 < 0.:
    mixing12 = opt.GetGlobalBestFit(conv.kPhypMixing12)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing12),": ",mixing12
  if mixing13 < 0.:
    mixing13 = opt.GetGlobalBestFit(conv.kPhypMixing13)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing13),": ",mixing13
  if mixing23 < 0.:
    mixing23 = opt.GetGlobalBestFit(conv.kPhypMixing23)
    print "Using default value for ",util.ParamAsString(conv.kPhypMixing23),": ",mixing23
  if dm12sq < 0.:
    dm12sq = opt.GetGlobalBestFit(conv.kPhypDm12sq)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm12sq),": ",dm12sq
  if dm23sq < 0.:
    dm23sq = opt.GetGlobalBestFit(conv.kPhypDm23sq)
    print "Using default value for ",util.ParamAsString(conv.kPhypDm23sq),": ",dm23sq
  
  phyp.SetCurrValue(conv.kPhypMixing12,  mixing12);
  phyp.SetCurrValue(conv.kPhypMixing13,  mixing13);
  phyp.SetCurrValue(conv.kPhypMixing23,  mixing23);
  phyp.SetCurrValue(conv.kPhypCPV     ,   deltacp);
  phyp.SetCurrValue(conv.kPhypDm12sq  ,    dm12sq);
  phyp.SetCurrValue(conv.kPhypDm23sq  ,    dm23sq);
  
  #~ print mixing12
  #~ print mixing13
  #~ print mixing23
  #~ print deltacp 
  #~ print dm12sq  
  #~ print dm23sq  
  
  # Disable BANFF tuning if requested
  if(nobanff):
    print "- Performing special systematic parameter variation which cancels the effect of BANFF tuning..."
    skexp.GetBANFF().RestoreNominalMC()
    skexp.GetBANFF().RestoreRHCNominalMC()
    
  #~ skexp.GetNuXSec().PrintSystematicVariations()
  
  # Make the dataset
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
  
  
  tex_output.write("\n\\section*{Job details}\n\n\\begin{itemize}\n")
  for param in params:
    tex_output.write("\\item[] "+util.ParamAsStringTex(param)+" = "+str(phyp.CurrValue(param))+"\n")
  
  tex_output.write("\end{itemize}")
  
  banff_status = ["Enabled","Disabled"][int(nobanff)]
  
  tex_output.write("\n\nBANFF tunings: "+ banff_status +"\n\n")
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
    
    # Count events
    for mode,(spec,proc) in modes_to_species_and_process.iteritems():
      
      #~ if util.Is1RESample (s) and (mode in [conv.kSKmNuMuNCoth_mu, conv.kSKmNuMuBarNCoth_mu, conv.kSKmNuENCoth_mu,  conv.kSKmNuEBarNCoth_mu ]): continue  
      #~ if util.Is1RMuSample(s) and (mode in [conv.kSKmNuMuNCoth_e , conv.kSKmNuMuBarNCoth_e , conv.kSKmNuENCoth_e ,  conv.kSKmNuEBarNCoth_e  ]): continue  
      
      spectrum = skexp.CurrPredictedErecoSpectrum(s,mode)
      if not spectrum:
        #~ print "Bad spectrum: ","\t",spec,"\t",proc
        event_counts_per_species[proc][spec] = "Unused"
        continue
      
      nev = spectrum.Integral()
      nev_tot += nev
      total_event_counts_per_species[spec] += nev
      total_event_counts_per_process[proc] += nev
      event_counts_per_species[proc][spec] = to3dp.format(nev)
      
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
      
      ev_counts.append(to3dp.format(total_event_counts_per_process[proc]).ljust(coln_width))
          
      tex_output.write(" & ".join(ev_counts) + " \\\\ \\hline\n")
    
    # Write totals
    tex_output.write("\\hline\n")
    tex_output.write(" "+"Total".ljust(col1_width)+" & ")
    string_totals = [to3dp.format(total_event_counts_per_species[spec]).ljust(coln_width) for spec in species]
    string_totals.append(to3dp.format(nev_tot).ljust(coln_width))
    tex_output.write(" & ".join(string_totals) + " \\\\ \\hline\n")
    
    tex_output.write("\n\\end{tabular}\n")
    total_spectrum  = skexp.CurrPredictedErecoSpectrum(s)
    total_events = total_spectrum.Integral()
    
    tex_output.write("\nTotal: "+to3dp.format(total_events)+"    ("+to3dp.format(nev_tot)+")\n\n")
    tex_output.write("POT: "+to3dp.format(skmc.AnalyzedSamplePOT(s)*1E-20)+"E20\n\n")
    
    if (isample%2): tex_output.write("\n\n\\clearpage\n\n")
    
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")
  
  parser.add_option('--mixing12', dest='mixing12',default=0.    ,type=float, help='Value of sinsqtheta12 or sinsq2theta12 depending on mixing param choice')
  parser.add_option('--mixing13', dest='mixing13',default=0.    ,type=float, help='Value of sinsqtheta13 or sinsq2theta13 depending on mixing param choice')
  parser.add_option('--mixing23', dest='mixing23',default=0.    ,type=float, help='Value of sinsqtheta23 or sinsq2theta23 depending on mixing param choice')
  parser.add_option('--dm12sq'  , dest='dm12sq'  ,default=7.5E-5,type=float, help='Value of small mass splitting')
  parser.add_option('--dm23sq'  , dest='dm23sq'  ,default=2.4E-3,type=float, help='Value of large mass splitting')
  parser.add_option('--deltacp' , dest='deltacp' ,default=0.    ,type=float, help='Value of delta CP')
  
  parser.add_option(     '--hierarchy'          , dest="hierarchy"     , default="normal"    , help='Mass hierarchy: "normal" or "inverted"')
  parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta", help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
  
  parser.add_option('-f', '--file', dest='output_file', default="validation.tex", help="Name for output TeX file.")
  
  parser.add_option('-c', '--config', dest='config', default="analysis_configuration.xml", help="Configuration file")
  parser.add_option('-l', '--label' , dest='label' , default=None, help="Label to put in file.")
  
  parser.add_option('-d', '--disable-banff', dest='nobanff', action="store_true",default=False)
  
  o,other = parser.parse_args()
  
  main(o.config,o.mixing12,o.mixing13,o.mixing23,o.dm12sq,o.dm23sq,o.deltacp,o.hierarchy,o.mixparamchoice,o.output_file,o.nobanff, o.label)
  
  #~ main("analysis_configuration_run1to5.xml",0.,0.,0.,7.5E-5,2.4E-3,0,"normal","sinsqtheta","validation.tex") 
  #~ main("analysis_configuration_run4.xml",0.8495,0.095,1.,7.54E-5,2.43E-3,0,"normal","sinsq2theta","validation.tex") # run1-4, PDG 2013?
  #~ main("analysis_configuration_run1to5.xml", 0.846,0.093,0.999,7.53E-5,2.44E-3,0.,"normal","sinsq2theta","validation.tex")
  
  if other:
    print "Warning - ignored arguments: ",other
