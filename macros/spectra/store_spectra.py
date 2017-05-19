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
processes      = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother", "NC1gamma"]
processes_used = ["CCQE", "CC1pi", "CCcoherent", "2p2h", "CCother", "NC1pi0", "NC1piPM", "NCcoherent", "NCother", "NC1gamma"]

samples   = []
if t2kopt.fFit_1Rmu:
    samples.append("1Rmu")
if t2kopt.fFit_1Re:
    samples.append("1Re")
if t2kopt.fFit_RHC_1Rmu:
    samples.append("RHC1Rmu")
if t2kopt.fFit_RHC_1Re:
    samples.append("RHC1Re")
if t2kopt.fFit_MultiRmu:
    samples.append("MultiRmu")
if t2kopt.fFit_MultiRe:
    samples.append("MultiRe")

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
modes_to_species_and_process[conv.kSKmOscNuECCQE]       = ("signue","CCQE")
modes_to_species_and_process[conv.kSKmOscNuECC1pi]      = ("signue","CC1pi")
modes_to_species_and_process[conv.kSKmOscNuECCcoh]      = ("signue","CCcoherent")
modes_to_species_and_process[conv.kSKmOscNuECC2p2h]     = ("signue","2p2h")
modes_to_species_and_process[conv.kSKmOscNuECCoth]      = ("signue","CCother")

def main(config_file, mixing12, mixing13, mixing23, dm12sq, dm23sq, deltacp,
    hierarchy, angle_scheme, octant, dmsqscheme, output_file):
    opt.gMassHierarchyAssumed = util.HierarchyFromString(hierarchy)
    opt.gOscMixingParamChoice = util.MixingTypeFromString(angle_scheme)
    opt.gOscOctantChoice = util.OctantFromString(octant)
    opt.gOscDmsqScheme = util.DmsqSchemeFromString(dmsqscheme)
  
    opt.SetDefaultParameterValues()
  
    phyp = valor.PhysHypothesis("3f_osc")

    skexp = valor.SuperkExpectation()
    skexp.Read(os.environ[env_var]+"/config/"+config_file)
    skmc = skexp.GetSuperkMC()
  
    phyp.SetCurrValue(conv.kPhypMixing12,  mixing12);
    phyp.SetCurrValue(conv.kPhypMixing13,  mixing13);
    phyp.SetCurrValue(conv.kPhypMixing23,  mixing23);
    phyp.SetCurrValue(conv.kPhypCPV     ,   deltacp);
    phyp.SetCurrValue(conv.kPhypDm12sq  ,    dm12sq);
    phyp.SetCurrValue(conv.kPhypDm23sq  ,    dm23sq);
  
    skexp.Generate(phyp)
 
    file = ROOT.TFile(output_file, "RECREATE")
    file.cd()
    for isample, sample_name in enumerate(samples):
        s = util.SKSampleFromString(sample_name)
        for mode, (spec, proc) in modes_to_species_and_process.iteritems():
            spectrum = skexp.CurrPredictedErecoVsQ2Spectrum(s,mode)
            if not spectrum:
                continue
            # Write spectrum to file
            spectrum.Write("spectrum_" + util.AsString(s) + "_" + util.AsString(mode))
    file.Close()
    
if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Make tables of validation spectra")
    parser.add_option('--mixing12', dest='mixing12', default=0., type=float, help='Value of sinsqtheta12')
    parser.add_option('--mixing13', dest='mixing13', default=0., type=float, help='Value of sinsqtheta13')
    parser.add_option('--mixing23', dest='mixing23', default=0., type=float, help='Value of sinsqtheta23')
    parser.add_option('--dm12sq', dest='dm12sq', default=7.53E-5, type=float, help='Value of small mass splitting')
    parser.add_option('--dm23sq', dest='dm23sq', default=2.509E-3, type=float, help='Value of large mass splitting')
    parser.add_option('--deltacp', dest='deltacp', default=0.0, type=float, help='Value of delta CP')
    parser.add_option('--hierarchy', dest="hierarchy", default="normal", help='Mass hierarchy: "normal" or "inverted"')
    parser.add_option('-m','--mixing-param-choice', dest="mixparamchoice", default="sinsqtheta", help='Mixing param coice: "sinsqtheta" or "sinsq2theta"')
    parser.add_option('-t','--mixing-param-octant', dest="mixparamoctant", default="two", help='Mixing param octant: "one" or "two"')
    parser.add_option('-a','--dmsq-scheme', dest="dmsqscheme", default="auto", help='"abs23", "abs31", "auto" or "fogli"')
    parser.add_option('-f', '--file', dest='output_file', default="spectra.root", help="Name for output root file.")
    parser.add_option('-c', '--config', dest='config', default="analysis_configuration_oa2017_run1-8_postfit.xml", help="Configuration file")

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults
  
    o, other = parser.parse_args()
  
    main(o.config, o.mixing12, o.mixing13, o.mixing23, o.dm12sq, o.dm23sq,
        o.deltacp, o.hierarchy, o.mixparamchoice, o.mixparamoctant, o.dmsqscheme,
        o.output_file)
    
    if other:
        print "Warning - ignored arguments: ",other
