#!/usr/bin/env python

import ROOT
import math
import array
import numpy as np
import sys
import os
import ctypes

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

env_var = "VALOR_EXPT"
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short =2**32 -1

# File naming convention
pattern = "[VERSION]_t2ksk19b.fqv4r0.[CURRENT].[TYPE]_[SELECTION].root"
models = dict()

def get_file_suffix(sample):
    if "1Rmu" in sample:
        return ""
    else:
        return "erectheta"
    

# This is a list containing parameters and the modes they apply to, in the format "param_mode",
# using raw (T2KReWeight/NDfit) naming
def get_systematics():
    # TODO: Add effective Eb splines
    return [
        "maqe_ccqe",
        "mecpdd_mec",
        "ca5_cc1pi", "ca5_ncpiz", "ca5_ncpipm",
        "manff_cc1pi", "manff_ncpiz", "manff_ncpipm",
        "bgscl_cc1pi", "bgscl_ncpiz", "bgscl_ncpipm",
        "bgsclbar_cc1pi", "bgsclbar_ncpiz", "bgsclbar_ncpipm",
        "FSIPIABS_ccqe", "FSIPIABS_cc1pi", "FSIPIABS_cccoh", "FSIPIABS_ccmisc",
            "FSIPIABS_ccmpi", "FSIPIABS_ccdis", "FSIPIABS_ncpiz", "FSIPIABS_ncpipm",
            "FSIPIABS_nccoh", "FSIPIABS_ncoth", "FSIPIABS_mec", "FSIPIABS_nc1gamma",
        "FSICEXLO_ccqe", "FSICEXLO_cc1pi", "FSICEXLO_cccoh", "FSICEXLO_ccmisc",
            "FSICEXLO_ccmpi", "FSICEXLO_ccdis", "FSICEXLO_ncpiz", "FSICEXLO_ncpipm",
            "FSICEXLO_nccoh", "FSICEXLO_ncoth", "FSICEXLO_mec", "FSICEXLO_nc1gamma",
        "FSIINELLO_ccqe", "FSIINELLO_cc1pi", "FSIINELLO_cccoh", "FSIINELLO_ccmisc",
            "FSIINELLO_ccmpi", "FSIINELLO_ccdis", "FSIINELLO_ncpiz", "FSIINELLO_ncpipm",
            "FSIINELLO_nccoh", "FSIINELLO_ncoth", "FSIINELLO_mec", "FSIINELLO_nc1gamma",
        "FSIPIPROD_ccqe", "FSIPIPROD_cc1pi", "FSIPIPROD_cccoh", "FSIPIPROD_ccmisc",
            "FSIPIPROD_ccmpi", "FSIPIPROD_ccdis", "FSIPIPROD_ncpiz", "FSIPIPROD_ncpipm",
            "FSIPIPROD_nccoh", "FSIPIPROD_ncoth", "FSIPIPROD_mec", "FSIPIPROD_nc1gamma",
        "FSIINELHI_ccqe", "FSIINELHI_cc1pi", "FSIINELHI_cccoh", "FSIINELHI_ccmisc",
            "FSIINELHI_ccmpi", "FSIINELHI_ccdis", "FSIINELHI_ncpiz", "FSIINELHI_ncpipm",
            "FSIINELHI_nccoh", "FSIINELHI_ncoth", "FSIINELHI_mec", "FSIINELHI_nc1gamma",
        "DISBY_ccdis",
        "MPiBY_ccmpi",
        "MPiAGKYXsec_ccmpi",
        "2p2hedeplowenu_mec",
        "2p2hedephienu_mec",
        "2p2hedeplowenubar_mec",
        "2p2hedephienubar_mec",
        "LowQ2Suppression0_ccqe",
        "LowQ2Suppression1_ccqe",
        "LowQ2Suppression2_ccqe",
        "LowQ2Suppression3_ccqe",
        "LowQ2Suppression4_ccqe"
    ]

def get_mode_map():
    valor_modes = ["CCQE","CCcoh","CCMEC","CCDIS","CCMpi","CCmisc","CC1piC","NCcoh","NCoth","NC1piC","NC1pi0","NC1Gamma"]
    raw_modes   = ["ccqe","cccoh","mec"  ,"ccdis","ccmpi","ccmisc","cc1pi" ,"nccoh","ncoth","ncpipm","ncpiz" ,"nc1gamma"]
    return list(zip(raw_modes, valor_modes))

def get_flavour_map():
    raw_flavours = ["numu_x_numu","numubar_x_numubar",
        "nue_x_nue","nuebar_x_nuebar",
        "numu_x_nue","numubar_x_nuebar"]
    valor_flavours = ["numu","numubar",
        "nue","nuebar",
        "numutonue","numubartonuebar"]
    return list(zip(raw_flavours, valor_flavours))

def get_param_map():
    # NOTE (2020-01-03, run 10):
    #   The mapping between the raw (i.e. T2KReWeight/NDfit) and VALOR/NIWG param names for the FSI params (FSI.. or FEF...)
    #   is not intuitive, but they are listed in table 9 of TN-325 (see https://t2k.org/docs/technotes/325/tn-325-v1).
    #   In fact, the raw names use some old conventions, and they will probably change to being the same as the
    #   VALOR/NIWG names next time splines are generated.
    # NOTE:
    #   Make sure that if any param name is a subset of another param name, then in the lists below, it appears after
    #   the param it is a subset of. Otherwise, the wrong param may be fetched.
    raw_params   = ["ca5", "bgsclbar"      , "bgscl"  , "maqe", "manff", "mecpdd"      ,"FSIPIABS", "FSICEXLO", "FSIINELLO", "FSIPIPROD", "FSIINELHI", "DISBY"    , "MPiBY"    , "MPiAGKYXsec" , "2p2hedeplowenu"  , "2p2hedephienu"    , "2p2hedeplowenubar"  , "2p2hedephienubar"    , "LowQ2Suppression0", "LowQ2Suppression1", "LowQ2Suppression2", "LowQ2Suppression3", "LowQ2Suppression4"]
    valor_params = ["CA5", "ISO_BKG_LowPPi", "ISO_BKG", "MAQE", "MARES", "2p2h_shape_O","FEFABS"  , "FEFCX"   , "FEFQE"    , "FEFINEL"  , "FEFQEH"   , "CC_BY_DIS", "CC_BY_MPi", "CC_AGKY_Mult", "2p2h_Edep_lowEnu", "2p2h_Edep_highEnu", "2p2h_Edep_lowEnubar", "2p2h_Edep_highEnubar", "Q2_norm_0"        , "Q2_norm_1"        , "Q2_norm_2"        , "Q2_norm_3"        , "Q2_norm_4"        ]
    return list(zip(raw_params, valor_params))

def get_beam(sample):
    if sample in ["FHC1Rmu","FHC1Re","FHCMultiRe"]:
        return "fhc"
    elif sample in ["RHC1Rmu","RHC1Re"]:
        return "rhc"
    else:
        print "Error: Unknown sample ", sample
        sys.exit()

def get_selection_map(sample):
    if "1Re" in sample:
        return "nueselec"
    elif "1Rmu" in sample:
        return "numuselec"
    elif "MultiRe" in sample:
        return "cc1piselec"
    else:
        print "Error: Unknown sample ", sample
        sys.exit()

def get_sdk_selection(sample):
    if "1Rmu" in sample:
        return "1rmu"
    elif "1Re" in sample:
        return "1re"
    else:
        return "multire"

def get_sdk_beam(sample):
    if "RHC" in sample:
        return "JPARC_RHC"
    else:
        return "JPARC_FHC"

def build_spline_name(sample, syst, flavour, mode, true_bin, reco_bin):
    beam = get_sdk_beam(sample)
    selection = get_sdk_selection(sample)
    new_syst = syst
    for (old, new) in get_param_map():
        new_syst = new_syst.replace(old, new)
    return "spl_{}_{}_{}{}_{}_{}_{}".format(beam, selection, flavour, mode, new_syst, true_bin, reco_bin)

def build_spline_array_name(syst, flavour, mode):
    new_syst = syst
    for (old, new) in get_param_map():
        new_syst = new_syst.replace(old, new)
    return "spl_{}{}_{}".format(flavour, mode, new_syst)

def get_fit_elements(config_file):
    elements = valor.CompoundFitElements()
    elements.Read(util.CompleteConfigPath(config_file), "nomc+nosplines")
    return elements

def set_model_map(model):
    if model.Name() == "superk__JPARC_FHC__1rmu":
        models["FHC1Rmu"] = model
    elif model.Name() == "superk__JPARC_RHC__1rmu":
        models["RHC1Rmu"] = model
    elif model.Name() == "superk__JPARC_FHC__1re":
        models["FHC1Re"] = model
    elif model.Name() == "superk__JPARC_RHC__1re":
        models["RHC1Re"] = model
    elif model.Name() == "superk__JPARC_FHC__multire":
        models["FHCMultiRe"] = model
    else:
        print "Unrecognised model ", model.Name()
        sys.exit()

# The raw spline file specified bins as true_reco1_reco2 (or true_reco1 if 2D), but in VALOR we put all
# true kinevars into 1 dimension and all reco into another, so we want to convert the reco bin number
# into VALOR format
def get_new_binning(sample, flavour, mode, true_bin, reco_bin1, reco_bin2=0):
    # Flavour and mode included just in case this is required for true binning
    # in the future. For now, the simple 1D true binning means the index we
    # receive is the same as the format required
    binning = models[sample].GetBinning()
    reco_binning = binning.XsecRecoBinning()
    reco_edges1 = reco_binning.BinEdgeArr(0)
    reco_edges2 = reco_binning.BinEdgeArr(1)
    reco_centre1 = 0.5 * (reco_edges1[reco_bin1] + reco_edges1[reco_bin1 + 1])
    # Even the muon samples are 2D internally, so reco_bin2 default value works
    reco_centre2 = 0.5 * (reco_edges2[reco_bin2] + reco_edges2[reco_bin2 + 1])
    reco_vec = ROOT.vector('double')()
    reco_vec.push_back(reco_centre1)
    reco_vec.push_back(reco_centre2)
    new_true_bin = true_bin
    new_reco_bin = reco_binning.FindBin(reco_vec)
    return new_true_bin, new_reco_bin

def get_new_binning_nbins(sample, flavour, mode):
    # Flavour included just in case this is required for true binning
    # in the future.
    reco_nbins = 1
    true_nbins = 1
    m = util.ReacModeFromString(flavour+mode)
    if m == conv.kmUndefined or m == python_undefined or m == python_undefined_short:
        print "ERROR: could not convert input mode [", mode, "] into a VALOR mode"
        util.Die(-1)
    binning = models[sample].GetBinning()
    reco_binning = binning.XsecRecoBinning()
    true_binning = binning.XsecTrueBinning(m)
    for i in range(0, reco_binning.NDimensions()):
        reco_nbins *= reco_binning.NBins(i)
    for i in range(0, true_binning.NDimensions()):
        true_nbins *= true_binning.NBins(i)
    return true_nbins, reco_nbins

def main(version, path, sample, config):
    modes = set()
    params = set()
    elements = get_fit_elements(config)
    compound_model = elements.FitModel()
    for i in range(compound_model.GetNSubModels()):
        selection_model = compound_model.GetSubModel(i)
        set_model_map(selection_model)

    current   = get_beam(sample)
    selection = get_selection_map(sample)

    output_file  = ROOT.TFile("splines_"+version+"_"+get_sdk_beam(sample)+"_"+get_sdk_selection(sample)+".root","RECREATE") # Output for this sample's splines

    for orig_type,valor_type in get_flavour_map():
        # Find the right input file
        infile = pattern
        indir  = path
        for old,new in [("[CURRENT]"  ,current),
                        ("[TYPE]"     ,orig_type),
                        ("[SELECTION]",selection),
                        ("[VERSION]"  ,version) ]:
            infile   = infile .replace(old,new)
            indir    = indir.replace(old,new)
  
        print "Processing file: ",infile
        path_in = indir+"/"+infile
        tf = ROOT.TFile(path_in)
        if not tf:
            print "Failed to open file:\n",path_in
  
        # Loop over spline types - systematic+mode
        for syst in get_systematics():
            par  = ""
            mode = ""
            if "erpa" in syst:
                ppar,coeff,mode = syst.split("_")
                par = ppar +"_"+ coeff
            else:  
                par,mode = syst.split("_")
            for (old,new) in get_mode_map():
                mode = mode.replace(old,new)

            mode_string = valor_type+mode
            m_enum = util.ReacModeFromString(mode_string)
            if m_enum == conv.kmUndefined or m_enum == python_undefined or m_enum == python_undefined_short:
                print "Skipping param+mode combination with undefined mode:", par, mode_string
                if "nc" in util.ToLowerCase(mode_string) and "to" in util.ToLowerCase(valor_type):
                    print "(note: skipping oscillated NC modes is expected)"
                continue
        
            # Get all the splines for this systematic+mode
            spline_dict = dict()
            for key in tf.GetListOfKeys():
                name = key.GetName()
                if not (syst in name and "_sp_" in name):
                    continue
                tokens = name.split("_")
                p = tokens[1]
                if p == "erpa":
                    p = "erpa_" + tokens[2]
                    m = tokens[3]
                else:
                    m = tokens[2]
                modes.add(m)
                params.add(p)
                spline_dict[name] = key.ReadObj()
      
            output_file.cd()

            # Setup the array that holds all splines for this systematic + mode
            nbins_true, nbins_reco = get_new_binning_nbins(sample, valor_type, mode)
            spline_array = ROOT.TObjArray(nbins_true * nbins_reco)

            # Loop over the splines - rename them and them to the spline array for this param/mode combination.
            for orig_spl_name, spline in spline_dict.iteritems():
                junk,bins = orig_spl_name.rsplit("_sp_")
       
                all_bins = map(int,bins.split("_"))
                if len(all_bins) == 2:
                    true_bin, reco_bin = map(int, bins.split("_"))
                    true_bin, reco_bin = get_new_binning(sample, valor_type, mode, true_bin, reco_bin)
                elif len(all_bins) == 3:
                    true_bin, reco_bin1, reco_bin2 = map(int, bins.split("_"))
                    true_bin, reco_bin = get_new_binning(sample, valor_type, mode, true_bin, reco_bin1, reco_bin2)
                else:
                    print "Unexpected binning ", bins
                    sys.exit()
              
                if reco_bin < 0:
                    continue
                new_spl_name = build_spline_name(sample, par, valor_type, mode, true_bin, reco_bin)
              
                nknots = spline.GetNp()
                knots = np.zeros(nknots,dtype=float)
              
                new_spline = ROOT.TSpline3(spline)
                new_spline.SetName(new_spl_name)

                idx = true_bin * nbins_reco + reco_bin
                if idx >= nbins_true * nbins_reco:
                    print "ERROR: Tried to add spline \'{}\'".format(new_spl_name), "to array", \
                        build_spline_array_name(par, valor_type, mode), ", but the index ({})".format(idx), \
                        "is out of range [0, {}]".format(nbins_true * nbins_reco), "; true bin =", true_bin, \
                        "(max = {})".format(nbins_true), "; reco bin =", reco_bin, "(max = {})".format(nbins_reco), \
                        "; input spline bins =", all_bins
                    sys.exit()
                #if build_spline_array_name(par, valor_type, mode) == "spl_numuCCQE_MAQE":
                #    print "Adding spline", new_spl_name, "to array", build_spline_array_name(par, valor_type, mode), "at idx =", idx
                spline_array.AddAt(new_spline, idx)

            print "Saving spline array named", build_spline_array_name(par, valor_type, mode), "\n"
            spline_array.Write(build_spline_array_name(par, valor_type, mode), 1)   # The 1 makes the array be written as a single object instead of saving each spline individually

            # The legacy T2K analysis did not distinguish between CC1piC and CC1pi0, but we do in the SDK. So
            # copy the newly made CC1piC spline into CC1pi0 as well (use the same spline for both modes)
            if mode == "CC1piC":
                spline_array.Write(build_spline_array_name(par, valor_type, "CC1pi0"), 1)

        tf.Close()
    output_file.Close()
    print "Params:"
    for p in params:
        print p
    print "Modes:"
    for m in modes:
        print m

    print "Done!"

  
if __name__=="__main__":
    import argparse
    samples = ["FHC1Rmu", "RHC1Rmu", "FHC1Re", "RHC1Re", "FHCMultiRe"]
    parser = argparse.ArgumentParser(usage="%prog [options]..",description="Reformat splines")
    parser.add_argument('-p', '--prefix', default='spline2018v1', help="Spline filename prefix.")
    parser.add_argument('-i', '--input', help="Input spline directory.")
    parser.add_argument('-s', '--sample', choices=samples, help="The SK sample.")
    parser.add_argument('-c', '--config', help="Analysis configuration file.")
    args = parser.parse_args()
    main(args.prefix + get_file_suffix(args.sample), args.input, args.sample, args.config)
 
