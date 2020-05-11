#!/usr/bin/env python

# Example usage:
# python $VALOR_EXPT/src/macros/preprocessing/dummy_spline_format.py --prefix spline2018v1 --input $VALOR_EXPT/data/dummy_data/xsec --config $VALOR_EXPT/config/analysis/dummy_wc_1sample.xml --sample 1Rmu

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
pattern = "[VERSION]_t2ksk.14a.neut5.3.2.13a_tuned_v1r1.[CURRENT].fine.[TYPE]_[SELECTION].root"
models = dict()

def get_file_suffix(sample):
    if "1Rmu" in sample:
        return ""
    else:
        return "erectheta"
    

def get_systematics():
    return ["maqe_ccqe"]

def get_mode_map():
    valor_modes = ["CCQE"]
    raw_modes = ["ccqe"]
    return list(zip(raw_modes, valor_modes))

def get_flavour_map():
    raw_flavours   = ["numu_x_numu", "nue_x_nue", "numu_x_nue"]
    valor_flavours = ["numu"       , "nue"      , "numutonue" ]
    return list(zip(raw_flavours, valor_flavours))

def get_param_map():
    raw_params = ["maqe"]
    valor_params = ["MAQE"]
    return list(zip(raw_params, valor_params))

def get_beam(sample):
    if sample in ["1Rmu","1Re","MultiRe"]:
        return "250ka"
    elif sample in ["RHC1Rmu","RHC1Re"]:
        return "m250ka"
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
    elements.Read(config_file, "nomc+nosplines")
    return elements

def set_model_map(model):
    if model.Name() == "superk__JPARC_FHC__1rmu":
        models["1Rmu"] = model
    elif model.Name() == "superk__JPARC_RHC__1rmu":
        models["RHC1Rmu"] = model
    elif model.Name() == "superk__JPARC_FHC__1re":
        models["1Re"] = model
    elif model.Name() == "superk__JPARC_RHC__1re":
        models["RHC1Re"] = model
    elif model.Name() == "superk__JPARC_FHC__multire":
        models["MultiRe"] = model
    else:
        print "Unrecognised model ", model.Name()
        sys.exit()

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
    samples = ["1Rmu", "RHC1Rmu", "1Re", "RHC1Re", "MultiRe"]
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--prefix', default='spline2018v1', help="Spline filename prefix.")
    parser.add_argument('-i', '--input', help="Input spline directory.")
    parser.add_argument('-s', '--sample', choices=samples, help="The SK sample.")
    parser.add_argument('-c', '--config', help="Analysis configuration file.")
    args = parser.parse_args()
    main(args.prefix + get_file_suffix(args.sample), args.input, args.sample, args.config)
 
