#!/usr/bin/env python

# A bunch of utilities used for testing the T2K run 10 spline files

import ROOT

#________________________________________________________________________________________
def print_spline_params_and_modes(tfile):
    import re
    import collections
    params_and_modes = collections.OrderedDict()
    keys = tfile.GetListOfKeys()
    for k in keys:
        name = k.GetName()
        if not "_sp_" in name:
            continue
        name_tmp0 = re.split("_sp_", name)[0]  # discards info on the bin numbers and part that labels it as a spline
        name_tmp1 = re.split("dev_", name_tmp0)[1]  # discards the tag at the start of the name
        param_mode = re.split("_(?!.*_)", name_tmp1)  # splits the param from the mode using the last occurrence of '_', so will not work if the mode name contains an underscore
        param = param_mode[0]
        mode = param_mode[1]
        if param not in params_and_modes:
            params_and_modes[param] = [mode]
        if mode not in params_and_modes[param]:
            params_and_modes[param].append(mode)
    for param in params_and_modes:
        print param
        for mode in params_and_modes[param]:
            print "   ", mode
        print ""

#________________________________________________________________________________________
def print_spline_nbins(tfile):
    import re
    import collections
    params_and_modes = collections.OrderedDict()
    keys = tfile.GetListOfKeys()
    for k in keys:
        name = k.GetName()
        if not "_sp_" in name:
            continue
        name_tmp0 = re.split("_sp_", name)[0]  # discards info on the bin numbers and part that labels it as a spline
        binning = re.split("_sp_", name)[1]
        bins = re.split("_", binning)
        name_tmp1 = re.split("dev_", name_tmp0)[1]  # discards the tag at the start of the name
        param_mode = re.split("_(?!.*_)", name_tmp1)  # splits the param from the mode using the last occurrence of '_', so will not work if the mode name contains an underscore
        param = param_mode[0]
        mode = param_mode[1]
        if param not in params_and_modes:
            params_and_modes[param] = dict()
        if mode not in params_and_modes[param]:
            params_and_modes[param][mode] = dict()
            params_and_modes[param][mode]["name"]  = param_mode
            params_and_modes[param][mode]["Etrue"] = 0
            params_and_modes[param][mode]["Ereco"] = 0
            params_and_modes[param][mode]["theta"] = 0
        params_and_modes[param][mode]["Etrue"] = max(int(bins[0])+1, params_and_modes[param][mode]["Etrue"])
        params_and_modes[param][mode]["Ereco"] = max(int(bins[1])+1, params_and_modes[param][mode]["Ereco"])
        if len(bins) == 3:
            params_and_modes[param][mode]["theta"] = max(int(bins[2])+1, params_and_modes[param][mode]["theta"])
        else:
            params_and_modes[param][mode]["theta"] = max(1, params_and_modes[param][mode]["theta"])
    print "printing: param_mode : Etrue bins ; Ereco bins ; theta bins"
    for param in params_and_modes:
        for mode in params_and_modes[param]:
            print params_and_modes[param][mode]["name"], ":", params_and_modes[param][mode]["Etrue"], ";", params_and_modes[param][mode]["Ereco"], ";", params_and_modes[param][mode]["theta"]
        print ""

#________________________________________________________________________________________
def get_spline_nbins(tfile):
    import re
    import collections
    params_and_modes = collections.OrderedDict()
    keys = tfile.GetListOfKeys()
    for k in keys:
        name = k.GetName()
        if not "_sp_" in name:
            continue
        name_tmp0 = re.split("_sp_", name)[0]  # discards info on the bin numbers and part that labels it as a spline
        binning = re.split("_sp_", name)[1]
        bins = re.split("_", binning)
        name_tmp1 = re.split("dev_", name_tmp0)[1]  # discards the tag at the start of the name
        param_mode = re.split("_(?!.*_)", name_tmp1)  # splits the param from the mode using the last occurrence of '_', so will not work if the mode name contains an underscore
        param = param_mode[0]
        mode = param_mode[1]
        if param not in params_and_modes:
            params_and_modes[param] = dict()
        if mode not in params_and_modes[param]:
            params_and_modes[param][mode] = dict()
            params_and_modes[param][mode]["name"]  = param_mode
            params_and_modes[param][mode]["Etrue"] = 0
            params_and_modes[param][mode]["Ereco"] = 0
            params_and_modes[param][mode]["theta"] = 0
        params_and_modes[param][mode]["Etrue"] = max(int(bins[0])+1, params_and_modes[param][mode]["Etrue"])
        params_and_modes[param][mode]["Ereco"] = max(int(bins[1])+1, params_and_modes[param][mode]["Ereco"])
        if len(bins) == 3:
            params_and_modes[param][mode]["theta"] = max(int(bins[2])+1, params_and_modes[param][mode]["theta"])
        else:
            params_and_modes[param][mode]["theta"] = max(1, params_and_modes[param][mode]["theta"])

    nbins_etrue = 0
    nbins_ereco = 0
    nbins_theta = 0
    for param in params_and_modes:
        for mode in params_and_modes[param]:
            nbins_etrue = max(nbins_etrue, params_and_modes[param][mode]["Etrue"])
            nbins_ereco = max(nbins_ereco, params_and_modes[param][mode]["Ereco"])
            nbins_theta = max(nbins_theta, params_and_modes[param][mode]["theta"])
    return [nbins_etrue, nbins_ereco, nbins_theta]

#________________________________________________________________________________________
def check_all_spline_binning():
    import sys
    runs = ["r9", "r10"]
    path = {}
    path["r9"] = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2018/xsec/splines2018v1/"
    path["r10"] = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/xsec/splines2020v3/"
    tag1 = {}
    tag1["r9"] = "spline2018v1"
    tag1["r10"] = "spline2019splinesfitqunPre"
    tag2 = {}
    tag2["r9"] = "_t2ksk.14a.neut5.3.2.13a_tuned_v1r1."
    tag2["r10"] = "_t2ksk19b.fqv4r0."
    currents = ["fhc","rhc"]
    current_name = {}
    current_name["fhc"] = {}
    current_name["rhc"] = {}
    current_name["fhc"]["r9"] = "250ka"
    current_name["rhc"]["r9"] = "m250ka"
    current_name["fhc"]["r10"] = "fhc"
    current_name["rhc"]["r10"] = "rhc"
    fine = {}
    fine["r9"] = "fine."
    fine["r10"] = ""
    flavs = ["numu_x_numu", "numu_x_nue", "numubar_x_numubar", "numubar_x_nuebar"]
    selecs = ["numu", "nue", "cc1pi"]
    for current in currents:
        for selec in selecs:
            for flav in flavs:
                filenames = {}
                files = {}
                bins = {}
                for run in runs:
                    tag = ""
                    if selec == "numu":
                        tag = tag1[run] + tag2[run]
                    else:
                        tag = tag1[run] + "erectheta" + tag2[run]
                    filenames[run] = path[run] + tag + current_name[current][run] + "." + fine[run] + flav + "_" + selec + "selec.root"
                    files[run] = ROOT.TFile(filenames[run], "read")
                    if not files[run]:
                        print "ERROR: could not open file:", filenames[run]
                        sys.exit()
                    bins[run] = get_spline_nbins(files[run])
                print current, flav, selec, "\n    run9 bins:", bins["r9"], "\n    run10 bins:", bins["r10"]

#________________________________________________________________________________________
def is_spline_flat(spline):
    nk = spline.GetNp()
    for ik in range(0, nk):
        kx = ROOT.Double(-999)
        ky = ROOT.Double(-999)
        spline.GetKnot(ik, kx, ky)
        if ROOT.TMath.Abs(ky - 1.0) > 1E-10:
            return False
    return True

#________________________________________________________________________________________
def print_nonflat_spline_names(tfile, param = ""):
    keys = tfile.GetListOfKeys()
    for k in keys:
        name = k.GetName()
        if not "_sp_" in name:
            continue
        if param and not param in name:
            continue
        if not is_spline_flat(tfile.Get(name)):
            print name

#________________________________________________________________________________________
def print_nonflat_valor_splines_for_param(tfile, param_name):
    for k in tfile.GetListOfKeys():
        name = k.GetName()
        if (not "spl_" in name) or (not param_name in name):
            continue
        arr = k.ReadObj()
        for spl in arr:
            if not is_spline_flat(spl):
                print name

#________________________________________________________________________________________
def print_valor_nonflat_spline_params_and_reac_types(tfile):
    import re
    import collections
    import sys
    params_and_modes = collections.OrderedDict()
    keys = tfile.GetListOfKeys()
    for k in keys:
        name = k.GetName()
        if not "spl_" in name:
            continue
        arr = k.ReadObj()
        any_nonflat_spline = False
        for spl in arr:
            if not is_spline_flat(spl):
                any_nonflat_spline = True
                break
        if not any_nonflat_spline:
            continue
        mode_param = re.split("spl_", name)[1]  # discards "spl_"
        mode  = re.split("_",mode_param)[0]
        param = mode_param.replace(mode+"_", "")
        reac_type = ""
        if "CC" in mode:
            reac_type = "CC" + re.split("CC",mode)[1]
        elif "NC" in mode:
            reac_type = "NC" + re.split("NC",mode)[1]
        else:
            "ERROR: could not parse mode named \'{}\'".format(mode)
            sys.exit()

        if param not in params_and_modes:
            params_and_modes[param] = [reac_type]
        if reac_type not in params_and_modes[param]:
            params_and_modes[param].append(reac_type)
    for param in params_and_modes:
        print param
        for reac_type in params_and_modes[param]:
            print "   ", reac_type
        print ""
