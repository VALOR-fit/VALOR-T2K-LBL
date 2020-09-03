#!/usr/bin/env python

import ROOT
import sys
from copy import deepcopy

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

valor_spline_file_tag = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/xsec/splines2020v6/valor_format/splines_spline2019splinesfitqunPre"
banff_filename  = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/banff/postfit_data_20200409/BANFF_postfit_OA2020_fixed2p2h_fixedIHalfLowPPi_20200409_valor_order_all.root"
ptheta_spline_file_tag = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/t2k_run1-10/new_sk_matrix_banff_postfit_data_20200409/syst_var/other_groups/p-theta_splines/RF/"

# Get p-theta mode enum with e.g. valor_modes_in_ptheta_order_all.index("CCcoherent")
valor_modes_in_ptheta_order_all = [
    "CCQE",
    "CCMEC",
    "CC1piC",
    "CCcoherent",
    "CCDIS",
    "CCMpi",
    "CCmisc",
    "NC1pi0",
    "NC1piC",
    "NCcoherent",
    "NC1Gamma",
    "NCother"
]

# Get p-theta param enum with e.g. 1 + valor_params_in_ptheta_order.index("2p2h_shape_O")
valor_params_in_ptheta_order = [
    "MAQE",
    "2p2h_shape_O",
    "2p2h_Edep_lowEnu",
    "2p2h_Edep_highEnu",
    "2p2h_Edep_lowEnubar",
    "2p2h_Edep_highEnubar",
    "CA5",
    "MARES",
    "ISO_BKG_LowPPi",
    "ISO_BKG",
    "CC_BY_DIS",
    "CC_BY_MPi",
    "CC_AGKY_Mult",
    "FEFQE",
    "FEFABS",
    "FEFCX",
    "FEFQEH",
    "FEFINEL",
    "Q2_norm_0",
    "Q2_norm_1",
    "Q2_norm_2",
    "Q2_norm_3",
    "Q2_norm_4",
    "Q2_norm_5",
    "Q2_norm_6",
    "Q2_norm_7",
    "EB_dial_O_nu",
    "EB_dial_O_nubar"
]

# p-theta's spline knot x values are (actual parameter value / scale from below list)
ptheta_xscales = {
    "MAQE"                 : 0.06,
    "2p2h_shape_O"         : 3.00,
    "2p2h_Edep_lowEnu"     : 1.00,
    "2p2h_Edep_highEnu"    : 1.00,
    "2p2h_Edep_lowEnubar"  : 1.00,
    "2p2h_Edep_highEnubar" : 1.00,
    "CA5"                  : 0.15,
    "MARES"                : 0.15,
    "ISO_BKG_LowPPi"       : 1.3,
    "ISO_BKG"              : 0.15,
    "CC_BY_DIS"            : 3.,
    "CC_BY_MPi"            : 3.,
    "CC_AGKY_Mult"         : 3.,
    "FEFQE"                : 0.313,
    "FEFABS"               : 0.859,
    "FEFCX"                : 1.101,
    "FEFQEH"               : 0.432,
    "FEFINEL"              : 0.305,
    "Q2_norm_0"            : 1.,
    "Q2_norm_1"            : 1.,
    "Q2_norm_2"            : 1.,
    "Q2_norm_3"            : 1.,
    "Q2_norm_4"            : 1.,
    "Q2_norm_5"            : 1.,
    "Q2_norm_6"            : 1.,
    "Q2_norm_7"            : 1.,
    "EB_dial_O_nu"         : 1.,
    "EB_dial_O_nubar"      : 1.
}

valor_sample_full_names = {
    "FHC1Rmu"    : "superk__JPARC_FHC__1rmu",
    "RHC1Rmu"    : "superk__JPARC_RHC__1rmu",
    "FHC1Re"     : "superk__JPARC_FHC__1re",
    "RHC1Re"     : "superk__JPARC_RHC__1re",
    "FHCMultiRe" : "superk__JPARC_FHC__multire"
}

sample_is_2d = {
    "FHC1Rmu"    : False,
    "RHC1Rmu"    : False,
    "FHC1Re"     : True,
    "RHC1Re"     : True,
    "FHCMultiRe" : True
}

valor_sample_to_ptheta_sample = {
    "FHC1Rmu"    : "numode_numu",
    "RHC1Rmu"    : "antinu_numu",
    "FHC1Re"     : "numode_nue",
    "RHC1Re"     : "antinu_nue",
    "FHCMultiRe" : "numode_cc1pi"
}


# p-theta xsec Etrue edges are shared between all samples and modes
ptheta_etrue_bin_edges = [0.0, 0.4, 0.5, 0.6, 0.8, 1.2, 2.0, 3.0, 5.0, 30.0]

# Usually to get the true binning, one just should do sskb.XsecTrueBinning(m).BinEdgeArr(0),
# but for some reason this is returning millions of junk entries onto the end of the list,
# which really slows the program down, so I'm hardcoding the binning instead.
valor_etrue_bin_edges_rightsign_nue = [
    0.0,	0.5,	0.7,	0.8,	1.5,	2.5,
    4.0,	30.0
]
valor_etrue_bin_edges_rightsign_numu_and_osc_nue = [
    0.0,	0.4,	0.5,	0.6,	0.7,	1.0,
    1.5,	2.5,	3.5,	5.0,	7.0,
    30.0
]
valor_etrue_bin_edges_wrongsign_nue = [
    0.0,	2.5,	30.0
]
valor_etrue_bin_edges_wrongsign_numu_and_osc_nue = [
    0.0,	0.7,	1.0,	1.5,	2.5,	30.0
]

# Erec edges for e-like samples
ptheta_erec_bin_edges = [
    0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450,
    0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950,
    1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450,
    1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950,
    2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450,
    2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950,
    3.000,  3.250,  3.500,  3.750,  4.000,  4.500,  5.000,  5.500,  6.000,  7.000,
    8.000,  9.000, 10.000, 30.000
]

valor_erec_bin_edges_mulike = [
    0.0,	0.3,	0.35,	0.4,	0.45,	0.5,
    0.55,	0.6,	0.65,	0.7,	0.75,
    0.8,	0.85,	0.9,	0.95,	1.0,
    1.1,	1.2,	1.3,	1.4,	1.5,
    1.6,	1.7,	1.8,	1.9,	2.0,
    2.1,	2.2,	2.3,	2.4,	2.5,
    2.6,	2.7,	2.8,	2.9,	3.0,
    3.5,	4.0,	5.0,	7.0,	10.0,
    30.0
]
valor_erec_bin_edges_elike = [
    0.0,	0.05,	0.1,	0.15,	0.2,	0.25,
    0.3,	0.35,	0.4,	0.45,	0.5,
    0.55,	0.6,	0.65,	0.7,	0.75,
    0.8,	0.85,	0.9,	0.95,	1.0,
    1.05,	1.1,	1.15,	1.2,	1.25
]

ptheta_theta_bin_edges = [
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 180.0
]

valor_theta_bin_edges_elike = [
    0.0,	10.0,	20.0,	30.0,	40.0,	50.0,
    60.0,	70.0,	80.0,	90.0,	100.0,
    110.0,	120.0,	130.0,	140.0,	180.0
]

valor_etrue_fit_bin_edges = [
    0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.325,
    0.350,  0.375,  0.400,  0.425,  0.450,  0.475,  0.500,  0.525,
    0.550,  0.575,  0.600,  0.625,  0.650,  0.675,  0.700,  0.725,
    0.750,  0.775,  0.800,  0.825,  0.850,  0.875,  0.900,  0.925,
    0.950,  0.975,  1.000,  1.050,  1.100,  1.150,  1.200,  1.250,
    1.300,  1.350,  1.400,  1.450,  1.500,  1.550,  1.600,  1.650,
    1.700,  1.750,  1.800,  1.850,  1.900,  1.950,  2.000,  2.050,
    2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450,
    2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,
    2.900,  2.950,  3.000,  3.100,  3.200,  3.300,  3.400,  3.500,
    4.000,  5.000,  7.000, 10.000,  30.000
]
valor_erec_fit_bin_edges_mulike = [
    0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450,
    0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950,
    1.000,  1.050,  1.100,  1.150,  1.200,  1.250,  1.300,  1.350,  1.400,  1.450,
    1.500,  1.550,  1.600,  1.650,  1.700,  1.750,  1.800,  1.850,  1.900,  1.950,
    2.000,  2.050,  2.100,  2.150,  2.200,  2.250,  2.300,  2.350,  2.400,  2.450,
    2.500,  2.550,  2.600,  2.650,  2.700,  2.750,  2.800,  2.850,  2.900,  2.950,
    3.000,  3.250,  3.500,  3.750,  4.000,  4.500,  5.000,  5.500,  6.000,  7.000,
    8.000,  9.000, 10.000, 30.000
]

valor_erec_fit_bin_edges_elike = [
    0.000,  0.050,  0.100,  0.150,  0.200,  0.250,  0.300,  0.350,  0.400,  0.450,
    0.500,  0.550,  0.600,  0.650,  0.700,  0.750,  0.800,  0.850,  0.900,  0.950,
    1.000,  1.050,  1.100,  1.150,  1.200,  1.250
]



# Overrides the nominal +/- 3sigma values for the following params
special_case_tweak_p3 = {
    "2p2h_shape_O"         : 1.0,
    "2p2h_Edep_lowEnu"     : 0.75,
    "2p2h_Edep_highEnu"    : 0.75,
    "2p2h_Edep_lowEnubar"  : 0.75,
    "2p2h_Edep_highEnubar" : 0.75,
    "EB_dial_O_nu"         : 15.0,
    "EB_dial_O_nubar"      : 15.0
}
special_case_tweak_m3 = {
    "2p2h_shape_O"         : -1.0,
    "2p2h_Edep_lowEnu"     :  0.0,
    "2p2h_Edep_highEnu"    :  0.0,
    "2p2h_Edep_lowEnubar"  :  0.0,
    "2p2h_Edep_highEnubar" :  0.0,
    "EB_dial_O_nu"         : -10.0,
    "EB_dial_O_nubar"      : -10.0
}

# Dictionary to keep track of how many times we skip drawing a spline comparison, along with the reason
skipped_draw = {
    "no_valor_spline_array"                  : 0,
    "etrue_edges_dont_align"                 : 0,
    "erec_edges_dont_align"                  : 0,
    "theta_edges_dont_align"                 : 0,
    "neither_valor_nor_ptheta_spline"        : 0,
    "no_ptheta_spline_but_flat_valor_spline" : 0,
    "no_ptheta_spline"                       : 0,
    "valor_MCT_bin_empty"                    : 0,
    "no_skip"                                : 0
}

# This is to consider the difference in weights between VALOR and p-theta, for all params, modes and bins.
# times them by the VALOR pre-param bin content and sum over bins. Then work out what the effect of the
# differences are on the rates to see if these are indeed causing the differences I see in the actual rate tables.
rates_valor_MCT_valor_weights  = {}
rates_valor_MCT_ptheta_weights = {}



#_________________________________________________________________________________
def main(config_file, tag, Erec_GeV, Etrue_GeV = -999, theta_deg = -999, sample_str = "FHC1Rmu",
         init_flav = "NuMu", final_flav = "NuMu", selected_param = "", match_edges = False):
    fit_elements = valor.CompoundFitElements()
    fit_elements.Read(util.CompleteConfigPath(config_file),"nosplines")
    model = fit_elements.FitModel()

    # Get a list of all spline params
    pool = valor.ParameterPool.Instance()
    filtered_params = [p for p in pool.GetParams() if p.GetType() == conv.kSystResponse]
    assert(len(filtered_params) == len(valor_params_in_ptheta_order))

    # Get the sample we're interested in
    submodel = None
    for imodel in range(0, model.GetNSubModels()):
        submodel = model.GetSubModel(imodel)
        if valor_sample_full_names[sample_str] == submodel.Name():
            break
    if not submodel:
        print "ERROR: Could not get submodel named \'{}\'".format(sample_str)
    beam = submodel.GetBeam()
    detector = submodel.GetDetector()
    sample = submodel.GetSample()
    sskb = submodel.GetBinning()

    valor_erec_bin_edges = valor_erec_bin_edges_mulike
    valor_theta_bin_edges = None
    if sample_is_2d[sample_str]:
        valor_erec_bin_edges  = valor_erec_bin_edges_elike
        valor_theta_bin_edges = valor_theta_bin_edges_elike


    # Write the VALOR MC templates into an array so that we can later check the nominal bin content before
    # systematic or oscillation parameters are applied
    array_MCTs = ROOT.TObjArray()
    array_MCTs.SetOwner(True)
    array_tag = "valor"
    submodel.Write(array_MCTs, array_tag, "mct_vf+SB")

    consider_all_erec_bins  = (Erec_GeV  == -999)
    consider_all_true_bins  = (Etrue_GeV == -999)
    consider_all_theta_bins = (theta_deg == -999 and sample_is_2d[sample_str])

    save_name = tag +"_" +sample_str +"_" +init_flav +"_to_" +final_flav
    if not consider_all_true_bins:
        save_name += "_Etrue=" +str(Etrue_GeV) +"GeV"
    if not consider_all_erec_bins:
        save_name += "_Erec=" +str(Erec_GeV) +"GeV"
    if sample_is_2d[sample_str]:
        if not consider_all_theta_bins:
            save_name += "_theta=" +str(theta_deg) +"degrees"
    if match_edges:
        save_name += "_matching_Etrue_edges_only"
    save_name += ".pdf"

    # Just get the reaction modes that correspond to the osc channel we're looking at.
    # Also, if we're looking at an oscillated channel, discard the NC reac types from the above hardcoded reac type list
    # (as these aren't tracked in VALOR)
    valor_modes_in_ptheta_order = valor_modes_in_ptheta_order_all
    if init_flav != final_flav:
        valor_modes_in_ptheta_order = [m for m in valor_modes_in_ptheta_order_all if "NC" not in m]
    modelist = submodel.GetModeList().Modes()
    filtered_modelist = [m for m in modelist
                         if (util.ReacModeToInitFlavor (m) == util.FlavourFromString(init_flav ) and
                             util.ReacModeToFinalFlavor(m) == util.FlavourFromString(final_flav) and
                             util.ReacModeToType(m) != conv.kmtCC1pi0)          # p-theta groups the CC1piC and CC1pi0 reaction types, while VALOR separates them. So just look at CC1piC
                        ]
    assert(len(filtered_modelist) == len(valor_modes_in_ptheta_order))

    valor_2d_tag = "_"
    if sample_is_2d[sample_str]:
        valor_2d_tag = "erectheta_"
    valor_filename  = valor_spline_file_tag +valor_2d_tag +util.BeamConfAsString(beam) +"_" +util.SampleAsString(sample) +".root"
    file_v = ROOT.TFile(valor_filename , "read")
    if not file_v.IsOpen():
        print "ERROR: Could not open", valor_filename
        sys.exit()
    file_banff = ROOT.TFile(banff_filename , "read")
    if not file_banff.IsOpen():
        print "ERROR: Could not open", banff_filename
        sys.exit()

    first_plot = True
    n_different_splines = 0
    for ip, param in enumerate(filtered_params):
        param_name = param.GetName()
        rates_valor_MCT_valor_weights [param_name] = {}
        rates_valor_MCT_ptheta_weights[param_name] = {}

        if selected_param != "" and selected_param not in param_name:
            continue

        ptheta_param_enum = valor_params_in_ptheta_order.index(param_name) +1

        banff_param_index = list(file_banff.Get("param_list")).index(param_name)
        prefit_sigma  = ROOT.TMath.Sqrt(file_banff.Get("prefit_cov")[banff_param_index][banff_param_index])
        postfit_sigma = ROOT.TMath.Sqrt(file_banff.Get("postfit_cov")[banff_param_index][banff_param_index])
        postfit_nom   = list(file_banff.Get("postfit_params"))[banff_param_index]
        untuned       = list(file_banff.Get("banff_untuned_central_values"))[banff_param_index]

        postfit_nom_spline_x = (postfit_nom - untuned) / prefit_sigma
        tweak_p3 =  postfit_nom + 3*postfit_sigma
        tweak_m3 =  postfit_nom - 3*postfit_sigma
        if param_name in special_case_tweak_p3:
            tweak_p3 = special_case_tweak_p3[param_name]
        if param_name in special_case_tweak_m3:
            tweak_m3 = special_case_tweak_m3[param_name]
        postfit_p3_spline_x  = (tweak_p3 - untuned) / prefit_sigma
        postfit_m3_spline_x  = (tweak_m3 - untuned) / prefit_sigma

        xscale_ptheta = ptheta_xscales[param_name]

        ptheta_filename = ptheta_spline_file_tag +valor_sample_to_ptheta_sample[sample_str] \
                          +"_binning/" + get_ptheta_osc_channel(init_flav, final_flav) +"_ntrw" + str(ptheta_param_enum) + ".root"
        file_p = ROOT.TFile(ptheta_filename, "read")
        if not file_p.IsOpen():
            print "ERROR: Could not open", ptheta_filename
            sys.exit()

        for im, m in enumerate(filtered_modelist):
            if not util.IsAppliedTo(param.GetApplyMask(), beam, sample, detector, m, conv.kntInclusive):
                continue

            reac_type = util.ReacModeToType(m)
            reac_type_name = util.ReacTypeAsString(reac_type)
            ptheta_mode_enum = valor_modes_in_ptheta_order.index(reac_type_name)
            rates_valor_MCT_valor_weights [param_name][reac_type_name] = 0
            rates_valor_MCT_ptheta_weights[param_name][reac_type_name] = 0

            valor_etrue_bin_edges = get_valor_etrue_bin_edges(beam, init_flav, final_flav)

            arr_title_v = "spl_" +util.ModeAsString(m) +"_" +param_name
            arr_v = file_v.Get(arr_title_v)
            if not arr_v:
                print "WARNING: Could not get VALOR array named", arr_title_v, "from file", valor_filename
                skipped_draw["no_valor_spline_array"] += 1
                continue

            # Get splines, convert to comparable format and plot
            for true_bin_v in range(0, len(valor_etrue_bin_edges)-1):
                if not consider_all_true_bins:
                    if not (Etrue_GeV >= valor_etrue_bin_edges[true_bin_v] and Etrue_GeV < valor_etrue_bin_edges[true_bin_v+1]):
                        continue

                true_bin_centre_v = 0.5 * (valor_etrue_bin_edges[true_bin_v] + valor_etrue_bin_edges[true_bin_v+1])
                true_bin_p  = get_bin_number(ptheta_etrue_bin_edges, true_bin_centre_v)

                if match_edges and (valor_etrue_bin_edges[true_bin_v] != ptheta_etrue_bin_edges[true_bin_p] or
                                    valor_etrue_bin_edges[true_bin_v+1] != ptheta_etrue_bin_edges[true_bin_p+1]):
                    skipped_draw["etrue_edges_dont_align"] += 1
                    continue

                # Loop over Erec bins, if needed
                nbins_erec = 1
                if consider_all_erec_bins:
                    nbins_erec = len(valor_erec_bin_edges) -1
                for erec_bin_v in range(0, nbins_erec):
                    if not consider_all_erec_bins:
                        erec_bin_v = get_bin_number(valor_erec_bin_edges, Erec_GeV)
                    else:
                        Erec_GeV = 0.5 * (valor_erec_bin_edges[erec_bin_v] + valor_erec_bin_edges[erec_bin_v+1])
                    reco_bin_p  = get_bin_number(ptheta_erec_bin_edges , Erec_GeV)

                    if match_edges and (valor_erec_bin_edges[erec_bin_v]   != ptheta_erec_bin_edges[reco_bin_p] or
                                        valor_erec_bin_edges[erec_bin_v+1] != ptheta_erec_bin_edges[reco_bin_p+1]):
                        skipped_draw["erec_edges_dont_align"] += 1
                        continue

                    # Loop over theta bins, if needed
                    nbins_theta = 1
                    if consider_all_theta_bins:
                        nbins_theta = len(valor_theta_bin_edges_elike) -1

                    for theta_bin_v in range(0, nbins_theta):
                        canv_name = "c_" +param_name +"_" +reac_type_name +"_" +str(true_bin_v) +\
                                    "_" +str(erec_bin_v) +"_" +str(theta_bin_v)
                        canv = ROOT.TCanvas(canv_name, canv_name, 800, 600)
                        if not consider_all_theta_bins:
                            theta_bin_v = get_bin_number(valor_theta_bin_edges, theta_deg)
                        else:
                            theta_deg = 0.5 * (valor_theta_bin_edges_elike[theta_bin_v] + valor_theta_bin_edges_elike[theta_bin_v+1])
                        theta_bin_p = get_bin_number(ptheta_theta_bin_edges, theta_deg)

                        if match_edges and (valor_theta_bin_edges[theta_bin_v]   != ptheta_theta_bin_edges[theta_bin_p] or
                                            valor_theta_bin_edges[theta_bin_v+1] != ptheta_theta_bin_edges[theta_bin_p+1]):
                            skipped_draw["theta_edges_dont_align"] += 1
                            continue

                        reco_bin_v = None
                        if sample_is_2d:
                            reco_bin_v = get_valor_bin_number_2D(valor_erec_bin_edges, valor_theta_bin_edges, Erec_GeV, theta_deg)
                        else:
                            reco_bin_v = get_bin_number(valor_erec_bin_edges , Erec_GeV)
                        spl_title_v = "spl_" +util.BeamConfAsString(beam) +"_" +util.SampleAsString(sample) +"_" \
                                      +util.ModeAsString(m) +"_" +param_name +"_" +str(true_bin_v) +"_" +str(reco_bin_v)
                        #print "Getting VALOR spline named", spl_title_v

                        # By convention, the p-theta Erec and theta bin numbers start at 1 for the TGraph names (but their Etrue bins start at 0)
                        gr_title_p = "gr_rw_erec_im" + str(ptheta_mode_enum) + "_enu" +str(true_bin_p)+ "_binx" +str(reco_bin_p+1)
                        if sample_is_2d[sample_str]:
                            gr_title_p = "gr_rw_etheta_im" + str(ptheta_mode_enum) + "_enu" +str(true_bin_p)+ "_binx" +str(reco_bin_p+1) +"_biny" +str(theta_bin_p+1)
                        #print "Getting p-theta spline named", gr_title_p

                        tg_p = file_p.Get(gr_title_p)
                        if not tg_p:
                            spl_v = arr_v.FindObject(spl_title_v)
                            if not spl_v:
                                print "WARNING: could not get p-theta tgraph", gr_title_p, "from file", file_p.GetName(), \
                                      "\n  Couldn't get the associated VALOR spline either \n\n"
                                skipped_draw["neither_valor_nor_ptheta_spline"] += 1
                            else:
                                if is_spline_flat(spl_v):
                                    skipped_draw["no_ptheta_spline_but_flat_valor_spline"] += 1
                                else:
                                    print "WARNING: could not get p-theta tgraph", gr_title_p, "from file", file_p.GetName(), \
                                          "\n  The associated VALOR spline is NOT flat, so we should be be getting a p-theta spline \n\n"
                                    skipped_draw["no_ptheta_spline"]
                            continue
                        new_tg_p = rescale_and_shift_ptheta_spline(tg_p, prefit_sigma, xscale_ptheta)

                        spl_v = arr_v.FindObject(spl_title_v)
                        if not spl_v:
                            print "ERROR: Could not get VALOR TSpline3 named", spl_title_v, "from file", valor_filename
                            sys.exit()
                        tg_v = convert_valor_spline_to_ptheta_graph(spl_v, prefit_sigma)
                        tg_v.SetLineColor(2)
                        tg_v.SetLineStyle(2)

                        new_tg_p.SetTitle("p-theta")
                        tg_v.SetTitle("valor")

                        mg_name = "mg_" +str(true_bin_v) +"_" +str(reco_bin_v)
                        mg = ROOT.TMultiGraph(mg_name, mg_name)
                        mg.Add(deepcopy(new_tg_p),"LP")
                        mg.Add(deepcopy(tg_v),"LP")
                        if len(list(mg.GetListOfGraphs())) != 2:
                            print "ERROR: Something went wrong constructing the multigraph. p-theta graph =", new_tg_p, "; valor graph =", tg_v
                            sys.exit()
                        mg.GetXaxis().SetTitle("(param - generated) / (prefit sigma)")
                        mg.GetYaxis().SetTitle("Weight")

                        flav_str = ""
                        if init_flav == final_flav:
                            flav_str = init_flav
                        else:
                            flav_str = init_flav +"_to_" +final_flav

                        title_box = ROOT.TPaveText(0.3, 0.95, 0.7, 1.0, "NDC")
                        title = sample_str +"_" +flav_str +"_" +reac_type_name+ "_" +param_name+ "_" +str(true_bin_v)+ "_" +str(reco_bin_v)
                        title_box.AddText(title)
                        title_box.SetTextSize(0.05)
                        title_box.SetFillColor(0)
                        title_box.SetFillStyle(0)
                        title_box.SetLineStyle(0)
                        title_box.SetLineColor(0)
                        title_box.SetBorderSize(0)

                        # Get the VALOR nominal bin content before systematic or oscillation parameters are applied
                        MCT_name = array_tag +"_mctnom_sb_vf" +"_" +submodel.Name() +"_" +util.ModeAsString(m)
                        MCT = array_MCTs.FindObject(MCT_name)
                        if not MCT:
                            print "ERROR: Something went wrong getting the VALOR MCT named", MCT_name
                            sys.exit()
                        valor_fit_binning_reco = sskb.RecoSmpBinning()
                        valor_fit_binning_true = sskb.TrueSmpBinning(m)
                        vec_reco = ROOT.vector('double')()
                        vec_reco.push_back(Erec_GeV)
                        if sample_is_2d[sample_str]:
                            vec_reco.push_back(theta_deg)
                        valor_fit_ibin_reco = valor_fit_binning_reco.FindBin(vec_reco)
                        # vec_true = ROOT.vector('double')()
                        # vec_true.push_back(true_bin_centre_v)
                        # valor_fit_ibin_true = valor_fit_binning_true.FindBin(vec_true)
                        valor_fit_ibin_true = get_bin_number(valor_etrue_fit_bin_edges, true_bin_centre_v)
                        valor_bin_content = MCT.GetEntry(valor_fit_ibin_reco, valor_fit_ibin_true)

                        # There's no point in comparing splines that would be applied to empty bins, as the resulting bin content would still be 0
                        if valor_bin_content < conv.kASmallNum:
                            skipped_draw["valor_MCT_bin_empty"] += 1
                            continue

                        # Count the effect of differences in this spline for this param, if all other params are not applied
                        wght_v = tg_v    .Eval(postfit_nom_spline_x, 0, "S")
                        wght_p = new_tg_p.Eval(postfit_nom_spline_x, 0, "S")
                        rates_valor_MCT_valor_weights [param_name][reac_type_name] += valor_bin_content * wght_v
                        rates_valor_MCT_ptheta_weights[param_name][reac_type_name] += valor_bin_content * wght_p
                        if ROOT.TMath.Abs(wght_v - wght_p) > 1E-5:
                            n_different_splines += 1
                            print "WARNING: Spline difference > 1E-5 at", title

                        ptheta_bin_txt = "p-theta: Et = ["+str(ptheta_etrue_bin_edges[true_bin_p])+","+str(ptheta_etrue_bin_edges[true_bin_p+1])+"]; Er = ["+str(ptheta_erec_bin_edges[reco_bin_p])+","+str(ptheta_erec_bin_edges[reco_bin_p+1])+"]"
                        valor_bin_txt  = "VALOR:   Et = ["+str(valor_etrue_bin_edges[true_bin_v]) +","+str(valor_etrue_bin_edges[true_bin_v+1]) +"]; Er = ["+str(valor_erec_bin_edges[erec_bin_v]) +","+str(valor_erec_bin_edges[erec_bin_v+1]) +"]"
                        if sample_is_2d[sample_str]:
                            ptheta_bin_txt += "; th = ["+str(ptheta_theta_bin_edges[theta_bin_p])+","+str(ptheta_theta_bin_edges[theta_bin_p+1])+"]"
                            valor_bin_txt  += "; th = ["+str(valor_theta_bin_edges[theta_bin_v]) +","+str(valor_theta_bin_edges[theta_bin_v+1]) +"]"
                        valor_content_txt = "Pre-param bin content = {:.1e}".format(valor_bin_content)
                        binning_box = ROOT.TPaveText(0.30, 0.75, 0.90, 0.95, "NDC")
                        binning_box.AddText(ptheta_bin_txt)
                        binning_box.AddText(valor_bin_txt)
                        binning_box.AddText(valor_content_txt)
                        binning_box.SetTextSize(0.05)
                        binning_box.SetFillColor(0)
                        binning_box.SetFillStyle(0)
                        binning_box.SetLineStyle(0)
                        binning_box.SetLineColor(0)
                        binning_box.SetBorderSize(0)

                        # Draw a shaded box to show the postfit nom +/- 3sigma region
                        min_y = min(min(new_tg_p.GetY()), min(tg_v.GetY()))
                        max_y = max(max(new_tg_p.GetY()), max(tg_v.GetY()))
                        hatch_interval  = ROOT.TBox(postfit_m3_spline_x , min_y, postfit_p3_spline_x, max_y)
                        hatch_interval.SetLineWidth(0)
                        hatch_interval.SetFillStyle(3003)
                        hatch_interval.SetFillColor(9)

                        # Draw a vertical line to show the postfit nominal value
                        postfit_nom_line = ROOT.TLine(postfit_nom_spline_x, min_y, postfit_nom_spline_x, max_y)
                        postfit_nom_line.SetLineColor(8)
                        postfit_nom_line.SetLineWidth(2)
                        postfit_nom_line.SetLineStyle(2)

                        leg = ROOT.TLegend(0.13, 0.80, 0.28, 0.95)
                        leg.AddEntry(new_tg_p        , "p-theta"     , "lp")
                        leg.AddEntry(tg_v            , "VALOR"       , "lp")
                        leg.AddEntry(postfit_nom_line, "Postfit nom.", "l" )
                        leg.AddEntry(hatch_interval  , "+/- 3#sigma" , "f" )
                        leg.SetTextFont(132)
                        leg.SetTextSize(0.03)

                        mg.Draw("A")
                        hatch_interval.Draw("same")
                        postfit_nom_line.Draw("same")
                        title_box.Draw()
                        binning_box.Draw()
                        leg.Draw()

                        save_name_use = save_name
                        if first_plot:
                            save_name_use += "("
                            first_plot = False
                        canv.Print(save_name_use)
                        canv.Clear()
                        skipped_draw["no_skip"] += 1
                        del canv
                        del mg

    # Close the PDF (this prints a blank page before closing, but is easier than looping over everything before the
    # main loop to determine the final used combination of ip, im and true_bin_v)
    canv_blank = ROOT.TCanvas("blank","blank",800,600)
    canv_blank.Print(save_name +")")

    if match_edges and consider_all_true_bins and consider_all_theta_bins:
        # Here we only want to consider bins where the Etrue edges match. This is so that the valor MCT
        # bin content will be ~= to that of p-theta, which would not be the case otherwise.
        # This allows us to see what rates differences are definitely caused by spline differences.
        print "Spline contributions to differences in", sample_str, init_flav, "-->", final_flav, \
              " event rates for single param+mode combinatons", \
              "(using VALOR MCT bin contents, but weights from both groups) (expressed as valor_rate/ptheta_rate -1):"
        for im, m in enumerate(filtered_modelist):
            reac_type = util.ReacModeToType(m)
            reac_type_name = util.ReacTypeAsString(reac_type)
            print "\t", reac_type_name
            for ip, param in enumerate(filtered_params):
                if not util.IsAppliedTo(param.GetApplyMask(), beam, sample, detector, m, conv.kntInclusive):
                    continue
                param_name = param.GetName()

                rv = rates_valor_MCT_valor_weights [param_name][reac_type_name]
                rp = rates_valor_MCT_ptheta_weights[param_name][reac_type_name]
                ratiovp = (rv - rp) *100.
                if rv != 0 and rp != 0:
                    ratiovp = (rv/rp -1) *100.
                print "\t\t", param_name, ":", ratiovp, "%"


    print "\n\nOut of", sum(skipped_draw.values()), "attempted spline comparions", sum(skipped_draw.values()) - skipped_draw["no_skip"], \
          "were not drawn, for the following reasons:"
    if skipped_draw["no_valor_spline_array"]:
        print skipped_draw["no_valor_spline_array"], "where the VALOR TObjarray (containing splines for a given sample+mode+param combination) could not be found. This shouldn't happen."
    if skipped_draw["etrue_edges_dont_align"]:
        print skipped_draw["etrue_edges_dont_align"], "where the VALOR and p-theta Etrue bin edges didn't align. This is fine if intentional."
    if skipped_draw["erec_edges_dont_align"]:
        print skipped_draw["erec_edges_dont_align"], "where the VALOR and p-theta Erec bin edges didn't align. This is fine if intentional."
    if skipped_draw["theta_edges_dont_align"]:
        print skipped_draw["theta_edges_dont_align"], "where the VALOR and p-theta theta bin edges didn't align. This is fine if intentional."
    if skipped_draw["neither_valor_nor_ptheta_spline"]:
        print skipped_draw["neither_valor_nor_ptheta_spline"], "where neither the VALOR nor the p-theta spline could be found. This shouldn't happen."
    if skipped_draw["no_ptheta_spline_but_flat_valor_spline"]:
        print skipped_draw["no_ptheta_spline_but_flat_valor_spline"], "where the p-theta spline couldn't be found, but the associated VALOR spline was flat. This is fine, as p-theta doesn't store flat splines."
    if skipped_draw["no_ptheta_spline"]:
        print skipped_draw["no_ptheta_spline"], "where the p-theta spline couldn't be found, but the associated VALOR spline exists and is NOT flat. This shouldn't happen."
    if skipped_draw["valor_MCT_bin_empty"]:
        print skipped_draw["valor_MCT_bin_empty"], "where the VALOR MC template bin was empty. This is fine, as the spline wouldn't affect this bin anyway."

    print "\nNumber of splines with different (> 1E-5) VALOR/p-theta weights at postfit nominal =", n_different_splines

    print "\n\nSaving output to file \'{}\'".format(save_name)
    print "Done!"

#_________________________________________________________________________________
def rescale_and_shift_ptheta_spline(tgraph, prefit_sigma, ptheta_scale):
    new_tgraph = ROOT.TGraph(tgraph.GetN())
    for ip in range(0, tgraph.GetN()):
        x = ROOT.Double(-999)
        y = ROOT.Double(-999)
        tgraph.GetPoint(ip,x,y)
        new_x = x * ptheta_scale / prefit_sigma
        new_y = y + 1.0
        new_tgraph.SetPoint(ip, new_x, new_y)
    return new_tgraph

#_________________________________________________________________________________
def convert_valor_spline_to_ptheta_graph(tspline, prefit_sigma):
    new_tgraph = ROOT.TGraph(tspline.GetNp())
    for ip in range(0, tspline.GetNp()):
        x = ROOT.Double(-999)
        y = ROOT.Double(-999)
        tspline.GetKnot(ip, x, y)
        new_x = x / prefit_sigma
        new_tgraph.SetPoint(ip, new_x, y)
    return new_tgraph

#_________________________________________________________________________________
def get_bin_number(edges, value):
    for i in range(0, len(edges)-1):
        if value >= edges[i] and value < edges[i+1]:
            return i
    print "ERROR: could not get bin number for x =", value, "in edges:", edges

#_________________________________________________________________________________
def get_valor_bin_number_2D(edges_erec, edges_theta, value_erec, value_theta):
    e = get_bin_number(edges_erec , value_erec )
    t = get_bin_number(edges_theta, value_theta)
    return (len(edges_erec)-1)*t + e


#_________________________________________________________________________________
def get_ptheta_osc_channel(init_flav, final_flav):
    if init_flav == final_flav:
        if init_flav == "NuMu":
            return "numu"
        elif init_flav == "NuE":
            return "nue"
        elif init_flav == "NuMuBar":
            return "numubar"
        elif init_flav == "NuEBar":
            return "nuebar"
    else:
        if init_flav == "NuMu" and final_flav == "NuE":
            return "nue_sig"
        elif init_flav == "NuMuBar" and final_flav == "NuEBar":
            return "nuebar_sig"

    print "Could not get the p-theta osc channel for init flav =", init_flav, ", final flav =", final_flav
    sys.exit()


#_________________________________________________________________________________
def get_valor_etrue_bin_edges(beam, init_flav, final_flav):
    rightsign = True
    if beam == conv.kbJPARC_FHC:
        if "Bar" in init_flav or "Bar" in final_flav:
            rightsign = False
    else:
        if "Bar" not in init_flav or "Bar" not in final_flav:
            rightsign = False

    numu_and_osc = True
    if "NuMu" in final_flav or ("NuMu" in init_flav and "NuE" in final_flav):
        numu_and_osc = True
    elif "NuE" in init_flav and "NuE" in final_flav:
        numu_and_osc = False
    else:
        print "ERROR: could not parse osc channel combination init_flav =", init_flav, ", final_flav =", final_flav
        sys.exit()

    if rightsign and numu_and_osc:
        return valor_etrue_bin_edges_rightsign_numu_and_osc_nue
    elif rightsign and not numu_and_osc:
        return valor_etrue_bin_edges_rightsign_nue
    elif not rightsign and numu_and_osc:
        return valor_etrue_bin_edges_wrongsign_numu_and_osc_nue
    elif not rightsign and not numu_and_osc:
        return valor_etrue_bin_edges_wrongsign_nue
    else:
        print "ERROR: Shouldn't be here!"
        sys.exit()

#________________________________________________________________________________________
def is_spline_flat(tspline):
    assert( type(tspline) == ROOT.TSpline3 )
    nk = tspline.GetNp()
    for ik in range(0, nk):
        kx = ROOT.Double(-999)
        ky = ROOT.Double(-999)
        tspline.GetKnot(ik, kx, ky)
        if ROOT.TMath.Abs(ky - 1.0) > 1E-10:
            return False
    return True

#_________________________________________________________________________________
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="%prog [options]..",description="Compare 1D Erec splines for VALOR and p-theta. Erec is fixed, Etrue is looped over. FHC 1Rmu mumu->numu only.")
    parser.add_argument('-c', '--config'        , dest='config'        , default="analysis/t2k_analysis_configuration_oa2018_run1-10_v1_postfit.xml", help="Configuration file")
    parser.add_argument(      '--tag'           , dest='tag'           , default="spline_comparison_valor_ptheta_JPARC_FHC_1rmu_numu", help="tag used to name output file")
    parser.add_argument('-r', '--Erec-GeV'      , dest="Erec_GeV"      , default='-999', help="Erec in GeV. If specified, just looks at one Erec bin, otherwise looks at all.")
    parser.add_argument('-a', '--theta-deg'     , dest="theta_deg"     , default='-999', help="theta rec in degrees. If specified, just looks at one theta bin, otherwise looks at all.")
    parser.add_argument('-t', '--Etrue-GeV'     , dest="Etrue_GeV"     , default='-999' , help="Etrue in GeV. Optional. If specified, just looks at one Erec bin, otherwise looks at all.")
    parser.add_argument('-s', '--sample'        , dest="sample"        , default='FHC1Rmu', choices=["FHC1Rmu","RHC1Rmu","FHC1Re","RHC1Re","FHCMultiRe"])
    parser.add_argument('-i', '--init-flav'     , dest="init_flav"     , default='NuMu', choices=["NuMu","NuE","NuMuBar","NuEBar"])
    parser.add_argument('-f', '--final-flav'    , dest="final_flav"    , default='NuMu', choices=["NuMu","NuE","NuMuBar","NuEBar"])
    parser.add_argument('-p', '--selected-param', dest="selected_param", default='', help="Write a param name here to only look at the splines for that param. Otherwise will look at all params.")
    parser.add_argument('-m', '--match-edges'   , dest="match_edges"   , default=False, action="store_true", help="If true, will only compare splines when the VALOR and p-theta Etrue bin edges are matching.")
    args = parser.parse_args()
    main(args.config, args.tag, float(args.Erec_GeV), float(args.Etrue_GeV), float(args.theta_deg), args.sample,
         args.init_flav, args.final_flav, args.selected_param, args.match_edges)
 
