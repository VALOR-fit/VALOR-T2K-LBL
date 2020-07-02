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

# Get p-theta mode enum with e.g. valor_modes_in_ptheta_order.index("CCcoherent")
valor_modes_in_ptheta_order = [
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

# Get p-theta param enum with e.g. valor_params_in_ptheta_order.index("2p2h_shape_O")
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

# Etrue edges specifically for FHC 1Rmu mumu->mumu
ptheta_etrue_bin_edges = [0.0, 0.4, 0.5, 0.6, 0.8, 1.2, 2.0, 3.0, 5.0, 30.0]
valor_etrue_bin_edges  = [0.0, 0.4,	0.5, 0.6, 0.7, 1.0, 1.5, 2.5, 3.5, 5.0, 7.0, 30.0]

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
valor_erec_bin_edges = [
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



#_________________________________________________________________________________
def main(param, reac_type, Erec_GeV):
    # fit_elements = valor.CompoundFitElements()
    # fit_elements.Read(util.CompleteConfigPath(config_file))
    # model = fit_elements.FitModel()

    ptheta_param_enum = valor_params_in_ptheta_order.index(param) +1
    ptheta_mode_enum = valor_modes_in_ptheta_order.index(reac_type)

    ptheta_filename = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/t2k_run1-10/new_sk_matrix_banff_postfit_data_20200409/syst_var/other_groups/p-theta_splines/RF/numode_numu_binning/numu_ntrw" + str(ptheta_param_enum) + ".root"
    valor_filename  = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/xsec/splines2020v6/valor_format/splines_spline2019splinesfitqunPre_JPARC_FHC_1rmu.root"
    banff_filename  = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/data/official/2020/banff/postfit_data_20200409/BANFF_postfit_OA2020_fixed2p2h_fixedIHalfLowPPi_20200409_valor_order_all.root"

    file_p = ROOT.TFile(ptheta_filename, "read")
    if not file_p.IsOpen():
        print "ERROR: Could not open", ptheta_filename
        sys.exit()
    file_v = ROOT.TFile(valor_filename , "read")
    if not file_v.IsOpen():
        print "ERROR: Could not open", valor_filename
        sys.exit()
    file_banff = ROOT.TFile(banff_filename , "read")
    if not file_banff.IsOpen():
        print "ERROR: Could not open", banff_filename
        sys.exit()

    banff_param_index = list(file_banff.Get("param_list")).index(param)
    prefit_sigma  = ROOT.TMath.Sqrt(file_banff.Get("prefit_cov")[banff_param_index][banff_param_index])
    postfit_sigma = ROOT.TMath.Sqrt(file_banff.Get("postfit_cov")[banff_param_index][banff_param_index])
    postfit_nom   = list(file_banff.Get("postfit_params"))[banff_param_index]
    untuned       = list(file_banff.Get("banff_untuned_central_values"))[banff_param_index]

    postfit_nom_spline_x = (postfit_nom - untuned) / prefit_sigma
    postfit_p3_spline_x  = (postfit_nom + 3*postfit_sigma - untuned) / prefit_sigma
    postfit_m3_spline_x  = (postfit_nom - 3*postfit_sigma - untuned) / prefit_sigma

    xscale_ptheta = ptheta_xscales[param]

    # first loop over the input p-theta file to see what range of Etrue bins we have splines for (when VALOR would have a flat spline, p-theta stores no spline)
    first_etrue_bin = -999
    last_etrue_bin = -999
    for true_bin_v in range(0, len(valor_etrue_bin_edges)-1):
        true_bin_centre_v = 0.5 * (valor_etrue_bin_edges[true_bin_v] + valor_etrue_bin_edges[true_bin_v+1])
        true_bin_p = get_bin_number(ptheta_etrue_bin_edges, true_bin_centre_v)
        reco_bin_p = get_bin_number(ptheta_erec_bin_edges , Erec_GeV)
        gr_title_p = "gr_rw_erec_im" + str(ptheta_mode_enum) + "_enu" +str(true_bin_p)+ "_binx" +str(reco_bin_p)
        tg_p = file_p.Get(gr_title_p)
        if tg_p:
            if first_etrue_bin == -999:
                first_etrue_bin = true_bin_v
            last_etrue_bin = true_bin_v

    # Get splines, convert to comparable format and plot
    for true_bin_v in range(0, len(valor_etrue_bin_edges)-1):
        canv = ROOT.TCanvas("c+"+str(true_bin_v),"c+"+str(true_bin_v),800,600)

        true_bin_centre_v = 0.5 * (valor_etrue_bin_edges[true_bin_v] + valor_etrue_bin_edges[true_bin_v+1])
        reco_bin_v = get_bin_number(valor_erec_bin_edges , Erec_GeV)
        arr_title_v = "spl_numu" +reac_type+ "_" + param
        spl_title_v = "spl_JPARC_FHC_1rmu_numu" +reac_type+ "_" +param+ "_" +str(true_bin_v)+ "_" +str(reco_bin_v)
        print "Getting VALOR spline named", spl_title_v

        true_bin_p = get_bin_number(ptheta_etrue_bin_edges, true_bin_centre_v)
        reco_bin_p = get_bin_number(ptheta_erec_bin_edges , Erec_GeV)
        gr_title_p = "gr_rw_erec_im" + str(ptheta_mode_enum) + "_enu" +str(true_bin_p)+ "_binx" +str(reco_bin_p)
        print "Getting p-theta spline named", gr_title_p

        tg_p = file_p.Get(gr_title_p)
        if not tg_p:
            print "WARNING: Could not get p-theta TGraph named", gr_title_p, "from file", ptheta_filename, ". Skipping this Etrue bin."
            continue
        new_tg_p = rescale_and_shift_ptheta_spline(tg_p, prefit_sigma, xscale_ptheta)

        spl_v = file_v.Get(arr_title_v).FindObject(spl_title_v)
        if not spl_v:
            print "ERROR: Could not get VALOR TSpline3 named", spl_title_v, "from file", valor_filename
            sys.exit()
        tg_v = convert_valor_spline_to_ptheta_graph(spl_v, prefit_sigma)
        tg_v.SetLineColor(2)
        tg_v.SetLineStyle(2)

        new_tg_p.SetTitle("p-theta")
        tg_v.SetTitle("valor")

        mg = ROOT.TMultiGraph("mg_"+str(true_bin_v), "mg_"+str(true_bin_v))
        mg.Add(deepcopy(new_tg_p),"LP")
        mg.Add(deepcopy(tg_v),"LP")
        if len(list(mg.GetListOfGraphs())) != 2:
            print "ERROR: Something went wrong constructing the multigraph. p-theta graph =", new_tg_p, "; valor graph =", tg_v
            sys.exit()
        mg.GetXaxis().SetTitle("(param - generated) / (prefit sigma)")
        mg.GetYaxis().SetTitle("Weight")

        title_box = ROOT.TPaveText(0.3, 0.95, 0.7, 1.0, "NDC")
        title_box.AddText(spl_title_v)
        title_box.SetTextSize(0.05)
        title_box.SetFillColor(0)
        title_box.SetFillStyle(0)
        title_box.SetLineStyle(0)
        title_box.SetLineColor(0)
        title_box.SetBorderSize(0)

        ptheta_bin_txt = "p-theta: Et = ["+str(ptheta_etrue_bin_edges[true_bin_p])+","+str(ptheta_etrue_bin_edges[true_bin_p+1])+"]; Er = ["+str(ptheta_erec_bin_edges[reco_bin_p])+","+str(ptheta_erec_bin_edges[reco_bin_p+1])+"]"
        valor_bin_txt  = "VALOR:   Et = ["+str(valor_etrue_bin_edges[true_bin_v]) +","+str(valor_etrue_bin_edges[true_bin_v+1]) +"]; Er = ["+str(valor_erec_bin_edges[reco_bin_v]) +","+str(valor_erec_bin_edges[reco_bin_v+1]) +"]"
        binning_box = ROOT.TPaveText(0.30, 0.80, 0.90, 0.95, "NDC")
        binning_box.AddText(ptheta_bin_txt)
        binning_box.AddText(valor_bin_txt)
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
        save_name = "spline_comparison_valor_ptheta_JPARC_FHC_1rmu_numu" +reac_type+ "_" +param+ "_Erec=" +str(Erec_GeV)+ "GeV.pdf"
        if true_bin_v == first_etrue_bin:
            save_name += "("
        if true_bin_v == last_etrue_bin:
            save_name += ")"
        canv.Print(save_name)

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
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(usage="%prog [options]..",description="Compare 1D Erec splines for VALOR and p-theta. Erec is fixed, Etrue is looped over. FHC 1Rmu mumu->numu only.")
    parser.add_argument('-p', '--valor-param', dest="param"    , default='MAQE' , help="VALOR name of the parameter to compare splines for.")
    parser.add_argument('-r', '--reac-type'  , dest="reac_type", default='CCQE' , help="VALOR name of the reaction type to compare splines for.")
    parser.add_argument('-e', '--Erec-GeV'   , dest="Erec_GeV" , default='0.555', help="Erec in GeV")
    args = parser.parse_args()
    main(args.param, args.reac_type, float(args.Erec_GeV))
 
