#!/usr/bin/env python

# Makes a data release from the T2K legacy analysis fit files

import ROOT
import sys

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()


tag = "fit_oa2018"
top_dir = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/t2k_run1-9d/run1-9d_final_fits_from_Andy_PDG2018_smeared/"
vars_list = ["mixing23_dmsq23", "mixing13_deltacp", "mixing23_deltacp", "deltacp"]
outvars = {
    "mixing23_dmsq23"  : "th23_dm2",
    "mixing13_deltacp" : "th13_dcp",
    "mixing23_deltacp" : "th23_dcp",
    "deltacp"          : "dcp"
}
filenames = {
    "mixing23_dmsq23"  : "fit_oa2018_[MO]_[REAC]_mixing23_dmsq23_data_smeared_contour.root",
    "mixing13_deltacp" : "fit_oa2018_[MO]_[REAC]_mixing13_deltacp_nexp_1_0.022_-1.601-job0_-_0_0_contour.root",
    "mixing23_deltacp" : "fit_oa2018_[MO]_[REAC]_mixing23_deltacp_nexp_1_0.528_-1.601-job0_-_0_0_contour.root",
    "deltacp"          : "fit_oa2018_[MO]_[REAC]_deltacp_nexp_1_-1.601_-999.000-job0_-_0_0_contour.root",
}

def make_data_release():
    util.plotting.SetupCustomPalette("valor_inverted_dark_body_radiator")

    fout = ROOT.TFile("data_release_analysis_a.root", "recreate")

    for use_reac in [True, False]:
        for vars in vars_list:
            for mo in ["normal", "inverted"]:
                add_objects(fout, vars, use_reac, mo)

    fout.Close()
    print "Done!"


#________________________________________________________________________________________________________
def add_objects(fout, vars, use_RC, mo):
    ROOT.TGaxis.SetMaxDigits(3)

    is_1D = ("_" not in vars)

    reac = "reactor"
    if not use_RC:
        reac = "t2k"

    out_ndim = "2D"
    if is_1D:
        out_ndim = "1D"
    out_reac = "woRC"
    if use_RC:
        out_reac = "wRC"
    out_mo = "NO"
    if mo == "inverted":
        out_mo = "IO"

    outname_hist = "A_{}_dchi2_{}_{}_{}"        .format(out_ndim, outvars[vars], out_reac, out_mo)
    outname_bf   = "A_{}_dchi2_{}_{}_{}_BestFit".format(out_ndim, outvars[vars], out_reac, out_mo)

    fname = top_dir +"/" +vars +"/data/" +filenames[vars]
    fname = fname.replace("[MO]"  , mo  )
    fname = fname.replace("[REAC]", reac)
    f = ROOT.TFile(fname, "read")
    if not f or not f.IsOpen():
        print "ERROR: Could not get file \'{}\'".format(f.GetName())
        sys.exit()

    # Get chisq hist
    hname = None
    if is_1D:
        hname = "h_deltachisq"
    else:
        hname = "dchisq_map_2D"
    h = f.Get(hname)
    if not h or type(h) == ROOT.TGraph:
        print "ERROR: Could not get hist \'{}\'".format(hname), "from file \'{}\'".format(f.GetName())
        sys.exit()

    if vars == "mixing23_dmsq23":
        # The smeared hists are far too large and need rebinning to make the file size acceptable
        h = rebin_smeared(h, outname_hist)

    if not is_1D:
        # Get rid of the stupid formatting from the legacy analysis
        h = reformat_2d_hist(h, outname_hist)
    else:
        h.SetName (outname_hist)
        h.SetTitle(outname_hist)

    # Get BF marker
    bf = None
    if not is_1D:
        bf = f.Get("best_fit_point_marker")
        if not bf or type(bf) == ROOT.TObject:
            print "ERROR: Could not get marker \'best_fit_point_marker\' from file \'{}\'".format(f.GetName())
            sys.exit()

    # Get confidence regions
    crs = []
    conf_levels = ["68","90"]
    if not is_1D:
        for cl in conf_levels:
            for icont in range(0,10):
                old_name = "g"+cl+"_valor_"+str(icont)
                new_name = "A_2D_dchi2_{}_{}_{}_{}pct_{}".format(outvars[vars], out_reac, out_mo, cl, icont)
                cr = f.Get(old_name)
                if not cr or type(cr) != ROOT.TGraph:
                    continue
                cr.SetName (new_name)
                cr.SetTitle(new_name)
                crs.append(cr)

    # Get FC critical values
    crit_graphs = []
    fc_conf_levels = ["68", "90", "955", "997"]
    old_fc_conf_names = {
        "68"  : "1sigma",
        "90"  : "90CL",
        "955" : "2sigma",
        "997" : "3sigma"
    }
    old_fc_mo_names = {
        "normal"   : "nh",
        "inverted" : "ih"
    }
    if vars == "deltacp" and use_RC:
        fcritname = "/hepstore/fbench/software/valor/VALOR-T2K-LBL/t2k_run1-9d/oa_long_paper/fc/legacy/FCgraphs.root"
        fcrit = ROOT.TFile(fcritname, "read")
        for cl in fc_conf_levels:
            old_name = "gFC{}_{}".format(old_fc_conf_names[cl], old_fc_mo_names[mo])
            new_name = "A_1D_dchi2_dcp_{}_{}_{}conf".format(out_reac, out_mo, cl)
            crit = fcrit.Get(old_name)
            if not crit or type(crit) == ROOT.TObject:
                print "ERROR: Could not get FC crit graph \'{}\' from file \'{}\'".format(old_name, f.GetName())
                sys.exit()
            crit.SetName (new_name)
            crit.SetTitle(new_name)
            crit_graphs.append(crit)
        fcrit.Close()

    fout.cd()
    h.Write(outname_hist)
    print "Saving", outname_hist
    if not is_1D:
        bf.Write(outname_bf)
        print "Saving", outname_bf
        for cr in crs:
            cr.Write()
            print "Saving", cr.GetName()
    if vars == "deltacp" and use_RC:
        for crit in crit_graphs:
            crit.Write()
            print "Saving", crit.GetName()

    f.Close()

#________________________________________________________________________________________________________
def reformat_2d_hist(h_old, outname_hist):
    nx = h_old.GetNbinsX()
    ny = h_old.GetNbinsY()
    lx = h_old.GetXaxis().GetBinLowEdge(1)
    ly = h_old.GetYaxis().GetBinLowEdge(1)
    ux = h_old.GetXaxis().GetBinUpEdge(nx)
    uy = h_old.GetYaxis().GetBinUpEdge(ny)
    h_new = ROOT.TH2D(outname_hist, outname_hist, nx, lx, ux, ny, ly, uy)
    assert(h_new)
    h_new.GetXaxis().SetTitle(h_old.GetXaxis().GetTitle())
    h_new.GetYaxis().SetTitle(h_old.GetYaxis().GetTitle())
    h_new.GetZaxis().SetTitle(h_old.GetZaxis().GetTitle())
    for ix in range(1, 1+h_old.GetNbinsX()):
        for iy in range(1, 1+h_old.GetNbinsY()):
            fill = h_old.GetBinContent(ix,iy)
            h_new.SetBinContent(ix,iy,fill)
    assert(h_new.Integral() == h_old.Integral())
    h_new.SetDirectory(0)
    return h_new

#________________________________________________________________________________________________________
def rebin_smeared(old_h, outname_hist):
    h = old_h.Rebin2D(9,9,outname_hist)
    h.SetDirectory(0)
    return h

#________________________________________________________________________________________________________
if __name__=="__main__":
    make_data_release()

