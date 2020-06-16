#!/usr/bin/env python

# A bunch of utilities used for testing the T2K run 10 spline files

import ROOT

#________________________________________________________________________________________
def print_BANFF_untuned_postfit_comparison(banff_filename, spline_filename):
    if ('elements' not in locals() or 'elements' not in globals()) and util.NameRegistry.Instance().PhysParamFromString("flux_numu00") == valor.conv.kPhypUndefined:
        elements = valor.CompoundFitElements()
        elements.Read(util.CompleteConfigPath("analysis/t2k_analysis_configuration_oa2018_run1-10_v1_postfit.xml"), "nomc,nosplines")
    pool = valor.ParameterPool.Instance()
    nr   =  util.NameRegistry.Instance()

    banff_file = ROOT.TFile(banff_filename, "read")         # This is the VALOR-format BANFF file
    param_list = banff_file.Get("param_list")
    banff_untuned_central_values = banff_file.Get("banff_untuned_central_values")
    postfit_params = banff_file.Get("postfit_params")
    prefit_cov = banff_file.Get("prefit_cov")
    postfit_cov = banff_file.Get("postfit_cov")

    spline_file = ROOT.TFile(spline_filename, "read")       # This is a VALOR-format spline file

    print "param \t\t untuned nom \t postfit nom \t prefit sig \t postfit sig \t spline limits \t\t (postnom-untuned)/presig"
    for ip in range(0,param_list.GetEntries()):
        name = str(param_list.At(ip))
        p = nr.PhysParamFromString(name)
        param = pool.GetParam(p)
        lims = "  -  "
        if param.GetType() == valor.conv.kSystResponse:
            spl = get_first_nonflat_valor_splines_for_param(spline_file, name)
            l = [spl.GetXmin()+banff_untuned_central_values[ip], spl.GetXmax()+banff_untuned_central_values[ip]]
            lims = "[{:.2f}, {:.2f}]".format(l[0],l[1])
        print name, "\t {:.2f}".format(banff_untuned_central_values[ip]), "\t\t {:.2f}".format(postfit_params[ip]), "\t\t {:.2f}".format(ROOT.TMath.Sqrt(prefit_cov[ip][ip])), "\t\t {:.2f}".format(ROOT.TMath.Sqrt(postfit_cov[ip][ip])), "\t\t", lims, "\t\t {:.2f}".format(
            (postfit_params[ip] - banff_untuned_central_values[ip])/ROOT.TMath.Sqrt(prefit_cov[ip][ip]))

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
def get_first_nonflat_valor_splines_for_param(tfile, param_name):
    for k in tfile.GetListOfKeys():
        name = k.GetName()
        if (not "spl_" in name) or (not param_name in name):
            continue
        arr = k.ReadObj()
        for spl in arr:
            if not is_spline_flat(spl):
                return spl
