#!/usr/bin/env python2

import ROOT
import os
import sys

ROOT.PyConfig.IgnoreCommandLineOptions = True

base_dir = os.environ["VALOR_FMWK"]
try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x %s/root/valor_interpreter.C"%base_dir);

ROOT.gROOT.SetBatch(1)

valor = ROOT.valor
util  = valor.util
conv  = valor.conv
opt   = valor.Opt.Instance()
env   = valor.Env.Instance()
rnd   = valor.Rnd.Instance()

#______________________________________________________________________________________________________________________________

def banff_table_header(tex):
    tex.write("\\begin{table}[htp] \n")
    tex.write("\\centering \n")
    tex.write("\\begin{tabular} { c c @{\\quad}c@{\\quad} c c }\n\\hline\\hline \n")
    tex.write("\multirow{2}{*}{Index} & \multirow{2}{*}{Parameter} & \multirow{2}{*}{Description} & \multirow{2}{*}{Best fit} & $1\sigma$ pre/postfit\\\\ \n")
    tex.write(" & &  &  & fractional error\\\\\n\\hline \n")

def sk_table_header(tex):
    tex.write("\\begin{table}[htp] \n")
    tex.write("\\centering \n")
    tex.write("\\begin{tabular} { c c @{\\quad}c@{\\quad} c }\n\\hline\\hline \n")
    tex.write("Index & Parameter & Description & $1\sigma$ fractional error\\\\ \n\\hline \n")

def table_footer(tex, caption = "", label = ""):
    tex.write("\\hline\\hline\n")
    tex.write("\\end{tabular}\n")
    tex.write("\\caption{" +caption +"}\n")
    tex.write("\\label{" +label +"}\n")
    tex.write("\\end{table}\n")
    
def thesis_format_name(name):
    if "FHC and RHC" in name:
        name = name.replace("FHC and RHC", "\\FHCandRHCName\\")
    elif "FHC" in name:
        name = name.replace("FHC", "\\FHCName\\")
    elif "RHC" in name:
        name = name.replace("RHC", "\\RHCName\\")
    if "1R$_{\\mu}$" in name:
        name = name.replace("1R$_{\\mu}$", "\\MuLikeSampleName\\")
    if "1R$_e$ + 1d.e." in name:
        name = name.replace("1R$_e$ + 1d.e.", "\\MultiRingELikeSampleName\\")
    if "MultiRing " in name:
        name = name.replace("MultiRing", "\\MultiRingELikeSampleName\\")
    if "1R$_{e}$" in name:
        name = name.replace("1R$_{e}$", "\\ELikeSampleName\\")
    if "1R$_e$" in name:
        name = name.replace("1R$_e$", "\\ELikeSampleName\\")
    return name

def main(output_file, config, banff_filename, thesis_format):

    # get_prefit_errors(config_pre)

    elements = valor.CompoundFitElements()
    elements.Read(util.CompleteConfigPath(config), "nosplines nomc")
    model = elements.FitModel()
    syst = model.GetSystematics()
    if not syst.Size():
        print "You haven't got any postfit systematics defined. This is never going to work."
        exit(-1)

    pool = valor.ParameterPool.Instance()
    params = list(pool.GetParams())

    file_banff = ROOT.TFile(env.TopLevelDataDir() +"/" +banff_filename , "read")
    if not file_banff.IsOpen():
        print "ERROR: Could not open", banff_filename
        sys.exit()
    banff_param_list = list(file_banff.Get("param_list"))
    banff_param_postfit_vals = list(file_banff.Get("postfit_params"))
    banff_prefit_cov  = file_banff.Get("prefit_cov")
    banff_postfit_cov = file_banff.Get("postfit_cov")

    xsec_tex = open(output_file + "_xsec.tex","w")
    flux_fhc_tex = open(output_file + "_fhc_flux.tex","w")
    flux_rhc_tex = open(output_file + "_rhc_flux.tex","w")
    sk_tex   = open(output_file + "_sk.tex"  ,"w")
    banff_table_header(xsec_tex)
    banff_table_header(flux_fhc_tex)
    banff_table_header(flux_rhc_tex)
    sk_table_header   (sk_tex)

    first_param_index = -999   # The first params in the pool are osc params, but we want to index from the first syst param
    for ip, param in enumerate(params):
        if not util.IsSystematic(param.GetID()):
            continue

        if first_param_index < 0:
            first_param_index = ip

        param_name     = param.GetName()
        tex_name_short = param.GetRooTeXName().replace("#","\\")
        description    = param.GetTeXName()

        if thesis_format:
            tex_name_short = thesis_format_name(tex_name_short)
            description    = thesis_format_name(description)

        categories = list(param.GetCategories())
        if "BANFF" in categories:
            banff_param_index = banff_param_list.index(param_name)
            prefit_sigma  = ROOT.TMath.Sqrt(banff_prefit_cov[banff_param_index][banff_param_index])
            postfit_sigma = ROOT.TMath.Sqrt(banff_postfit_cov[banff_param_index][banff_param_index])
            postfit_nom   = banff_param_postfit_vals[banff_param_index]

            row = "%d & %s & %s & %.3f & %.3f / %.3f\\\\ \n" % (ip-first_param_index, tex_name_short, description, postfit_nom, prefit_sigma, postfit_sigma)
            if "flux" in categories:
                if "FHC" in description:
                    flux_fhc_tex.write(row)
                else:
                    flux_rhc_tex.write(row)
            elif "int" in categories:
                xsec_tex.write(row)
            else:
                print "Could not sub-categorise BANFF parameter \'{}\'".format(param_name)
                util.Die(-1)
        elif "SKDetFSI" in categories:
            sigma = syst.OneSigmaError(param.GetID())
            sk_tex.write("%d & %s & %s & %.3f\\\\ \n" % (ip-first_param_index, tex_name_short, description, sigma))
        else:
            print "Could not categorise parameter \'{}\'".format(param_name)
            util.Die(-1)

    table_footer(xsec_tex    , "Summary of cross section systematic parameters. Those with equal prefit and postfit errors were intentionally left unconstrained.", "tab:params_summary_xsec")
    table_footer(flux_fhc_tex, "Summary of $\\nu$-mode flux systematic normalisation parameters."                                                                 , "tab:params_summary_flux_numode")
    table_footer(flux_rhc_tex, "Summary of $\\bar{\\nu}$-mode flux systematic normalisation parameters."                                                          , "tab:params_summary_flux_nubarmode")
    table_footer(sk_tex      , "Summary of SK detector + SI + PN systematic normalisation parameters."                                                            , "tab:params_summary_sk")

    xsec_tex.close()
    flux_fhc_tex.close()
    flux_rhc_tex.close()
    sk_tex.close()

    print "Done!"

if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")
  
  parser.add_option('-f', '--file'          , dest='output_file'  , default="size_of_systematics"                                                                                                    , help="Name for output TeX file.")
  parser.add_option('-c', '--config-postfit', dest='config'       , default="analysis/t2k_analysis_configuration_oa2020_run1-10_v1_postfit.xml"                                                      , help="Configuration file")
  parser.add_option('-b', '--banff-file'    , dest='banff_file'   , default="official/2020/banff/postfit_data_20200409/BANFF_postfit_OA2020_fixed2p2h_fixedIHalfLowPPi_20200409_valor_order_all.root", help="File containing the prefit and postfit BANFF matrices")
  parser.add_option('-t', '--thesis-format' , dest='thesis_format', default=False, action="store_true"                                                                                               , help="If true, replaces the beam current and sample names with macros (because I've given them non-standard names in my thesis).")
  
  o,other = parser.parse_args()
  main(o.output_file, o.config, o.banff_file, o.thesis_format)
  
  if other:
    print "Warning - ignored arguments: ",other
