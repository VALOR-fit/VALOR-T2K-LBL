#!/usr/bin/python

import ROOT
import math
import array
import sys
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

env_var = "T2KOSCANA"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()
arridx= util.ArrIdx.Instance()

samples   = ["1Rmu","1Re","RHC1Rmu","RHC1Re","MultiRmu","MultiRe"]

col1_width = 50
col_width = 20

line_width = 1

to3dp = "{0:.3f}"
to6dp = "{0:.6f}"

#______________________________________________________________________________________________________________________________

def main(sample,output_file,config_file,config_file_noflux,config_file_prefit):

  
  nominal_modes        = dict()
  NA61_flux_modes      = dict()
  prefit_tuned_modes   = dict()
  postfit_tuned_modes  = dict()
  
  osc_nominal_modes        = dict()
  osc_NA61_flux_modes      = dict()
  osc_prefit_tuned_modes   = dict()
  osc_postfit_tuned_modes  = dict()
  
  all_dicts_unosc = [nominal_modes, NA61_flux_modes, prefit_tuned_modes, postfit_tuned_modes]
  all_dicts_osc   = [osc_nominal_modes, osc_NA61_flux_modes, osc_prefit_tuned_modes, osc_postfit_tuned_modes]
  
  min_energy_spectrum = 0.0
  min_energy_ratio    = 0.
  max_energy = 1.25

  
  #  Make no flux SKExp
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file_noflux)
  
  #~ skexp.GetBANFF().SetTweaksToNominal()
  # skexp.GetBANFF().RestoreNominalMC()
  # skexp.GetBANFF().RestoreRHCNominalMC()
  #~ skexp.GetBANFF().RestoreFHCpRHCNominalMC()
  skexp.GetBANFF().Tweak(0,0,"untuned")
  
  skexp.GetSuperkMC().FillTemplateArrays(False)
  
  fparm_unosc,fparm_osc = make_fit_param_lists()
  print fparm_osc.GetParListName()
  skexp.Generate(fparm_unosc)
  nominal_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  nominal_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    nominal_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
    
  skexp.Generate(fparm_osc)
  osc_nominal_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  osc_nominal_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    osc_nominal_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
  
  del skexp
  
  
  # Now do the NA61 tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  #~ skexp.GetBANFF().SetTweaksToNominal()
  # skexp.GetBANFF().RestoreNominalMC()
  # skexp.GetBANFF().RestoreRHCNominalMC()
  #~ skexp.GetBANFF().RestoreFHCpRHCNominalMC()
  skexp.GetBANFF().Tweak(0,0,"untuned")
  
  skexp.GetSuperkMC().FillTemplateArrays(False)
  
  fparm_unosc,fparm_osc = make_fit_param_lists()
  
  skexp.Generate(fparm_unosc)
  NA61_flux_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  flux_tuned_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    NA61_flux_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
    
  skexp.Generate(fparm_osc)
  osc_NA61_flux_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  osc_flux_tuned_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    osc_NA61_flux_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral()
  
  del skexp
  
  # NIWG tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file_prefit)
  fparm_unosc,fparm_osc = make_fit_param_lists()
  
  skexp.Generate(fparm_unosc)
  prefit_tuned_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    prefit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
    
  skexp.Generate(fparm_osc)
  osc_prefit_tuned_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    osc_prefit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
  
  del skexp
  
  # postfit tune
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ[env_var]+"/config/"+config_file)
  
  fparm_unosc,fparm_osc = make_fit_param_lists()
  
  skexp.Generate(fparm_unosc)
  postfit_tuned_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  postfit_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    postfit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
    
  skexp.Generate(fparm_osc)
  osc_postfit_tuned_modes["total"] = skexp.CurrPredictedErecoVsQ2VArray(sample).Integral()
  osc_postfit_hist = hist_to_50MeV(skexp.CurrPredictedErecoSpectrum(sample))
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    osc_postfit_tuned_modes[m] = skexp.CurrPredictedErecoSpectrum(sample,m).Integral("")
  
  del skexp
  
  # Start writing the TeX
  tex_output = open(output_file,"w")
  
  tex_output.write("""\documentclass[10pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\usepackage{geometry}

\\begin{document}

""")
  tex_output.write("\\begin{table}\n")
  tex_output.write(util.AsString(sample))
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*4)+" | }\n\\hline\n")
  
  cols = ["","$N_{SK}$","$N_{SK}$","$N_{SK}$","$N_{SK}$"]
  widths = [col1_width]+([col_width]*(len(cols)-1))
  print widths
    
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  
  cols = ["","Nominal MC","Flux-tuned MC","BANFF prefit MC","BANFF postfit MC"]
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  cols = ["Total"] + [to6dp.format(d["total"]) for d in all_dicts_unosc]
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  tex_output.write("\\hline \n")
  
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    cols = [util.ModeAsStringTex(m)] + [to6dp.format(d[m]) for d in all_dicts_unosc]
    entries = [x.ljust(w) for x,w in zip(cols,widths)]
    tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  tex_output.write("\n\\end{tabular}")
  tex_output.write("\n\\caption{Unoscillated}\n")
  tex_output.write("\\end{table}\n\n")
  
  # Oscillated
  tex_output.write("\\begin{table}\n")
  #~ tex_output.write("\\centering\n")
  tex_output.write("\n\n\\begin{tabular} { | c ||"+("c | "*4)+" | }\n\\hline\n")
  
  cols = ["","$N_{SK}$","$N_{SK}$","$N_{SK}$","$N_{SK}$"]
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \n")
  
  cols = ["","Nominal MC","Flux-tuned MC","BANFF prefit MC","BANFF postfit MC"]
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  cols = ["Total"] + [to6dp.format(d["total"]) for d in all_dicts_osc]
  entries = [x.ljust(w) for x,w in zip(cols,widths)]
  tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  tex_output.write("\\hline \n")
  
  for im in xrange(conv.knm):
    m = util.GetMode(sample,im)
    if m == conv.kSKmUndefined or m > conv.knm: continue
    cols = [util.ModeAsStringTex(m)] + [to6dp.format(d[m]) for d in all_dicts_osc]
    entries = [x.ljust(w) for x,w in zip(cols,widths)]
    tex_output.write("& ".join(entries) + " \\\\ \\hline \n")
  
  tex_output.write("\n\\end{tabular}")
  tex_output.write("\n\\caption{Oscillated}\n")
  tex_output.write("\\end{table}\n\n")
  
  
  tex_output.write("\\end{document}\n")
  
  print "Table output written to ",output_file
  
  
  
  # Now make some figures
  
  line_at_one = ROOT.TLine(min_energy_ratio,1,max_energy,1)
  line_at_one.SetLineColor(ROOT.kBlue)
  
  # First - flux tuning effect compared to nominal
  canv = ROOT.TCanvas("c","c",1280,1024)
  canv.Divide(2,2)
  
  canv.cd(1)
  if nominal_hist.GetXaxis().GetXmax() > max_energy:
    nominal_hist.GetXaxis().SetRangeUser(min_energy_spectrum,max_energy)
  nominal_hist.SetLineColor(ROOT.kBlue)
  nominal_hist.SetLineWidth(line_width)
  nominal_hist.SetMaximum(1.1*max([nominal_hist.GetMaximum(),flux_tuned_hist.GetMaximum()]))
  nominal_hist.SetTitle("")
  nominal_hist.Draw()
  flux_tuned_hist.SetLineColor(ROOT.kRed)
  flux_tuned_hist.SetLineWidth(line_width)
  flux_tuned_hist.SetTitle("")
  flux_tuned_hist.Draw("same")
  #fluxnom_leg = ROOT.TLegend(0.4,0.6,0.9,0.9)
  fluxnom_leg = ROOT.TLegend(0.7,0.7,1.0,1.0)
  fluxnom_leg.AddEntry(nominal_hist,"Nominal")
  #fluxnom_leg.AddEntry(flux_tuned_hist,"Flux tuning applied")
  fluxnom_leg.AddEntry(flux_tuned_hist,"Flux tuning NA61")
  fluxnom_leg.Draw()
  
  canv.cd(2)
  if osc_nominal_hist.GetXaxis().GetXmax() > max_energy:
    osc_nominal_hist.GetXaxis().SetRangeUser(min_energy_spectrum,max_energy)
  osc_nominal_hist.SetLineColor(ROOT.kBlue)
  osc_nominal_hist.SetLineWidth(line_width)
  osc_nominal_hist.SetMaximum(1.1*max([osc_nominal_hist.GetMaximum(),osc_flux_tuned_hist.GetMaximum()]))
  osc_nominal_hist.SetTitle("")
  osc_nominal_hist.Draw()
  osc_flux_tuned_hist.SetLineColor(ROOT.kRed)
  osc_flux_tuned_hist.SetLineWidth(line_width)
  osc_flux_tuned_hist.SetTitle("")
  osc_flux_tuned_hist.Draw("same")
  fluxnom_leg.Draw()
  
  canv.cd(3)
  ratio = ROOT.TH1D(flux_tuned_hist)
  ratio.Divide(nominal_hist)
  ratio.GetYaxis().SetTitle("Flux tuned/nominal")
  ratio.SetLineColor(ROOT.kBlack)
  if ratio.GetXaxis().GetXmax() > max_energy:
    ratio.GetXaxis().SetRangeUser(min_energy_ratio,max_energy)
  ratio_dev = max([abs(1.-x) for x in hist_to_list(ratio)[5:]])
  ratio.SetMinimum(1-1.1*ratio_dev)
  ratio.SetMaximum(1+1.1*ratio_dev)
  ratio.Draw()
  line_at_one.Draw("same")
  
  canv.cd(4)
  osc_ratio = ROOT.TH1D(osc_flux_tuned_hist)
  osc_ratio.Divide(osc_nominal_hist)
  osc_ratio.GetYaxis().SetTitle("Flux tuned/nominal")
  osc_ratio.SetLineColor(ROOT.kBlack)
  if osc_ratio.GetXaxis().GetXmax() > max_energy:
    osc_ratio.GetXaxis().SetRangeUser(min_energy_ratio,max_energy)
  osc_ratio_dev = max([abs(1.-x) for x in hist_to_list(osc_ratio)[5:]])
  osc_ratio.SetMinimum(1-1.1*osc_ratio_dev)
  osc_ratio.SetMaximum(1+1.1*osc_ratio_dev)
  osc_ratio.Draw()
  line_at_one.Draw("same")
  
  canv.cd(0)
  canv.Print("flux_tuning_effect.eps")
  
  del canv
  
  # Second - BANFF tuning effect compared to flux tuned
  canv = ROOT.TCanvas("cc","cc",1280,1024)
  canv.Divide(2,2)
  
  pf_ratio = ROOT.TH1D(postfit_hist)                 # do these copies up here before we set the limits
  pf_ratio.Divide(flux_tuned_hist)
  osc_pf_ratio = ROOT.TH1D(osc_postfit_hist)
  osc_pf_ratio.Divide(osc_flux_tuned_hist)
  
  canv.cd(1)
  if postfit_hist.GetXaxis().GetXmax() > max_energy:
    postfit_hist.GetXaxis().SetRangeUser(min_energy_spectrum,max_energy)
  postfit_hist.SetLineColor(ROOT.kBlack)
  postfit_hist.SetLineWidth(line_width)
  postfit_hist.SetMaximum(1.1*max([postfit_hist.GetMaximum(),flux_tuned_hist.GetMaximum()]))
  postfit_hist.SetTitle("")
  postfit_hist.Draw()
  flux_tuned_hist.SetLineColor(ROOT.kRed)
  flux_tuned_hist.SetLineWidth(line_width)
  flux_tuned_hist.SetTitle("")
  flux_tuned_hist.Draw("same")
  #postfit_leg = ROOT.TLegend(0.4,0.6,0.9,0.9)
  postfit_leg = ROOT.TLegend(0.7,0.7,1.0,1.0)
  #postfit_leg.AddEntry(flux_tuned_hist,"Flux tuning applied")
  postfit_leg.AddEntry(flux_tuned_hist,"Flux tuning NA61")
  #postfit_leg.AddEntry(postfit_hist,"BANFF tune applied")
  postfit_leg.AddEntry(postfit_hist,"BANFF tune")
  postfit_leg.Draw()
  
  canv.cd(2)
  if osc_postfit_hist.GetXaxis().GetXmax() > max_energy:
    osc_postfit_hist.GetXaxis().SetRangeUser(min_energy_spectrum,max_energy)
  osc_postfit_hist.SetLineColor(ROOT.kBlack)
  osc_postfit_hist.SetLineWidth(line_width)
  osc_postfit_hist.SetMaximum(1.1*max([osc_postfit_hist.GetMaximum(),osc_flux_tuned_hist.GetMaximum()]))
  osc_postfit_hist.SetTitle("")
  osc_postfit_hist.Draw()
  osc_flux_tuned_hist.SetLineColor(ROOT.kRed)
  osc_flux_tuned_hist.SetLineWidth(line_width)
  osc_flux_tuned_hist.SetTitle("")
  osc_flux_tuned_hist.Draw("same")
  postfit_leg.Draw()
  
  canv.cd(3)
  pf_ratio.GetYaxis().SetTitle("BANFF tuned / flux tuned")
  pf_ratio.SetLineColor(ROOT.kBlack)
  pf_ratio.SetLineWidth(line_width)
  pf_ratio.SetTitle("")
  if pf_ratio.GetXaxis().GetXmax() > max_energy:
    pf_ratio.GetXaxis().SetRangeUser(min_energy_ratio,max_energy)
  pf_ratio_dev = max([abs(1.-x) for x in hist_to_list(pf_ratio)[5:]])
  pf_ratio.SetMinimum(1-1.1*pf_ratio_dev)
  pf_ratio.SetMaximum(1+1.1*pf_ratio_dev)
  pf_ratio.GetYaxis().SetRangeUser(0.8,1.2)
  pf_ratio.Draw()
  line_at_one.Draw("same")
  
  canv.cd(4)
  osc_pf_ratio.GetYaxis().SetTitle("BANFF tuned / flux tuned")
  osc_pf_ratio.SetLineColor(ROOT.kBlack)
  osc_pf_ratio.SetLineWidth(line_width)
  osc_pf_ratio.SetTitle("")
  if osc_pf_ratio.GetXaxis().GetXmax() > max_energy:
    osc_pf_ratio.GetXaxis().SetRangeUser(min_energy_ratio,max_energy)
  osc_pf_ratio_dev = max([abs(1.-x) for x in hist_to_list(osc_pf_ratio)[5:]])
  osc_pf_ratio.SetMinimum(1-1.2*osc_pf_ratio_dev)
  osc_pf_ratio.SetMaximum(1+1.2*osc_pf_ratio_dev)
  osc_pf_ratio.GetYaxis().SetRangeUser(0.8,1.2)
  osc_pf_ratio.Draw()
  line_at_one.Draw("same")
  
  canv.cd(0)
  canv.Print("banff_tuning_effect.eps")
  
  del canv
  

def hist_to_50MeV(h):
  hout = ROOT.TH1D(h)
  for ix in xrange(1,hout.GetNbinsX()+1):
    hout.SetBinContent(ix, 0.05 * h.GetBinContent(ix) / h.GetBinWidth(ix));
  hout.GetYaxis().SetTitle("Events per 0.050 GeV");
  return hout
  
def  make_fit_param_lists():
    
  #mixing12 = opt.GetGlobalBestFit(conv.kPhypMixing12)
  #mixing13 = conv.kT2KBestFitTheta13_NH
  #mixing23 = opt.GetGlobalBestFit(conv.kPhypMixing23)
  #dm12sq   = opt.GetGlobalBestFit(conv.kPhypDm12sq)
  #dm23sq   = opt.GetGlobalBestFit(conv.kPhypDm23sq)
  #deltacp  = conv.kT2KBestFitDeltaCP_NH
  
  # fparm_osc = util.NamedParameterList("nuosc3f_numudisapp_fix4+nuisance")
  fparm_osc = util.NamedParameterList("nuosc3f+nuisance")
  
  fparm_osc.SetInitValue(conv.kPhypMixing12,  0.306);
  fparm_osc.SetInitValue(conv.kPhypMixing13,  0.0217);
  fparm_osc.SetInitValue(conv.kPhypMixing23,  0.528);
  fparm_osc.SetInitValue(conv.kPhypCPV     ,  -1.601);
  fparm_osc.SetInitValue(conv.kPhypDm12sq  ,  7.5E-5);
  fparm_osc.SetInitValue(conv.kPhypDm23sq  ,  2.509E-3);
  
  # fparm_unosc = util.NamedParameterList("nuosc3f_numudisapp_fix4+nuisance")
  fparm_unosc = util.NamedParameterList("nuosc3f+nuisance")
  
  fparm_unosc.SetInitValue(conv.kPhypMixing12,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypMixing13,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypMixing23,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypCPV     ,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypDm12sq  ,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypDm23sq  ,  0.0)
  
  return fparm_unosc,fparm_osc
  
def hist_to_list(h):
  nbins = h.GetXaxis().GetNbins()
  return [h.GetBinContent(ibin) for ibin in xrange(1,nbins+1)]
  
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of applied systematics")  
  parser.add_option('-f', '--file', dest='output_file', default="effects_of_tuning.tex", help="Name for output TeX file.")
  
  parser.add_option('-s', '--sample' , dest='sample', default="1Rmu", help="Sample (1Rmu, 1Re, RHC1Rmu, RHC1Re, MultiRmu, MultiRe)")
  
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  parser.add_option('-p', '--config-prefit' , dest='pconfig', default="analysis_configuration_prefit.xml", help="Prefit config file")
  parser.add_option('-n', '--config-no-flux', dest='fconfig', default="analysis_configuration_nofluxtune.xml", help="Configuration file without flux tune")
  
  
  o,other = parser.parse_args()
  sample = util.SKSampleFromString(o.sample)
  main(sample,o.output_file, o.config, o.fconfig, o.pconfig)
  
  if other:
    print "Warning - ignored arguments: ",other
