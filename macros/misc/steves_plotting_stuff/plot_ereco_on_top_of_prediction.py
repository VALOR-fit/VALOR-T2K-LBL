#!/usr/bin/env python
import ROOT
import os
import math

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

to2dp = "{0:.2f}"
to3dp = "{0:.3f}"
to4sf = '{0:.3e}'

orig_plot_opt = ""

text_entry_height = 0.08

def main(file_patterns,tag="data"):
  
  max_energy = 7
  
  sample = conv.kSKsRHCFCFV1RingMuLike
  
  objects_to_keep = []
  
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ['T2KOSCANA']+"/config/analysis_configuration_2.315E20.xml")
  
  fparm_unosc = util.NamedParameterList("nuosc3f_numudisapp_fix4+nuisance")
  fparm_unosc.SetInitValue(conv.kPhypMixing12,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypMixing13,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypMixing23,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypCPV     ,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypDm12sq  ,  0.0)
  fparm_unosc.SetInitValue(conv.kPhypDm23sq  ,  0.0)
  skexp.Generate(fparm_unosc);
  unosc_spectrum = skexp.CurrPredictedErecoSpectrum(sample);
  
  unosc_spectrum.GetXaxis().SetRangeUser(0,0.5)
  n_unosc_low_energy  = unosc_spectrum.Integral()
  unosc_spectrum.GetXaxis().SetRangeUser(0.5,1.0)
  n_unosc_mid_energy  = unosc_spectrum.Integral()
  unosc_spectrum.GetXaxis().SetRangeUser(1.0,30)
  n_unosc_high_energy = unosc_spectrum.Integral()
  
  
  skobs = valor.SuperkObservation()
  skobs.Read(os.environ['T2KOSCANA']+"/inputs/data.xml")
  
  data_spectrum = skobs.ErecoSpectrum(sample)
  data_spectrum.SetLineColor(ROOT.kBlack)
  data_spectrum.SetLineWidth(2)
  
  n_data = data_spectrum.Integral()
  
  data_spectrum.GetXaxis().SetRangeUser(0,0.5)
  n_data_low_energy  = data_spectrum.Integral()
  data_spectrum.GetXaxis().SetRangeUser(0.5,1.0)
  n_data_mid_energy  = data_spectrum.Integral()
  data_spectrum.GetXaxis().SetRangeUser(1.0,30)
  n_data_high_energy = data_spectrum.Integral()
  
  
  data_spectrum.GetXaxis().SetRangeUser(0,max_energy)
  
  assert(abs(data_spectrum.Integral()-n_data)<conv.kASmallNum)
  
  data_max = data_spectrum.GetMaximum()
  data_max = data_max+math.sqrt(data_max)
  
  ch = ROOT.TChain("fitnt")
  for pattern in file_patterns:
    ch.Add(pattern)
  assert(ch.GetEntries())
  ch.GetEntry(0)
  
  canvas = ROOT.TCanvas("c","c",1280,720)
  
  for ientry,entry in enumerate(ch):
    fparm = entry.best_fit_params
    
    # Set up hierarchy
    if ("31" in fparm.Name(conv.kPhypDm23sq)):
      opt.gOscDmsqScheme = conv.kNuDmsq12andAbs31
    elif ("32" in fparm.Name(conv.kPhypDm23sq)):
      opt.gOscDmsqScheme = conv.kNuDmsq12andAbs23
    else:
      print "Something went wrong, probably unsupported Fogli-Lisi"
      exit(-1)
    
    skexp.Generate(fparm)
    histogram = skexp.CurrPredictedErecoSpectrum(sample)
    histogram.SetMaximum(1.1*max(data_max,histogram.GetMaximum()))
    histogram.SetLineColor(ROOT.kRed)
    histogram.SetLineWidth(2)
    n_pred = histogram.Integral()
    
    canvas.Clear()
    
    legend = ROOT.TLegend(0.7,0.7,0.89,0.89)
    legend.AddEntry(data_spectrum,"Data")
    legend.AddEntry(histogram,"Best-fit")
    
    pave = ROOT.TPaveText(0.5,0.56,0.68,0.87,"ndc")
    pave.SetShadowColor(0)
    pave.SetFillColor(0)
    pave.SetFillStyle(0)
    pave.SetLineColor(0)
    pave.AddText("N_{data} = "+`ROOT.TMath.Nint(n_data)`)
    pave.AddText("N_{bf} = "+to3dp.format(histogram.Integral()))
    
    pave2 = ROOT.TPaveText(0.71,0.63,0.89,0.68,"ndc")
    pave2.SetShadowColor(0)
    pave2.SetFillColor(0)
    pave2.SetFillStyle(0)
    pave2.SetLineColor(0)
    pave2.AddText("#chi^{2} = "+to3dp.format(fparm.GetBestFitStatistic()))
    
    pave3 = ROOT.TPaveText(0.65,0.40,0.89,0.54,"ndc")
    pave3.SetShadowColor(0)
    pave3.SetFillColor(0)
    pave3.SetFillStyle(0)
    pave3.SetLineColor(0)
    pave3.AddText(util.ParamAsStringPaperRooTex(conv.kPhypMixing23Bar)+" = "+to3dp.format(fparm.BestFitValue(conv.kPhypMixing23Bar)))
    pave3.AddText(util.ParamAsStringPaperRooTex(conv.kPhypDm23sqBar  )+" = "+to4sf.format(fparm.BestFitValue(conv.kPhypDm23sqBar))+" eV^{2}")
    
    data_spectrum.Draw("e")
    histogram.Draw("same")
    legend.Draw()
    pave.Draw()
    pave2.Draw()
    pave3.Draw()
    canvas.Print(tag+"_spectrum_"+`ientry`+".eps")
    
    
    histogram.GetXaxis().SetRangeUser(0,0.5)
    n_bf_low_energy  = histogram.Integral()
    histogram.GetXaxis().SetRangeUser(0.5,1.0)
    n_bf_mid_energy  = histogram.Integral()
    histogram.GetXaxis().SetRangeUser(1.0,30)
    n_bf_high_energy = histogram.Integral()
  
  
    print "N predicted: ",n_pred
    print "N data: ",ROOT.TMath.Nint(n_data)
    print
  
  print "Low  E ratio: ",to2dp.format(100.*n_data_low_energy /n_unosc_low_energy ),"%"
  print "Mid  E ratio: ",to2dp.format(100.*n_data_mid_energy /n_unosc_mid_energy ),"%"
  print "High E ratio: ",to2dp.format(100.*n_data_high_energy/n_unosc_high_energy),"%"

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] pattern1 [pattern2] ...",description="Make plots.")
  parser.add_option('-t', dest='tag',default="data"  ,help='The filename')
  
  options,patterns = parser.parse_args()
  
  main(patterns,options.tag)
