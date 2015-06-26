#!/usr/bin/python
import ROOT
import os
import math
import array

ROOT.PyConfig.IgnoreCommandLineOptions = True

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

coarse_binning = [0.0,0.4,0.7,1.0,1.5,2.0,3.0,4.0,5.0,30]

def main(file_patterns,config,tag="data",rebin=False,nopave=False,altlabels=True):
  
  max_energy = 5
  title_size = 0.08
  label_size = 0.07
  
  sample = conv.kSKsRHCFCFV1RingMuLike
  
  plot_type = "groups"
  if altlabels:
    plot_type = "bayes_vs_freq"
  
  objects_to_keep = []
  
  skexp = valor.SuperkExpectation()
  skexp.Read(os.environ['T2KOSCANA']+"/config/"+config)
  
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
  unosc_spectrum.GetXaxis().SetRangeUser(1.0,30.)
  n_unosc_high_energy = unosc_spectrum.Integral()
  unosc_spectrum.GetXaxis().SetRangeUser(0.,30.)
  
  
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
  
  if (rebin):
    data_spectrum = rebin_histogram_coarser(data_spectrum,coarse_binning)
    
  data_spectrum.GetXaxis().SetRangeUser(0,max_energy)
  
  assert(abs(data_spectrum.Integral()-n_data)<conv.kASmallNum)
  
  data_max = data_spectrum.GetMaximum()
  data_max = data_max+math.sqrt(data_max)
  
  data_spectrum.GetXaxis().SetTitleOffset(0.8)
  data_spectrum.GetYaxis().SetTitleOffset(0.6)
  data_spectrum.GetXaxis().SetTitleSize(title_size)
  data_spectrum.GetYaxis().SetTitleSize(title_size)
  data_spectrum.GetXaxis().SetLabelSize(label_size)
  data_spectrum.GetYaxis().SetLabelSize(label_size)
  
  ch = ROOT.TChain("fitnt")
  ch.Add(file_patterns[0])
  assert(ch.GetEntries())
  ch.GetEntry(0)
  
  canvas = ROOT.TCanvas("c","c",1280,720)
  SetBLMargins(canvas)
  SetRMargins (canvas,False)
  
  ch.GetEntry(0)
  
  fparm = ch.best_fit_params
  
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
  
  if (rebin): histogram = rebin_histogram_coarser(histogram,coarse_binning)
  histogram.GetXaxis().SetRangeUser(0,max_energy)
  
  canvas.Clear()
  
  # That's the VALOR contour done, let's get MaCh3
  m3_file = ROOT.TFile(file_patterns[1])
  assert(m3_file)
  mach3_hist = m3_file.Get("hist_1d")
  #~ mach3_hist = m3_file.Get("bestfit")
  assert(mach3_hist)
  mach3_hist.SetLineColor(ROOT.kBlue)
  mach3_hist.SetLineWidth(2)
  
  # And P-theta
  pt_file = ROOT.TFile(file_patterns[2])
  assert(pt_file)
  pthet_hist = pt_file.Get("hNumuPDF")
  #~ pthet_hist = pt_file.Get("coarse_binning_spectrum")
  assert(pthet_hist)
  pthet_hist.SetLineColor(ROOT.kGreen+2)
  pthet_hist.SetLineWidth(2)
  
  if rebin and pthet_hist.GetXaxis().GetNbins() > len(coarse_binning):
    pthet_hist = rebin_histogram_coarser(pthet_hist,coarse_binning)
  
  legend = ROOT.TLegend(0.68,0.55,0.93,0.93)
  legend.AddEntry(data_spectrum,"Data"    )
  if altlabels:
    legend.AddEntry(histogram    ,"Frequentist","l")
    legend.AddEntry(mach3_hist   ,"Bayesian"   ,"l")
  else:
    legend.AddEntry(histogram    ,"VALOR"   ,"l")
    legend.AddEntry(mach3_hist   ,"MaCh3"   ,"l")
    legend.AddEntry(pthet_hist   ,"P-theta" ,"l")
  
  legend.SetBorderSize(0)
  legend.SetFillStyle (0)
  

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
  pave3.AddText("#bar{sin^{2}#theta_{23}} = "+to3dp.format(fparm.BestFitValue(conv.kPhypMixing23Bar)))
  pave3.AddText("#bar{#Delta m^{2}_{23}} = " +to4sf.format(fparm.BestFitValue(conv.kPhypDm23sqBar))+" eV^{2}")
  
  data_spectrum.GetYaxis().SetTitle("Events per bin")
  
  data_spectrum.Draw("e")
  histogram .Draw("hist same")
  mach3_hist.Draw("hist same")
  if not altlabels:
    pthet_hist.Draw("hist same")
  histogram .Draw("hist same")
  data_spectrum.Draw("e same")
  legend.Draw()
  if not nopave:
    pave.Draw()
    pave2.Draw()
    pave3.Draw()
  canvas.Print(tag+"_spectrum_"+plot_type+".eps")
  
  
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
  
  data_spectrum.GetXaxis().SetRangeUser(0,30)
  histogram.GetXaxis().SetRangeUser(0,30)
  mach3_hist.GetXaxis().SetRangeUser(0,30)
  pthet_hist.GetXaxis().SetRangeUser(0,30)
  print
  print "Event count (data):    ",data_spectrum.Integral()
  print "Event count (VALOR):   ",histogram.Integral()
  print "Event count (MaCH3):   ",mach3_hist.Integral()
  print "Event count (p-theta): ",pthet_hist.Integral()
  
  # Plot per 100 MeV
  histogram.GetXaxis().SetRangeUser(0.,30)
  modify_histogram_to_events_per_width(data_spectrum,   0.5)
  modify_histogram_to_events_per_width(histogram    ,   0.5)
  modify_histogram_to_events_per_width(mach3_hist   ,   0.5)
  modify_histogram_to_events_per_width(pthet_hist   ,   0.5)
  data_spectrum.GetYaxis().SetTitle("Events per 0.5 GeV")
  data_spectrum.GetXaxis().SetRangeUser(0,max_energy)
  
  canvas2 = ROOT.TCanvas("c2","c2",1280,720)
  canvas2.cd()
  data_spectrum.Draw("e")
  histogram .Draw("hist same")
  mach3_hist.Draw("hist same")
  if not altlabels:
    pthet_hist.Draw("hist same")
  histogram .Draw("hist same")
  data_spectrum.Draw("e same")
  legend.Draw()
  if not nopave:
    pave.Draw()
    pave2.Draw()
    pave3.Draw()
  canvas.Print(tag+"_spectrum_"+plot_type+"_per_energy"+".eps")
  

def modify_histogram_to_events_per_width(h,width=0.1):
  for ibin in xrange(1,h.GetXaxis().GetNbins()+1):
    lower_edge = h.GetXaxis().GetBinLowEdge(ibin)
    upper_edge = h.GetXaxis().GetBinUpEdge (ibin)
    error = h.GetBinError(ibin)
    h.SetBinContent(ibin,h.GetBinContent(ibin) * width/(upper_edge-lower_edge))
    h.SetBinError(ibin, error * width/(upper_edge-lower_edge))
  return
  
def rebin_histogram_coarser(h,new_binning):
  arr = array.array('d',new_binning)
  new_hist = ROOT.TH1D(h.GetName(),h.GetTitle(),len(new_binning)-1,arr)
  
  for ibin in xrange(1,h.GetXaxis().GetNbins()+1):
    new_hist.Fill(h.GetBinCenter(ibin),h.GetBinContent(ibin))
  
  copy_hist_style(h,new_hist)
  new_hist.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
  new_hist.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
  
  return new_hist
  
def copy_hist_style(hfrom,hto):
  hto.SetLineColor  (hfrom.GetLineColor())
  hto.SetLineStyle  (hfrom.GetLineStyle())
  hto.SetLineWidth  (hfrom.GetLineWidth())
  
  hto.SetFillColor  (hfrom.GetFillColor())
  hto.SetFillStyle  (hfrom.GetFillStyle())
  
  hto.SetMarkerColor(hfrom.GetMarkerColor())
  hto.SetMarkerStyle(hfrom.GetMarkerStyle())
  hto.SetMarkerSize (hfrom.GetMarkerSize ())
  
def SetBLMargins(canv):
  canv.SetBottomMargin(0.16)
  canv.SetLeftMargin  (0.11)
  return
  
def SetRMargins(canv,needZ = False):
  if needZ: 
    canv.SetRightMargin(0.1)
  else: 
    canv.SetRightMargin(0.05)
  return

def make_frame(canvas,h,xmin,xmax,ymin=None,ymax=None):
  if ymin is None: ymin = h.GetMinimum()
  if ymax is None: ymax = h.GetMaximum()
  frame = canvas.DrawFrame(xmin,ymin,
                           xmax,ymax)
                           
  for oldaxis,newaxis in [(h.GetXaxis(),frame.GetXaxis()),
                          (h.GetYaxis(),frame.GetYaxis())]:
    newaxis.SetTitle       (oldaxis.GetTitle())
    newaxis.SetTitleSize   (oldaxis.GetTitleSize())
    newaxis.SetTitleOffset (oldaxis.GetTitleOffset())
    
    newaxis.SetLabelSize   (oldaxis.GetLabelSize())
    newaxis.SetLabelOffset (oldaxis.GetLabelOffset())
    
  return frame

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] valor mach3 ptheta",description="Make plots.")
  parser.add_option('-t', dest='tag'       ,default="data"  ,help='The filename')
  parser.add_option('-r', dest='rebin'     ,default=False,action="store_true"  ,help="Use the coarse binning")
  parser.add_option('-n', dest='nopave'    ,default=False,action="store_true"  ,help="Don't draw the paves")
  parser.add_option('-l', dest='altlabels' ,default=False,action="store_true"  ,help="Do Frequentist-vs-Bayesian rather than Group-by-Group")
  parser.add_option('-c', '--config'        , dest='config' , default="analysis_configuration.xml", help="Configuration file")
  
  options,patterns = parser.parse_args()
  
  assert(len(patterns)==3)
  
  main(patterns, options.config, options.tag, options.rebin, options.nopave,options.altlabels)
