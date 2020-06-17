#!/usr/bin/env python

# Script to make plots of relative error envelope sizes
# Should be able to compare (roughly) N distributions, although you might need to specify more colors etc
# Takes output of valor_dune-tweak_pdf_N

#Steve Dennis, November 2016

import ROOT
import math
import os

ROOT.PyConfig.IgnoreCommandLineOptions = True

R = ROOT

max_energy = 10.
linewidth = 2
title_offset = 1.1
left_margin   = 0.18
bottom_margin = 0.20
right_margin  = 0.05
use_log_plots = False
title_size = 0.08
label_size = 0.07
font = 132 # Times Medium, scalable rotatable
canvas_x = 1280
canvas_y = 720

fill_colors = [ROOT.kGray+1, ROOT.kRed, ROOT.kGreen+2, R.kBlue]
mark_colors = [ROOT.kGray+1, ROOT.kRed, ROOT.kGreen+2, R.kBlue]

#~ prefit_style = 3244 # hatched
prefit_style  = 1001 # Solid
postfit_style = 1001 # Solid

fill_styles = [1001] * 4
line_styles = [R.kDashed,R.kSolid,R.kDotted]

to3dp = "{0:.3f}"

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

#~ colors = [ROOT.kBlack, R.kRed, R.kBlue, ROOT.kMagenta+1,  ROOT.kGreen+2, ROOT.kOrange+7, ROOT.kCyan+2,]

plot_opt = ""

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C");
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

def main(config, filenames, labels, tag, write_tex, img_format,ratio_mult):
  files = []
  for f in filenames:
    tf = ROOT.TFile(f)
    if not tf:
      raise DUNEAnaError("Could not open input file: "+f)
    files.append(tf)
  
  print files
  
  envelopes(config, files, labels, tag, write_tex, img_format,ratio_mult)
  

def envelopes(config, files, labels, tag, write_tex, img_format, ratio_mult):
  
  all_hists = []
    
  envelope_images = []
  ratio_images = []
  envelope_captions = []
  ratio_captions = []
  
  elem = valor.CompoundFitElements()
  elem.Read(util.CompleteConfigPath(config),"nomc nosplines")
  cm = elem.FitModel()
  
  icanv = 0
  
  for imodel in xrange(cm.GetNSubModels()):
    icanv+=1
    
    m = cm.GetSubModel(imodel)
    s = m.GetSample()
    b = m.GetBeam()
    samp_hists = []
    
    hname_prefix = "env_sb_%s_%s" % (util.BeamConfAsString(b) , util.SampleAsString(s))
    hname_nom = hname_prefix + "_nom"
    hname_p1  = hname_prefix + "_p1"
    hname_m1  = hname_prefix + "_m1"
    
    out_prefix = "%s_sb_%s_%s" % (tag, util.BeamConfAsString(b) , util.SampleAsString(s))
    
    nom_hists   = []
    plus_hists  = []
    minus_hists = []
    hists = []
    
    for ihist,f in enumerate(files):
      nom = get_hist(f,hname_nom)
      nom.SetLineColor(R.kBlack)
      nom.SetLineStyle(1)
      nom.SetLineWidth(linewidth)
      nom_hists  .append(nom)
      plus_hists .append(get_hist(f,hname_p1))
      minus_hists.append(get_hist(f,hname_m1))
    
    hists = nom_hists+plus_hists+minus_hists
    
    for h in hists:
      samp_hists.append(h)
      all_hists.append(h)
      h.UseCurrentStyle()
      
      h.GetXaxis().SetTitleSize(title_size)
      h.GetYaxis().SetTitleSize(title_size)
      h.GetXaxis().SetLabelSize(label_size)
      h.GetYaxis().SetLabelSize(label_size)
      h.GetXaxis().SetTitle("E_{reco}")
      h.GetYaxis().SetTitle("Events per bin")
      h.SetTitle("")
      h.GetXaxis().SetRangeUser(0.,max_energy)
      h.GetYaxis().SetTitleOffset(title_offset)
    
    for h in plus_hists+minus_hists:
      h.SetLineStyle(0)
      h.SetLineWidth(0)
      h.GetXaxis().SetRangeUser(0.,max_energy)
      
    legend = ROOT.TLegend(0.7,0.7,0.93,0.92);
    legend.SetFillStyle(0);
    legend.SetBorderSize(0);
    #~ legend.SetTextSize(0.04);
   
    limit = max([abs(h.GetMinimum()) for h in samp_hists] + [abs(h.GetMaximum()) for h in samp_hists])
    
    # Hack so that the legend stays above the line
    val_at_legend_start = max([abs(h.GetBinContent(h.GetXaxis().FindBin(max_energy*0.7))) for h in samp_hists])
    limit = max(limit,val_at_legend_start / 0.65)
    
    # Start plotting
    cname = "c"+`icanv`+"_env"
    canv = ROOT.TCanvas(cname,cname,canvas_x,canvas_y)
    canv.SetLeftMargin  (left_margin)
    canv.SetBottomMargin(bottom_margin)
    canv.SetRightMargin (right_margin )
    if use_log_plots: canv.SetLogx()
    
    for h in samp_hists:
      h.SetMinimum(0.)
      h.SetMaximum(1.05*limit)
      
    env_hists = zip(nom_hists,plus_hists,minus_hists)
    env_graphs = [ROOT.TGraphAsymmErrors(h.GetXaxis().GetNbins()) for h in plus_hists]
    
    for (hn,hp,hm),g in zip(env_hists,env_graphs):
      for i in xrange(hp.GetXaxis().GetNbins()):
        low_edge  = hp.GetBinLowEdge(i+1)
        high_edge = low_edge + hp.GetBinWidth(i+1)
        centre = hn.GetBinContent(i+1)
        g.SetPoint      (i,hp.GetBinCenter(i+1),centre)
        g.SetPointEYlow (i,abs(hm.GetBinContent(i+1)-centre))
        g.SetPointEYhigh(i,abs(hp.GetBinContent(i+1)-centre))
        g.SetPointEXlow (i,hp.GetBinWidth(i+1)*0.5)
        g.SetPointEXhigh(i,hp.GetBinWidth(i+1)*0.5)
        
      g.GetXaxis().SetRangeUser(0,max_energy)
      CopyHistAxisFormatting(hp,g)
      g.GetYaxis().SetTitle("Ratio of # Total Events")
    
    plus_hists[0].Draw("AXIS")
    
    for ig,env_graph in enumerate(env_graphs):
      env_graph.SetFillStyle(fill_styles[ig])
      env_graph.SetFillColor(fill_colors[ig])
      env_graph.SetMarkerColor(mark_colors[ig])
      env_graph.SetMarkerStyle(ROOT.kFullSquare)
      env_graph.SetLineWidth(0)
      env_graph.SetLineStyle(0)
      env_graph.SetLineWidth(0)
      env_graph.Draw("same2")
    
    for nom in nom_hists:
      nom.Draw("same")
    
    canv.RedrawAxis()
    
    legend_hists = []
    for nom,graph in zip(nom_hists,env_graphs):
      legend_hist = ROOT.TH1D(nom)
      legend_hist.SetLineStyle(nom.GetLineStyle())
      legend_hist.SetLineColor(nom.GetLineColor())
      legend_hist.SetLineWidth(nom.GetLineWidth())
      legend_hist.SetFillStyle(graph.GetFillStyle())
      legend_hist.SetFillColor(graph.GetFillColor())
      legend_hists.append(legend_hist)
    
    all_hists.extend([legend_hists])
    
    for label,hist in zip(labels,legend_hists):
      legend.AddEntry(hist  ,label,"FL")
    
    legend.Draw()
    envelope_fn = out_prefix+"_envelope."+img_format
    canv.Print(envelope_fn)
    envelope_images.append(envelope_fn)
    cap = "%s %s %s"%(tag,util.BeamConfAsString(b),util.SampleAsString(s))
    envelope_captions.append(cap.replace("_"," "))

#______________SCRIPT NOT UPDATED BEYOND THIS POINT___________________________________________

    continue

#______________SCRIPT NOT UPDATED BEYOND THIS POINT___________________________________________

    # Now do it with ratios
    ratios_plus   = [R.TH1D(hp) for hp in plus_hists]
    ratios_minus  = [R.TH1D(hm) for hm in minus_hists]
    
    central_hist = nom_hists[0] # Should this be the specific one for this histogram, or the first (generally prefit) one
    ratios = ratios_plus + ratios_minus
    for ratio in ratios:
      ratios.Divide(central_hist)
      
    for i in xrange(central_hist.GetXaxis().GetNbins()):
      if abs(central_hist.GetBinContent(i+1))<conv.kASmallNum:
        for ratio in ratios:
          ratio.SetBinContent(i+1,1)
    
    ratio_graph_prefit  = ROOT.TGraphAsymmErrors(ratio_plus_prefit  .GetXaxis().GetNbins())
    ratio_graph_postfit = ROOT.TGraphAsymmErrors(ratio_plus_postfit .GetXaxis().GetNbins())
    
    ratio_graphs = [ratio_graph_prefit  ,
                    ratio_graph_postfit , ]
    
    if ratio_mult is not None:
      eff_mult = ratio_mult
    else:
      eff_mult = 1.
      
    fill_ratio_graph_from_hist(ratio_graph_prefit , ratio_minus_prefit , ratio_plus_prefit , 1.)
    fill_ratio_graph_from_hist(ratio_graph_postfit, ratio_minus_postfit, ratio_plus_postfit, eff_mult)
      
    ratio_graph_prefit .SetFillColor(prefit_color)
    ratio_graph_postfit.SetFillColor(postfit_color)
    ratio_graph_prefit .SetLineColor(prefit_color)
    ratio_graph_postfit.SetLineColor(postfit_color)
    
    ratio_graph_prefit .SetLineWidth(3)
    ratio_graph_postfit.SetLineWidth(3)
    ratio_graph_prefit .SetFillStyle(1001)
    ratio_graph_postfit.SetFillStyle(1001)
    
    cname2 = "c"+`icanv`+"_ratio"
    canv2 = ROOT.TCanvas(cname2,cname2,canvas_x,canvas_y)
    canv2.SetLeftMargin  (left_margin)
    canv2.SetBottomMargin(bottom_margin)
    canv2.SetRightMargin (right_margin )
    if use_log_plots: canv2.SetLogx()
    canv2.cd()
    
    ratio_max = 1.3*ratio_plus_prefit .GetMaximum()
    
    #~ ratio_graph_postfit.SetMaximum(ratio_max)
    #~ ratio_graph_postfit.SetMinimum(2-ratio_max)
    ratio_graph_prefit.GetYaxis().SetRangeUser(0.,2.5)
    
    ratio_graph_prefit .Draw("a2")
    ratio_graph_postfit.Draw("2")
    canv2.RedrawAxis()
    
    legend2 = ROOT.TLegend(0.5,0.85,0.93,0.94);
    legend2.SetFillStyle(0);
    legend2.SetBorderSize(0);
    legend2.SetNColumns(2)
    legend2.AddEntry(ratio_graph_prefit  ,prefit_label ,"F")
    if ratio_mult is not None:
      legend2.AddEntry(ratio_graph_postfit ,"%s#times%g"%(postfit_label,eff_mult),"FL")
    else:
      legend2.AddEntry(ratio_graph_postfit ,postfit_label,"FL")
    
    
    legend2.Draw("same")
    canv2.Update() # Need to update here so line gets right maximum
    
    line_at_1 = ROOT.TLine(0,1,canv2.GetUxmax(),1)
    line_at_1.SetLineWidth(2)
    line_at_1.Draw("same")
    
    canv2.Update()
    ratio_fn = out_prefix+"_envelope_ratio."+img_format
    canv2.Print(ratio_fn)
    ratio_images.append(ratio_fn)
    cap = "%s %s %s"%(tag,util.BeamConfAsString(b),util.SampleAsString(s))
    ratio_captions.append(cap.replace("_"," "))
    
    #~ # Print out some event rates
    #~ nom_events     = prefit_nom  .Integral()
    #~ p1_events_pre  = hplus_prefit  .Integral()
    #~ m1_events_pre  = hminus_prefit .Integral()
    #~ p1_events_post = hplus_postfit .Integral()
    #~ m1_events_post = hminus_postfit.Integral()
    #~ 
    #~ p1_percent_pre  = 100.*( (p1_events_pre -nom_events)/nom_events )
    #~ m1_percent_pre  = 100.*( (m1_events_pre -nom_events)/nom_events )
    #~ p1_percent_post = 100.*( (p1_events_post-nom_events)/nom_events )
    #~ m1_percent_post = 100.*( (m1_events_post-nom_events)/nom_events )
    #~ 
    #~ items = [name+"_"+par_set,to3dp.format(nom_events)]
    #~ items.append(to3dp.format(p1_events_pre ) +" ("+to3dp.format(p1_percent_pre )+"%)")
    #~ items.append(to3dp.format(m1_events_pre ) +" ("+to3dp.format(m1_percent_pre )+"%)")
    #~ items.append(to3dp.format(p1_events_post) +" ("+to3dp.format(p1_percent_post)+"%)")
    #~ items.append(to3dp.format(m1_events_post) +" ("+to3dp.format(m1_percent_post)+"%)")
    #~ print "\t".join(items)
    
  # End of model loop
  
  if (write_tex):
    tex = start_latex_file("%s_envelopes.tex"%(tag))
    add_figures_to_tex(tex,envelope_images,envelope_captions)
    tex.write("\n\n\\clearpage\n\n")
    add_figures_to_tex(tex,ratio_images,ratio_captions)
    end_latex_file(tex)
    
  return

class DUNEAnaError(Exception):
  pass
  
def CopyHistAxisFormatting(hfrom,hto):
  for (afrom,ato) in [(hfrom.GetXaxis(),hto.GetXaxis()),
                      (hfrom.GetYaxis(),hto.GetYaxis())]:
    ato.SetTitle      (afrom.GetTitle      ())
    ato.SetTitleOffset(afrom.GetTitleOffset())
    ato.SetTitleSize  (afrom.GetTitleSize  ())
    ato.SetLabelSize  (afrom.GetLabelSize  ())
    ato.SetLabelOffset(afrom.GetLabelOffset())
    
  return
  
def get_hist(rootfile,hname):
  h = rootfile.Get(hname)
  if not h:
    print "Failed to get histogram: ",hname ; exit()
  return ROOT.TH1D( h )
  
def start_latex_file(name):
  f = open(name,"w")
  f.write(r"""\documentclass[10pt]{article}
\usepackage[utf8]{inputenc}
\usepackage{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{float}
\usepackage{setspace}

\setlength\parindent{0pt}

\begin{document}

\begin{doublespacing}

""")
  return f
  
def end_latex_file(f):
  f.write(r"""
\end{doublespacing}
\end{document}""")
  f.close()

def add_figures_to_tex(f,figure_names,captions,plots_per_row=2):
  figure_width  = str(0.99/float(plots_per_row))+"\\textwidth"
  
  n_figs = len(figure_names)
  for i in xrange(n_figs):
    if (i%plots_per_row == 0): 
      f.write(" \\begin{figure*}\n")
    image = figure_names[i]
    caption = captions[i]
    label =  "fig:"+image.split("/")[-1].rsplit(".",1)[0]
    f.write("\n \\begin{minipage}[t]{"+figure_width+"}\n")
    f.write("  \\includegraphics[width=0.95\\textwidth]{"+image+"}\n")
    f.write("  \\caption{"+caption+"}\n")
    f.write("  \\label{"+label+"}\n")
    f.write("\\end{minipage}")
    if i%plots_per_row or i==n_figs-1:
      f.write("\\end{figure*}\n")
  return
  
def fill_ratio_graph_from_hist(g,hm,hp,multiplier=1.):
  for i in xrange(hp.GetXaxis().GetNbins()):
    low_edge  = hm.GetBinLowEdge(i+1)
    high_edge = low_edge + hm.GetBinWidth(i+1)
    g.SetPoint      (i,hm.GetBinCenter(i+1),1)
    g.SetPointEYlow (i,multiplier*abs(hm.GetBinContent(i+1)-1))
    g.SetPointEYhigh(i,multiplier*abs(hp.GetBinContent(i+1)-1))
    g.SetPointEXlow (i,hm.GetBinWidth(i+1)*0.5 - 0.005)
    g.SetPointEXhigh(i,hm.GetBinWidth(i+1)*0.5 - 0.005)
          
  g.SetFillStyle(1001)
  g.SetLineStyle(0)
  g.SetLineWidth(0)
  g.GetXaxis().SetRangeUser(0,max_energy)
  CopyHistAxisFormatting(hp,g)
  g.GetYaxis().SetTitle("Ratio of # Total Events")
  return

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] files",description="Make envelopes of error effects on spectrum")
  parser.add_option('-c','--config', dest='config'  ,default="fit/dunend_fit.xml",help='The config file we get the sample list from.')
  parser.add_option('-l','--label', dest='labels'  ,default="Prefit,Postfit"  ,help='The label for the envelopes')
  parser.add_option('-f','--format', dest='imgformat',default='eps',help='The image format')
  parser.add_option('-t','--tag', dest='tag',default='env',help='Prefix for output filenames')
  parser.add_option('-w','--write-tex', dest='tex',action="store_true",default=False,help='Write TeX file with images in.')
  parser.add_option('-r','--ratio-mult', dest='ratiomult',type=float,default=None,help='Multiply the postfit envelope size in the ratio plots by this number.')
  
  options,files = parser.parse_args()
  
  labels = options.labels.split(",")
  assert(len(labels)==len(files))
  
  main(options.config, files, labels ,options.tag,options.tex,options.imgformat,options.ratiomult)
