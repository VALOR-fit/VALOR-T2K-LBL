#!/usr/bin/python

# Plots contours
# Steve Dennis, March 2015

import ROOT
import math
import array
import numpy as np
import sys
import os

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
R=ROOT

#~ R.gROOT.Reset()
#~ R.gROOT.SetStyle("T2Kcol")
#~ R.gROOT.ForceStyle(True)
  
R.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

opt.gOscMixingParamChoice = conv.kNuMixSinsqTheta
opt.SetDefaultParameterValues()

image_format = "eps"
conf90 = 4.61
conf68 = 2.30

to3dp = "{0:.3f}"
to6dp = "{0:.6f}"

p1 = conv.kPhypMixing23Bar
p2 = conv.kPhypDm23sqBar

# First file is stat-only, second is stat+syst
def main(contour_file, bf_normal =None, bf_inverted =None, tag=None,draw_kdi=True,draw_68=False):
  
  ROOT.gStyle.SetPalette(56);
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  
  title_size = 0.08
  label_size = 0.07
  
  line_style_68 = ROOT.kDashed
  line_style_90 = ROOT.kSolid
  
  kdi_file = ROOT.TFile(os.environ["T2KOSCANA"]+"/inputs/contours/KDI_run1-4_nuis4_sintheta_normal_hist.root")
  kdi_graph = kdi_file.Get("contour90_constdchisq")
  kdi_points = kdi_file.Get("bestfit")
  kdi_graph.SetLineColor(ROOT.kGreen+3)
  kdi_points.SetMarkerColor(ROOT.kGreen+3)
  kdi_graph.SetMarkerColor(ROOT.kGreen+3)
  
  
  f_nh = ROOT.TFile(contour_file)
  h_nh = f_nh.Get("h_chisq_dmsq_mixing23_Ave_"+"nh")
  assert(h_nh)
  h_nh.Print()
  min_nh = h_nh.GetMinimum()
  h_nh.SetTitle("Normal")
  
  f_ih = f_nh
    
  h_ih = f_ih.Get("h_chisq_dmsq_mixing23_Ave_"+"ih")
  assert(f_ih)
  h_ih.Print()
  min_ih = h_ih.GetMinimum()
  h_ih.SetTitle("Inverted")
  
  print min_nh
  print min_ih
  
  # Set title offset
  h_nh.GetXaxis().SetTitleOffset(0.8)
  h_ih.GetXaxis().SetTitleOffset(0.8)
  h_nh.GetYaxis().SetTitleOffset(0.9)
  h_ih.GetYaxis().SetTitleOffset(0.9)
    
  h_nh.GetXaxis().SetTitleSize(title_size)
  h_ih.GetXaxis().SetTitleSize(title_size)
  h_nh.GetYaxis().SetTitleSize(title_size)
  h_ih.GetYaxis().SetTitleSize(title_size)
  h_nh.GetXaxis().SetLabelSize(label_size)
  h_ih.GetXaxis().SetLabelSize(label_size)
  #~ h_nh.GetYaxis().SetLabelSize(label_size)
  #~ h_ih.GetYaxis().SetLabelSize(label_size)
  
  # Precalculate the best-fit points from the histogram
  bf_bin_ih = h_ih.GetMinimumBin()
  bf_bin_nh = h_nh.GetMinimumBin()
  
  bf_chisq_ih = h_ih.GetBinContent(bf_bin_ih)
  bf_chisq_nh = h_nh.GetBinContent(bf_bin_nh)
  
  x_bin_ih = ROOT.Long() ; y_bin_ih = ROOT.Long()
  x_bin_nh = ROOT.Long() ; y_bin_nh = ROOT.Long()
  dummy = ROOT.Long()
  
  h_ih.GetBinXYZ(bf_bin_ih,x_bin_ih,y_bin_ih,dummy)
  h_nh.GetBinXYZ(bf_bin_nh,x_bin_nh,y_bin_nh,dummy)
  
  bf_x_ih = h_ih.GetXaxis().GetBinCenter(x_bin_ih)
  bf_x_nh = h_nh.GetXaxis().GetBinCenter(x_bin_nh)
  bf_y_ih = h_ih.GetYaxis().GetBinCenter(y_bin_ih)
  bf_y_nh = h_nh.GetYaxis().GetBinCenter(y_bin_nh)
  
  
  # If we have files, go and get the bestfits from those instead.
  if bf_normal is not None:
    tf_normal = ROOT.TFile(bf_normal)
    tr = tf_normal.Get("fitnt")
    tr.GetEntry(0)
    bf_x_nh = tr.best_fit_params.BestFitValue(p1)
    bf_y_nh = tr.best_fit_params.BestFitValue(p2)
    bf_chisq_nh = tr.best_fit_params.GetBestFitStatistic()

  if bf_inverted is not None:
    tf_inverted = ROOT.TFile(bf_inverted)
    tr = tf_inverted.Get("fitnt")
    tr.GetEntry(0)
    bf_x_ih = tr.best_fit_params.BestFitValue(p1)
    bf_y_ih = tr.best_fit_params.BestFitValue(p2)
    bf_chisq_ih = tr.best_fit_params.GetBestFitStatistic()
      
  
  h_delta_nh = GetDeltaHist(h_nh,bf_chisq_nh)
  h_delta_nh.SetMinimum(0)
  
  h_contour_nh = ROOT.TH2D(h_delta_nh)
  h_contour_nh.SetContour(2)
  h_contour_nh.SetContourLevel(1,conf90)
  h_contour_nh.SetContourLevel(0,conf68)
  h_contour_nh.SetLineWidth(1)
  h_contour_nh.SetLineColor(R.kBlack)
  
  h_delta_ih = GetDeltaHist(h_ih,bf_chisq_ih)
  h_delta_ih.SetMinimum(0)
  
  h_contour_ih = ROOT.TH2D(h_delta_ih)
  h_contour_ih.SetContour(2)
  h_contour_ih.SetContourLevel(1,conf90)
  h_contour_ih.SetContourLevel(0,conf68)
  h_contour_ih.SetLineWidth(1)
  h_contour_ih.SetLineStyle(2)
  h_contour_ih.SetLineColor(R.kBlack)
  
  contours_ih = GetContours(h_contour_ih)
  contours_nh = GetContours(h_contour_nh)
  
  marker_ih = ROOT.TMarker(bf_x_ih,bf_y_ih,ROOT.kFullStar)
  marker_nh = ROOT.TMarker(bf_x_nh,bf_y_nh,ROOT.kFullStar)
  
  if (bf_chisq_ih > bf_chisq_nh):
    bf_is_normal = True
    marker_nh.SetMarkerStyle(ROOT.kOpenStar)
  else:
    bf_is_normal = False
    marker_ih.SetMarkerStyle(ROOT.kOpenStar)
    
  marker_ih.SetMarkerColor(ROOT.kViolet);
  marker_nh.SetMarkerColor(ROOT.kBlue+1);
  marker_ih.SetMarkerSize(1.5);
  marker_nh.SetMarkerSize(1.5);
  kdi_points.SetMarkerSize(1.5)
  
  CopyMarkerToHist(marker_ih,h_contour_ih)
  CopyMarkerToHist(marker_nh,h_contour_nh)
  CopyMarkerToHist(kdi_points,kdi_graph)
  
  
  # Dump graphs for Kirsty
  outfile = ROOT.TFile("valor_fit_contours.root","RECREATE")
  outfile.cd()
  
  h_delta_nh.Write("valor_nh_dchisq")
  h_delta_ih.Write("valor_ih_dchisq")
  
  contours_nh[0][0].Write("valor_nh_68_contour")
  contours_nh[1][0].Write("valor_nh_90_contour")
  contours_ih[0][0].Write("valor_ih_68_contour")
  contours_ih[1][0].Write("valor_ih_90_contour")
  
  for c in contours_nh[0]:
    c.SetLineColor(ROOT.kRed)
    c.SetLineStyle(h_contour_nh.GetLineStyle())
  for c in contours_nh[1]:
    c.SetLineColor(ROOT.kBlack)
    c.SetLineStyle(h_contour_nh.GetLineStyle())
  for c in contours_ih[0]:
    c.SetLineColor(ROOT.kRed)
    c.SetLineStyle(h_contour_ih.GetLineStyle())
  for c in contours_ih[1]:
    c.SetLineColor(ROOT.kBlack)
    c.SetLineStyle(h_contour_ih.GetLineStyle())
  
  bf_nh_graph = ROOT.TGraph(1)
  bf_ih_graph = ROOT.TGraph(1)
  bf_nh_graph.SetTitle("Bestfit: Normal")
  bf_ih_graph.SetTitle("Bestfit: Inverted")
  bf_nh_graph.SetPoint(0,bf_x_nh,bf_y_nh)
  bf_ih_graph.SetPoint(0,bf_x_ih,bf_y_ih)
  bf_nh_graph.Write("valor_nh_bestfit")
  bf_ih_graph.Write("valor_ih_bestfit")
  outfile.Close()
  
  print "Bestfit point (normal)  : ",bf_x_nh,bf_y_nh,"\twith chisq: ",bf_chisq_nh
  print "Bestfit point (inverted): ",bf_x_ih,bf_y_ih,"\twith chisq: ",bf_chisq_ih
  
  if bf_is_normal:
    print "\n",to3dp.format(bf_x_nh), " & ",to6dp.format(bf_y_nh)," & Normal"
  else:
    print "\n",to3dp.format(bf_x_ih), " & ",to6dp.format(bf_y_ih)," & Inverted"
    
  print "\n\n"
  
  # Start drawing
  canv = ROOT.TCanvas("canv1","canv1",1280,720)
  SetBLMargins(canv)
  SetRMargins(canv,False)
  canv.cd()
  
  # Set Axis labels
  #~ if (h_nh.GetXaxis().GetTitle()!=h_ih.GetXaxis().GetTitle()):
  new_title = h_nh.GetYaxis().GetTitle() + " or "+util.ParamAsStringPaperRooTex(conv.kPhypDm23sq)
  h_contour_nh.GetYaxis().SetTitle(new_title)
  h_contour_ih.GetYaxis().SetTitle(new_title)
  #~ if (h_nh.GetYaxis().GetTitle()!=h_ih.GetYaxis().GetTitle()):
    #~ new_title = h_nh.GetYaxis().GetTitle() + " / "+h_ih.GetYaxis().GetTitle()
    #~ h_contour_nh.GetYaxis().SetTitle(new_title)
    #~ h_contour_ih.GetYaxis().SetTitle(new_title)
  
  h_contour_nh.GetYaxis().SetRangeUser(1E-3,4.5E-3)
  h_contour_nh.Draw("AXIS")
  #~ for icont,cont_list in enumerate(contours_ih + contours_nh):
  for icont,cont_list in enumerate(contours_nh):
    if icont%2==0 and not draw_68: continue
    for cont in cont_list:
      cont.Draw("same")
      
  # Draw markers
  #~ marker_ih.Draw()
  marker_nh.Draw()
  
  # Draw MINOS
  if draw_kdi:
    kdi_graph.Draw("l same")
    kdi_points.Draw("p same")
  
  # Make legend
  
  legend = ROOT.TLegend(0.55,0.6,0.93,0.93)
  legend.AddEntry(h_contour_nh,"T2K Run 5-6","lp")
  #~ legend.AddEntry(h_contour_ih,"IH","lp")
  if draw_kdi:
    legend.AddEntry(kdi_graph,"T2K Run 1-4","lp")
  legend.SetFillStyle(0)
  legend.SetBorderSize(0)
  legend.Draw()
  
  if tag is None:
    plot_name = "profiled_nh_vs_ih"
  else:
    plot_name = tag
  
  canv.Print(plot_name+"."+image_format)
  
  canv_nh = ROOT.TCanvas("canv_nh","canv_nh",1280,720)
  SetBLMargins(canv_nh)
  SetRMargins(canv_nh,True)
  canv_nh.cd()
  h_delta_nh.SetContour(100)
  h_delta_nh.Draw("colz")
  marker_nh.Draw()
  if draw_68:
    for g in contours_nh[0]:
      g.SetLineStyle(2)
      g.SetLineColor(R.kBlue)
      g.Draw("same")
    
  for g in contours_nh[1]:
    g.SetLineStyle(1)
    g.SetLineColor(R.kBlue)
    g.Draw("same")
  
  if draw_kdi:  
   kdi_graph.Draw("l same")
   kdi_points.Draw("p same")
   
  canv_nh.Print(plot_name+"_nh."+image_format)
  
  #~ canv_ih = ROOT.TCanvas("canv_ih","canv_ih",1280,720)
  #~ canv_ih.cd()
  #~ SetBLMargins(canv_ih)
  #~ SetRMargins(canv_ih,True)
  #~ h_delta_ih.SetContour(100)
  #~ h_delta_ih.Draw("colz")
  #~ marker_ih.Draw()
  #~ 
  #~ if draw_68:
    #~ for g in contours_ih[0]:
      #~ g.SetLineStyle(2)
      #~ g.SetLineColor(R.kBlue)
      #~ g.Draw("same")
      #~ 
  #~ for g in contours_ih[1]:
    #~ g.SetLineStyle(1)
    #~ g.SetLineColor(R.kBlue)
    #~ g.Draw("same")
  #~ if draw_kdi:  
    #~ kdi_graph.Draw("l same")
    #~ kdi_points.Draw("p same")
  #~ canv_ih.Print(plot_name+"_ih."+image_format)

def GetContours(h):
  dummy_canv = ROOT.TCanvas()
  dummy_canv.cd()
  h.Draw("cont list")
  dummy_canv.Update()
  obj_array = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
  del dummy_canv
  tlists         = [obj_array[i]         for i in xrange(obj_array.GetEntries())]
  contour_levels = [h.GetContourLevel(i) for i in xrange(obj_array.GetEntries())]
  to_return = []
  for tl,level in zip(tlists,contour_levels):
    cont_list = [ROOT.TGraph(tl[i]) for i in xrange(tl.GetEntries())]
    for c in cont_list:
      c.SetName ("contour_"+to3dp.format(level))
      c.SetTitle(str(h.GetTitle())+" Contour "+to3dp.format(level))
    to_return.append(cont_list)
  return to_return
  
def CopyMarkerToHist(marker,histogram):
  histogram.SetMarkerSize (marker.GetMarkerSize ())
  histogram.SetMarkerStyle(marker.GetMarkerStyle())
  histogram.SetMarkerColor(marker.GetMarkerColor())
  
def GetDeltaHist(h,minval=None):
  h_delta = ROOT.TH2D(h)
  h_delta.Reset()
  if minval is None: minval = h.GetMinimum()
  for ix in xrange(h_delta.GetXaxis().GetNbins()):
    for iy in xrange(h_delta.GetYaxis().GetNbins()):
      h_delta.SetBinContent(ix+1,iy+1,h.GetBinContent(ix+1,iy+1)-minval)
      
  return h_delta

def SetBLMargins(canv):
  canv.SetBottomMargin(0.16)
  canv.SetLeftMargin  (0.14)
  return
  
def SetRMargins(canv,needZ = False):
  if needZ: 
    canv.SetRightMargin(0.1)
  else: 
    canv.SetRightMargin(0.05)
  return
  
if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make contours")
  
  parser.add_option('--bf-nh' , dest='bf_normal' ,  default=None, help="File containing true best-fit (syst, NH)")
  parser.add_option('--bf-ih' , dest='bf_inverted', default=None, help="File containing true best-fit (syst, IH)")
  parser.add_option('-t','--tag' , dest='tag', default=None, help="Tag for output")
  parser.add_option('-k','--kdi' , dest='draw_kdi', default=False,action="store_true", help="Should we add a MINOS contour")
  parser.add_option('--no-68' , dest='no68', default=False,action="store_true", help="Should we remove the 68% confidence region")
  
  o,args = parser.parse_args()
  assert(len(args)==1)
  
  main(args[0], o.bf_normal, o.bf_inverted, o.tag,o.draw_kdi, not o.no68)

