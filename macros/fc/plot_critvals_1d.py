#!/usr/bin/env python

#
# For each FC point (for fixed FC hierarchy):
#   Collect the FC Delta chisq value for each FC toy.
#   Use this calculate critical values
# Assumes that this is for delta CP with 9 FC points, but this should be fairly easy to make more generic, if needed
#
# Andy Chappell
# Francis Bench <F.Bench \at liverpool.ac.uk>
# August 2019

import ROOT
import numpy as np
import sys
import os

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)

try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")

valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

_conf_levels = [conv.kCL1Sigma, conv.kCL90, conv.kCL2Sigma, conv.kCL3Sigma]
_colours     = [ROOT.kOrange-4  , ROOT.kPink-5 , ROOT.kAzure+8    , ROOT.kGreen+4]
_line_styles = [1, 2]
_n_points = 11
_ymax = 10

def plot_critical_values(output_file, input_file_nh, input_file_ih, legacy, prelim):
  canvas = ROOT.TCanvas("canvas","canvas")
  canvas.DrawFrame(-3.6, 0.0, 3.6, 5.0)
  legend = ROOT.TLegend(0.92,0.60,1.0,0.95)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetFillColorAlpha(0,0)
  legend.SetTextSize(0.04)
  g_shade = dict()

  for n_idx, cl in enumerate(_conf_levels):
    for m_idx, (input_file, hier_name) in enumerate(zip([input_file_nh, input_file_ih],["normal","inverted"])):
      if not input_file:
        continue
      f = ROOT.TFile(input_file, "READ")
      graph    = f.Get(get_graph_name(cl, None, hier_name, legacy))
      graph_p1 = f.Get(get_graph_name(cl, "p1", hier_name, legacy))
      graph_m1 = f.Get(get_graph_name(cl, "m1", hier_name, legacy))

      if not graph:
        print "ERROR: Could not get graph named", get_graph_name(cl, None, hier_name, legacy)
        sys.exit()
      if not graph_p1:
        print "ERROR: Could not get graph named", get_graph_name(cl, "p1", hier_name, legacy)
        sys.exit()
      if not graph_m1:
        print "ERROR: Could not get graph named", get_graph_name(cl, "m1", hier_name, legacy)
        sys.exit()

      xx = np.zeros(_n_points)
      yy = np.zeros(_n_points)
      yy_min = np.zeros(_n_points)
      yy_max = np.zeros(_n_points)
      for i in range(0, _n_points):
        x = ROOT.Double(0)
        y = ROOT.Double(0)
        graph.GetPoint(i, x, y)
        xx[i] = x
        yy[i] = y
        graph_p1.GetPoint(i, x, y)
        yy_max[i] = y
        graph_m1.GetPoint(i, x, y)
        yy_min[i] = y
      full_name = util.ConfLevelAsString(cl) + hier_name
      g_shade[full_name] = ROOT.TGraph(2 * _n_points)
      for i in range(0, _n_points):
        g_shade[full_name].SetPoint(i, xx[i], yy_min[i])
        g_shade[full_name].SetPoint(_n_points + i, xx[_n_points - i - 1], yy_max[_n_points - i - 1])

      g_shade[full_name].SetFillColorAlpha(_colours[n_idx], 0.5)
      g_shade[full_name].SetFillStyle(1000)
      g_shade[full_name].SetTitle("")
      g_shade[full_name].GetXaxis().SetTitle("#delta_{CP} (Radians)")
      g_shade[full_name].GetYaxis().SetTitle("#Delta #chi^{2}_{crit}")
      if n_idx == 0 and m_idx == 0:
        g_shade[full_name].GetYaxis().SetRangeUser(0, _ymax)
        g_shade[full_name].GetYaxis().SetTitleOffset(0.8)
        #g_shade[full_name].GetYaxis().CenterTitle()
        #g_shade[full_name].GetXaxis().CenterTitle()
        g_shade[full_name].Draw("af")
      else:
        g_shade[full_name].Draw("f same")
      graph.SetLineColor(_colours[n_idx])
      graph.SetLineStyle(_line_styles[m_idx])
      graph.Draw("lp same")
      if m_idx == 0:
        legend.AddEntry(graph, util.ConfLevelAsStringRooTex(cl), "l")
      if n_idx == len(_conf_levels)-1 and m_idx == 0:
        legend.AddEntry(graph, "NO", "l")
      elif n_idx == len(_conf_levels)-1 and m_idx == 1:
        legend.AddEntry(graph, "IO", "l")

      f.Close()

  legend.Draw()

  title = ROOT.TLatex(0.0, 1.01 * 5.0, "T2K Run1-9 Preliminary")
  title.SetTextSize(0.05)
  if prelim:
    title.Draw()

  canvas.SaveAs(output_file + ".pdf")
  print "Done!"

#_______________________________________________________________________
def get_graph_name(conf_level_enum, suffix, hier_name, legacy):
  name = "gFC"
  if legacy:
    name += util.ConfLevelAsString(conf_level_enum).replace("pct","CL")
    if hier_name == "normal":
      name += "_nh"
    else:
      name += "_ih"
  else:
    name += "_" + util.ConfLevelAsString(conf_level_enum)
  if suffix:
    name += "_" + suffix
  return name

#_______________________________________________________________________
if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Plot FC Delta chisq values")

  # Input/Output
  parser.add_option('-o', '--out-tag'   , dest='output_tag', default="", help="Tag used to name the output file")
  parser.add_option('-n', '--file-in-nh', dest='file_in_nh', default="", help="Name of the input NH file containing 1D TGraphs of critical values")
  parser.add_option('-i', '--file-in-ih', dest='file_in_ih', default="", help="Name of the input IH file containing 1D TGraphs of critical values")
  parser.add_option('-l', '--legacy'    , dest='legacy'    , action="store_true", default=False, help="For reading FC files produced by the legacy T2K analysis.")
  parser.add_option('-p', '--prelim'    , dest='prelim'    , action="store_true", default=False, help="Draws a text box indicating that the plot is preliminary.")

  o,other = parser.parse_args()

  if not (o.file_in_nh or o.file_in_ih):
    print "ERROR: You need to input an NH or IH file, or both"
    sys.exit()

  plot_critical_values(o.output_tag, o.file_in_nh, o.file_in_ih, o.legacy, o.prelim)

  if other:
    print "Warning - ignored arguments: ", other
