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

# For comparison with critical values produced by the p-theta group:
_graph_name_pattern_ptheta = "gCritical_[MO]_[CL]_[ERR]"
_ptheta_col = 16
_cl_valor_to_ptheta_index = {
  conv.kCL68     : 0,
  conv.kCL90     : 1,
  conv.kCL95     : "NONEXISTENT",
  conv.kCL99     : "NONEXISTENT",
  conv.kCL1Sigma : 0,
  conv.kCL2Sigma : 2,
  conv.kCL3Sigma : 3
}
_mo_valor_to_ptheta = {
  "normal"   : "NH",
  "inverted" : "IH"
}

def plot_critical_values(output_file, input_file_nh, input_file_ih, input_file_ptheta, legacy, prelim):
  canvas = ROOT.TCanvas("canvas","canvas")
  canvas.DrawFrame(-3.6, 0.0, 3.6, 5.0)
  legend = ROOT.TLegend(0.92,0.60,1.0,0.95)
  legend.SetBorderSize(0)
  legend.SetFillStyle(0)
  legend.SetFillColorAlpha(0,0)
  legend.SetTextSize(0.04)
  g_shade        = dict()
  g_shade_ptheta = dict()

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

      make_error_shade_graph(g_shade, graph, graph_p1, graph_m1, cl, hier_name, n_idx, m_idx, False)
      graph.SetLineColor(_colours[n_idx])
      graph.SetLineStyle(_line_styles[m_idx])
      graph.SetMarkerSize(1)
      graph.Draw("lp same")
      if m_idx == 0:
        legend.AddEntry(graph, util.ConfLevelAsStringRooTex(cl), "l")

      # If requested, plot a comparison to p-theta's critical values
      if input_file_ptheta != "":
        file_ptheta = ROOT.TFile(input_file_ptheta, "READ")
        ptheta_cl_index = _cl_valor_to_ptheta_index[cl]
        if ptheta_cl_index == "NONEXISTENT":
          print "WARNING: Could not find p-theta graphs for CL \'{}\'".foramt(util.ConfLevelAsStringRooTex(cl))
        else:
          graph_name_ptheta    = get_ptheta_graph_name(ptheta_cl_index, hier_name, "Central")
          graph_name_ptheta_p1 = get_ptheta_graph_name(ptheta_cl_index, hier_name, "Upper")
          graph_name_ptheta_m1 = get_ptheta_graph_name(ptheta_cl_index, hier_name, "Lower")
          graph_ptheta    = file_ptheta.Get(graph_name_ptheta)
          graph_ptheta_p1 = file_ptheta.Get(graph_name_ptheta_p1)
          graph_ptheta_m1 = file_ptheta.Get(graph_name_ptheta_m1)
          if not graph_ptheta or type(graph_ptheta) == ROOT.TObject:
            print "WARNING: Could not find p-theta comparison TGraph named \'{}\', so will skip this one.".format(graph_name_ptheta)
          else:
            make_error_shade_graph(g_shade_ptheta, graph_ptheta, graph_ptheta_p1, graph_ptheta_m1, cl, hier_name, n_idx, m_idx, True)
            graph_ptheta.SetLineColor(_ptheta_col)
            graph_ptheta.SetLineStyle(_line_styles[m_idx])
            graph_ptheta.Draw("l same")
            if n_idx == len(_conf_levels)-1 and m_idx == 0:
              legend.AddEntry(graph_ptheta, "p-#theta", "l")

      # Add the MO legend entries right at the end
      if n_idx == len(_conf_levels)-1 and m_idx == 0:
        legend.AddEntry(graph, "NO", "l")
      elif n_idx == len(_conf_levels)-1 and m_idx == 1:
        legend.AddEntry(graph, "IO", "l")

      f.Close()

  legend.Draw()

  title = ROOT.TLatex(0.50, 0.96, "T2K Run1-10 Preliminary")
  title.SetNDC(True)
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
def make_error_shade_graph(g_shade, graph, graph_p1, graph_m1, cl, hier_name, n_idx, m_idx, is_comparison):
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

  if not is_comparison:
    g_shade[full_name].SetFillColorAlpha(_colours[n_idx], 0.5)
  else:
    g_shade[full_name].SetFillColorAlpha(_ptheta_col, 0.35)
  g_shade[full_name].SetFillStyle(1000)
  g_shade[full_name].SetTitle("")
  g_shade[full_name].GetXaxis().SetTitle("#delta_{CP} (Radians)")
  g_shade[full_name].GetYaxis().SetTitle("#Delta #chi^{2}_{crit}")
  if n_idx == 0 and m_idx == 0 and not is_comparison:
    g_shade[full_name].GetYaxis().SetRangeUser(0, _ymax)
    g_shade[full_name].GetYaxis().SetTitleOffset(0.8)
    #g_shade[full_name].GetYaxis().CenterTitle()
    #g_shade[full_name].GetXaxis().CenterTitle()
    g_shade[full_name].Draw("af")
  else:
    g_shade[full_name].Draw("f same")

#_______________________________________________________________________
def get_ptheta_graph_name(ptheta_cl_index, hier_name, err_name):
  graph_name_ptheta = _graph_name_pattern_ptheta.replace("[CL]" , str(ptheta_cl_index))
  graph_name_ptheta = graph_name_ptheta         .replace("[MO]" , _mo_valor_to_ptheta[hier_name])
  graph_name_ptheta = graph_name_ptheta         .replace("[ERR]", err_name)
  return graph_name_ptheta

#_______________________________________________________________________
if __name__ == "__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Plot FC Delta chisq values")

  # Input/Output
  parser.add_option('-o', '--out-tag'       , dest='output_tag'    , default="", help="Tag used to name the output file")
  parser.add_option('-n', '--file-in-nh'    , dest='file_in_nh'    , default="", help="Name of the input NH file containing 1D TGraphs of critical values")
  parser.add_option('-i', '--file-in-ih'    , dest='file_in_ih'    , default="", help="Name of the input IH file containing 1D TGraphs of critical values")
  parser.add_option(      '--file-in-ptheta', dest='file_in_ptheta', default="", help="Name of the input file for comparing to the p-theta group")
  parser.add_option('-l', '--legacy'        , dest='legacy'        , action="store_true", default=False, help="For reading FC files produced by the legacy T2K analysis.")
  parser.add_option('-p', '--prelim'        , dest='prelim'        , action="store_true", default=False, help="Draws a text box indicating that the plot is preliminary.")

  o,other = parser.parse_args()

  if not (o.file_in_nh or o.file_in_ih):
    print "ERROR: You need to input an NH or IH file, or both"
    sys.exit()

  plot_critical_values(o.output_tag, o.file_in_nh, o.file_in_ih, o.file_in_ptheta, o.legacy, o.prelim)

  if other:
    print "Warning - ignored arguments: ", other
