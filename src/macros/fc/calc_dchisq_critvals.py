#!/usr/bin/env python

#
# For each FC point (for fixed FC hierarchy):
#   Collect the FC Delta chisq value for each FC toy.
#   Use this calculate critical values
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# August 2019

import ROOT
import sys
import os

sys.path.append(os.path.abspath(os.environ["VALOR_FMWK"] + "/src/valor/Utilities/"))

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)

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

_n_hier = 2
_Dchisq_list_name = "Dchi2_toy_vec"

_hier_to_index = dict()
_hier_to_index["normal"]   = 0
_hier_to_index["inverted"] = 1

_conf_levels = [conv.kCL1Sigma, conv.kCL90, conv.kCL2Sigma]


#_______________________________________________________________________

def main(output_tag, in_tag,
         npt_fc_par1, start_pt_fc_par1, min_val_fc_par1, max_val_fc_par1,
         npt_fc_par2, start_pt_fc_par2, min_val_fc_par2, max_val_fc_par2,
         fc_2d):

  # Set up TGraphs to hold the critical values
  graphs    = dict()
  graphs_p1 = dict()
  graphs_m1 = dict()
  for conf_level_enum in _conf_levels:
    if not fc_2d:
      graphs   [conf_level_enum] = ROOT.TGraph()
      graphs_p1[conf_level_enum] = ROOT.TGraph()
      graphs_m1[conf_level_enum] = ROOT.TGraph()
    else:
      graphs   [conf_level_enum] = ROOT.TGraph2D()
      graphs_p1[conf_level_enum] = ROOT.TGraph2D()
      graphs_m1[conf_level_enum] = ROOT.TGraph2D()
    graphs   [conf_level_enum].SetName (get_graph_name(conf_level_enum, None))
    graphs_p1[conf_level_enum].SetName (get_graph_name(conf_level_enum, "p1"))
    graphs_m1[conf_level_enum].SetName (get_graph_name(conf_level_enum, "m1"))
    graphs   [conf_level_enum].SetTitle(get_graph_name(conf_level_enum, None))
    graphs_p1[conf_level_enum].SetTitle(get_graph_name(conf_level_enum, "p1"))
    graphs_m1[conf_level_enum].SetTitle(get_graph_name(conf_level_enum, "m1"))

  if not fc_2d:
    if not (npt_fc_par2 == 1 and start_pt_fc_par2 == 0):
      print "WARNING: The option to run 2D FC was not supplied, but the number of points for par 2 was != 1 or " \
            "the start point for par 2was != 0. Setting to defaults."
      npt_fc_par2        = 1
      start_pt_fc_par2   = 0

  # Loop over FC grid points and calculate critical values
  for ifc_p2 in range(start_pt_fc_par2, start_pt_fc_par2 + npt_fc_par2):
    for ifc_p1 in range(start_pt_fc_par1, start_pt_fc_par1 + npt_fc_par1):

      input_filename = get_filename(in_tag, fc_2d, ifc_p1, ifc_p2)
      input_file = ROOT.TFile(input_filename, "read")
      print "Reading from file", input_filename

      dchisq_vals = list(input_file.Get(_Dchisq_list_name))
      if not dchisq_vals:
        print "ERROR: Could not read chisq list from file", input_filename
        sys.exit()
      if len(dchisq_vals) == 0:
        print "ERROR: The dchisq vector was empty"
        sys.exit()
      ntoys = len(dchisq_vals)

      # Calculate the grid point value(s) and print
      fc_grid_p1_val = fc_pt_to_par_value(ifc_p1, npt_fc_par1, min_val_fc_par1, max_val_fc_par1)
      fc_grid_p2_val = -999
      if fc_2d:
        fc_grid_p2_val = fc_pt_to_par_value(ifc_p2, npt_fc_par2, min_val_fc_par2, max_val_fc_par2)

      print_str = "For x = " + str(fc_grid_p1_val)
      if fc_2d:
        print_str += ", y = " + str(fc_grid_p2_val)
      print_str += ":"
      print print_str

      # TODO: copy the cpv -pi value onto the +pi value?
      for conf_level_enum in _conf_levels:
        # Calculate the critical Delta chisq
        conf_level = util.ConfLevelToFrac(conf_level_enum)
        critval = get_critical_dchisq(dchisq_vals, conf_level, 0.0)

        # Calculate the errors on the critical value
        err_conf_level = ROOT.TMath.Sqrt(ntoys * conf_level * (1.-conf_level) )
        critval_p1 = get_critical_dchisq(dchisq_vals, conf_level, +err_conf_level)
        critval_m1 = get_critical_dchisq(dchisq_vals, conf_level, -err_conf_level)

        # Store the critical values in a TGraph
        if not fc_2d:
          graphs   [conf_level_enum].SetPoint( graphs   [conf_level_enum].GetN(), fc_grid_p1_val, critval    )
          graphs_p1[conf_level_enum].SetPoint( graphs_p1[conf_level_enum].GetN(), fc_grid_p1_val, critval_p1 )
          graphs_m1[conf_level_enum].SetPoint( graphs_m1[conf_level_enum].GetN(), fc_grid_p1_val, critval_m1 )
        else:
          graphs   [conf_level_enum].SetPoint( graphs   [conf_level_enum].GetN(), fc_grid_p1_val, fc_grid_p2_val, critval    )
          graphs_p1[conf_level_enum].SetPoint( graphs_p1[conf_level_enum].GetN(), fc_grid_p1_val, fc_grid_p2_val, critval_p1 )
          graphs_m1[conf_level_enum].SetPoint( graphs_m1[conf_level_enum].GetN(), fc_grid_p1_val, fc_grid_p2_val, critval_m1 )

        # Print the critical value
        print "\t", util.ConfLevelAsString(conf_level_enum), "critical value =", str(critval), \
                     "+", str(ROOT.TMath.Abs(critval - critval_p1)), \
                     "-", str(ROOT.TMath.Abs(critval - critval_m1))

      input_file.Close()
      del input_file
      del dchisq_vals

  # Sort 1D graphs only
  for conf_level_enum in _conf_levels:
    if not fc_2d:
      graphs   [conf_level_enum].Sort()
      graphs_p1[conf_level_enum].Sort()
      graphs_m1[conf_level_enum].Sort()

  # Save to an output file
  output_file_name = output_tag + ".root"
  print "Saving critical value graphs to file", output_file_name
  output_file = ROOT.TFile(output_file_name, "recreate")
  for conf_level_enum in _conf_levels:
    graphs   [conf_level_enum].Write()
    graphs_p1[conf_level_enum].Write()
    graphs_m1[conf_level_enum].Write()

  output_file.Close()
  del output_file

  print "Done!"

#_______________________________________________________________________
def get_graph_name(conf_level_enum, suffix):
  name = "gFC_" + util.ConfLevelAsString(conf_level_enum)
  if suffix:
    name += "_" + suffix
  return name

#_______________________________________________________________________
def get_filename(tag, fc_2d, ifc_p1, ifc_p2):
  if not fc_2d:
    return tag.replace("[FC]", str(ifc_p1))
  else:
    return tag.replace("[FC]", str(ifc_p1) + "_" + str(ifc_p2))

#_______________________________________________________________________
def get_critical_dchisq(dchisq_vals, conf_level, sigma):
  dchisq_vals.sort()  # In ascending order
  ntoys_conf_level = conf_level * len(dchisq_vals)  # The number of toys contained within this confidence level
  index = int(ntoys_conf_level + sigma)

  if index >= len(dchisq_vals):
    print "ERROR: the index is outside of the Delta chisq list range.\nconf_level =", \
      conf_level, ", sigma =", sigma, ", index =", index, ", length =", len(dchisq_vals), \
      "\n\nDelta chisq values =", dchisq_vals
    sys.exit()

  return dchisq_vals[index]

#_______________________________________________________________________
# Calculate the parameter value corresponding to the FC grid point
def fc_pt_to_par_value(ifc, npt_fc, min_val_fc, max_val_fc):
  step_fc = 0
  if npt_fc > 1:
    step_fc = (max_val_fc - min_val_fc) / (npt_fc - 1)
  if step_fc == 0 and npt_fc > 1:
    print "ERROR: FC step size is 0 with > 1 FC grid point"
    sys.exit()
  return min_val_fc + step_fc * ifc

#_______________________________________________________________________
if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Calculate FC critical values for a single mass hierarchy")

  # Input/Output
  parser.add_option('-o', '--out-tag'     , dest='output_tag', default=""     , help="Tag used to name the output file(s)")
  parser.add_option('-i', '--in-tag'      , dest='in_tag'        ,
                    default="FC_NH_Dchisq_distribution__FCpt_[FC]__10_toys.root",
                    help="Input file tag, with the FC point no. and no. of toys replaced with "
                         "\'[FC]\' and \'[NTOYS]\', respectively")

  # FC grid points
  parser.add_option('--npt-fc-par1'     , dest='npt_fc_par1'     , type=int  , default=9        , help="For param 1: The number of FC grid points to calculate Delta chisq values at")
  parser.add_option('--start-pt-fc-par1', dest='start_pt_fc_par1', type=int  , default=0        , help="For param 1: Start at this FC point number")
  parser.add_option('--min-val-fc-par1' , dest='min_val_fc_par1' , type=float, default=-conv.kpi, help="For param 1: The first FC grid point is at this value")
  parser.add_option('--max-val-fc-par1' , dest='max_val_fc_par1' , type=float, default= conv.kpi, help="For param 1: The last FC grid point is at this value")
  parser.add_option('--npt-fc-par2'     , dest='npt_fc_par2'     , type=int  , default=1        , help="For param 2: The number of FC grid points to calculate Delta chisq values at")
  parser.add_option('--start-pt-fc-par2', dest='start_pt_fc_par2', type=int  , default=0        , help="For param 2: Start at this FC point number")
  parser.add_option('--min-val-fc-par2' , dest='min_val_fc_par2' , type=float, default=0        , help="For param 2: The first FC point is at this value")
  parser.add_option('--max-val-fc-par2' , dest='max_val_fc_par2' , type=float, default=0        , help="For param 2: The last FC point is at this value")

  # Other
  parser.add_option('--fc-2d'       , dest='fc_2d'   , action="store_true"    , help="Specify this flag if you are doing 2D FC")

  o,other = parser.parse_args()

  main(o.output_tag, o.in_tag,
       o.npt_fc_par1, o.start_pt_fc_par1, o.min_val_fc_par1, o.max_val_fc_par1,
       o.npt_fc_par2, o.start_pt_fc_par2, o.min_val_fc_par2, o.max_val_fc_par2,
       o.fc_2d)

  if other:
    print "Warning - ignored arguments: ", other
