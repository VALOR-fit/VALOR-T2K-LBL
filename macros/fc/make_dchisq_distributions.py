#!/usr/bin/env python

#
# For each FC grid point:
#   Collect the marginal chisq values for all FC toys, marginalisation grid points and both mass hierarchies.
#   Then use this information to calculate the Delta chisq values for each FC toy.
#
# Francis Bench <F.Bench \at liverpool.ac.uk>
# August 2019

import ROOT
import sys
import os
import math

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
_chisq_marg_list_name = "chi2_best_fit_toy"
_hier_to_index = dict()
_hier_to_index["normal"]   = 0
_hier_to_index["inverted"] = 1


#_______________________________________________________________________

def main(output_tag, in_tag,
         npt_fc_par1, start_pt_fc_par1, min_val_fc_par1, max_val_fc_par1, pts_fc_par1,
         npt_fc_par2, start_pt_fc_par2, min_val_fc_par2, max_val_fc_par2, pts_fc_par2,
         ntoys_fc, fc_2d, fc_hier, do_brazilian,
         npt_marg_par1, start_pt_marg_par1, min_val_marg_par1, max_val_marg_par1,
         npt_marg_par2, start_pt_marg_par2, min_val_marg_par2, max_val_marg_par2):

  if not fc_2d:
    if npt_fc_par2 != 1 or start_pt_fc_par2 != 0 or npt_marg_par2 != 1 or start_pt_marg_par2 != 0:
      print "WARNING: The option to run 2D FC was not supplied, but the number of points for par 2 was != 1 or " \
            "the start point for par 2was != 0. Setting to defaults."
      npt_fc_par2        = 1
      start_pt_fc_par2   = 0
      npt_marg_par2      = 1
      start_pt_marg_par2 = 0

  # Loop over FC grid points and save the results to a separate file for each
  for ifc_p2 in range(start_pt_fc_par2, start_pt_fc_par2 + npt_fc_par2):
    for ifc_p1 in range(start_pt_fc_par1, start_pt_fc_par1 + npt_fc_par1):

      # Loop over mass hierachy, marginalisation grid points and toys to fill a large chisq vector
      chisq_all = [[[[-1.0]*_n_hier for i in range(npt_marg_par2)] for j in range(npt_marg_par1)] for k in range(ntoys_fc)]

      for hierarchy in ["normal", "inverted"]:
        index_hier = _hier_to_index[hierarchy]

        for imarg_p2 in range(start_pt_marg_par2, start_pt_marg_par2 + npt_marg_par2):
          for imarg_p1 in range(start_pt_marg_par1, start_pt_marg_par1 + npt_marg_par1):

            input_filename = get_filename(in_tag, fc_2d, hierarchy, ifc_p1, ifc_p2, imarg_p1, imarg_p2)
            print "Reading from file", input_filename
            input_file = ROOT.TFile(input_filename, "read")

            chisq_vals = list(input_file.Get(_chisq_marg_list_name))
            if not chisq_vals:
              print "ERROR: Could not read chisq list from file", input_filename
              sys.exit()
            if len(chisq_vals) > ntoys_fc:
              print "ERROR: Tried to read", ntoys_fc, "FC toys from a list of", len(chisq_vals), "toys"
              sys.exit()

            for itoy in range(0, ntoys_fc):
              chisq_all[itoy][imarg_p1 - start_pt_marg_par1][imarg_p2 - start_pt_marg_par2][index_hier] = chisq_vals[itoy]

            input_file.Close()
            del input_file
            del chisq_vals

      # Now that we've filled chisq_all, calculate the FC Delta chisq for each FC toy
      Dchisq_toy = ROOT.vector('double')(ntoys_fc,-1.0)   # Use a ROOT vector instead of a python list so that we can save it to the output ROOT file
       
      if do_brazilian:
        #A vector to store brazilian dchisq
        Dchisq_BR_toy = [[[[-1.0]*_n_hier for i in range(npt_marg_par2)] for j in range(npt_marg_par1)] for k in range(ntoys_fc)]

      index_hier = _hier_to_index[fc_hier]
      index_marg_p1 = fc_pt_to_marg_pt(ifc_p1, npt_fc_par1, min_val_fc_par1, max_val_fc_par1, pts_fc_par1,
                                       npt_marg_par1, start_pt_marg_par1, min_val_marg_par1, max_val_marg_par1)
      index_marg_p2 = fc_pt_to_marg_pt(ifc_p1, npt_fc_par2, min_val_fc_par2, max_val_fc_par2, pts_fc_par2,
                                       npt_marg_par2, start_pt_marg_par2, min_val_marg_par2, max_val_marg_par2)
      for itoy in range(0, ntoys_fc):
        chisq_fix = chisq_all[itoy][index_marg_p1][index_marg_p2][index_hier]
        chisq_bf  = get_chisq_min(chisq_all[itoy])
        Dchisq_toy[itoy] = chisq_fix - chisq_bf
        if do_brazilian:
        #loop over MH, fc toys and marginalization points to calculate dchisq[toy][marg_point][MH] for Brazilian plots
          for hierarchy in ["normal", "inverted"]:
            index_hier = _hier_to_index[hierarchy]
            for imarg_p2 in range(start_pt_marg_par2, start_pt_marg_par2 + npt_marg_par2):
              for imarg_p1 in range(start_pt_marg_par1, start_pt_marg_par1 + npt_marg_par1):
                Dchisq_BR_toy[itoy][imarg_p1][imarg_p2][index_hier] = chisq_all[itoy][imarg_p1][imarg_p2][index_hier] - chisq_bf

      if do_brazilian:
        # Define Delta chisq Brazillian histograms for both MH:
        Dchisq_brazilian_hist_nh = ROOT.TH1D("hDchisq_brazilian_hist_nh", "hDchisq_brazilian_hist_nh", 150, 0, 100)#math.ceil(max(Dchisq_toy)))
        Dchisq_brazilian_hist_ih = ROOT.TH1D("hDchisq_brazilian_hist_ih", "hDchisq_brazilian_hist_ih", 150, 0, 100)#math.ceil(max(Dchisq_toy)))

      # For testing, fill a Delta chisq FC hist
      Dchisq_hist = ROOT.TH1D("hDchisq", "hDchisq", 100, 0, math.ceil(max(Dchisq_toy)))
      for itoy in range(0, ntoys_fc):
        Dchisq_hist.Fill(Dchisq_toy[itoy])
        if do_brazilian:
          #Fill Brazillian dchisq distributions for both MH 
          for hierarchy in ["normal", "inverted"]:
            index_hier = _hier_to_index[hierarchy]
            for imarg_p2 in range(start_pt_marg_par2, start_pt_marg_par2 + npt_marg_par2):
              for imarg_p1 in range(start_pt_marg_par1, start_pt_marg_par1 + npt_marg_par1):
                if index_hier==0:
                  Dchisq_brazilian_hist_nh.Fill(Dchisq_BR_toy[itoy][imarg_p1][imarg_p2][index_hier])
                else:
                  Dchisq_brazilian_hist_ih.Fill(Dchisq_BR_toy[itoy][imarg_p1][imarg_p2][index_hier])


      # Save to an output file
      output_file_name = output_tag + "__FCpt_" + str(ifc_p1)
      if fc_2d:
        output_file_name += "_" + str(ifc_p2)
      output_file_name += "__" + str(ntoys_fc) + "_toys.root"
      print "Saving distribution at FC point (", ifc_p1, ",", ifc_p2, ") to file", output_file_name
      output_file = ROOT.TFile(output_file_name, "recreate")
      output_file.WriteObject(Dchisq_toy , "Dchi2_toy_vec")
      output_file.WriteObject(Dchisq_hist, "Dchi2_hist")
      if do_brazilian:
        output_file.WriteObject(Dchisq_brazilian_hist_nh,"Brazilian_dchi2_NH")
        output_file.WriteObject(Dchisq_brazilian_hist_ih,"Brazilian_dchi2_IH")

      output_file.Close()
      del output_file
      del Dchisq_toy
      del Dchisq_hist

  print "Done!"

#_______________________________________________________________________
def get_filename(tag, fc_2d, hierarchy, ifc_p1, ifc_p2, imarg_p1, imarg_p2):
  input_filename = tag.replace("[HIER]", hierarchy)
  if not fc_2d:
    input_filename = input_filename.replace("[FC]", str(ifc_p1))
  else:
    input_filename = input_filename.replace("[FC]", str(ifc_p1) + "_" + str(ifc_p2))
  if not fc_2d:
    input_filename = input_filename.replace("[MARG]", str(imarg_p1))
  else:
    input_filename = input_filename.replace("[MARG]", str(imarg_p1) + "_" + str(imarg_p2))
  return input_filename

#_______________________________________________________________________
# Find out which marginalisation grid point this FC grid point lies closest to
def fc_pt_to_marg_pt(ifc, npt_fc, min_val_fc, max_val_fc, pts_fc_str, npt_marg, start_pt_marg, min_val_marg, max_val_marg):

  val_fc = None

  if not pts_fc_str:
    # Calculate and check step sizes
    step_fc = 0
    if npt_fc > 1:
      step_fc = (max_val_fc - min_val_fc) / (npt_fc - 1)
    if step_fc == 0 and npt_fc > 1:
      print "ERROR: FC step size is 0 with > 1 FC grid point"
      sys.exit()

    # Calculate the parameter value corresponding to the FC grid point
    val_fc = min_val_fc + step_fc * ifc
  else:
    pts_fc_float = [float(i) for i in pts_fc_str]
    if ifc >= len(pts_fc_float) or ifc < 0:
      print "ERROR: Tried to access custom FC point #\'{}\', but it was out of range. length = {}".format(ifc, len(pts_fc_float))
      sys.exit()
    val_fc = pts_fc_float[ifc]

  step_marg = 0
  if npt_marg > 1:
    step_marg = (max_val_marg - min_val_marg) / (npt_marg - 1)
  if step_marg == 0 and npt_marg > 1:
    print "ERROR: Marginalisation step size is 0 with > 1 marginalisation grid point"
    sys.exit()
  if step_marg == 0 and npt_marg == 1:
    return start_pt_marg

  # Find the marginalisation grid point that this is closest to
  for imarg in range(start_pt_marg, start_pt_marg + npt_marg):
    marg_grid_pt_val = min_val_marg + step_marg * imarg
    low  = marg_grid_pt_val - 0.5 * step_marg
    high = marg_grid_pt_val + 0.5 * step_marg
    if val_fc >= low and val_fc < high:
      return imarg

  print "ERROR: Could not find the marginalisation grid point closest to the FC grid point (at", val_fc, ")"
  sys.exit()
  
#_______________________________________________________________________
# The input param chisq_list_toy should be a 3D list with indices [index_marg_p1][index_marg_p2][index_hier]
def get_chisq_min(chisq_list_toy):
  min = 99999999.9
  for i in range(0, len(chisq_list_toy)):
    for j in range(0, len(chisq_list_toy[i])):
      for k in range(0, len(chisq_list_toy[i][j])):
        if chisq_list_toy[i][j][k] < min:
          min = chisq_list_toy[i][j][k]

  # Check to see if an element was unfilled (as they're initialised to -1)
  if min < 0:
    print "ERROR: A chisq list element was not filled!"
    sys.exit()

  return min

#_______________________________________________________________________
if __name__=="__main__":
  
  import argparse
  parser = argparse.ArgumentParser(usage="%prog [options]..",description="Calculate FC Delta chisq values")

  # Input/Output
  parser.add_argument('-o', '--out-tag'     , dest='output_tag', default=""     , help="Tag used to name the output file(s)")
  parser.add_argument('-i', '--in-tag'      , dest='in_tag'        ,
                    default="fit_run1-9d_NH_FC_[HIER]_cpv_reactor__FCpt_[FC]__margpt_[MARG]_result.root",
                    help="Input file tag, with the hierarchy, FC point no. and marginalisation point no. replaced with "
                         "\'[HIER]\', \'[FC]\' and \'[MARG]\', respectively")

  # FC grid points
  parser.add_argument('--npt-fc-par1'     , dest='npt_fc_par1'     , type=int  , default=9        , help="For param 1: The number of FC grid points to calculate Delta chisq values at")
  parser.add_argument('--start-pt-fc-par1', dest='start_pt_fc_par1', type=int  , default=0        , help="For param 1: Start at this FC point number")
  parser.add_argument('--min-val-fc-par1' , dest='min_val_fc_par1' , type=float, default=-conv.kpi, help="For param 1: The first FC grid point is at this value")
  parser.add_argument('--max-val-fc-par1' , dest='max_val_fc_par1' , type=float, default= conv.kpi, help="For param 1: The last FC grid point is at this value")
  parser.add_argument("--pts-fc-par1"     , dest="pts_fc_par1"     , nargs="+" , default=None     , help="For param 1: If the FC grid points aren't evenly spaced, specify them all here instead of using --min... and --max...")
  parser.add_argument('--npt-fc-par2'     , dest='npt_fc_par2'     , type=int  , default=1        , help="For param 2: The number of FC grid points to calculate Delta chisq values at")
  parser.add_argument('--start-pt-fc-par2', dest='start_pt_fc_par2', type=int  , default=0        , help="For param 2: Start at this FC point number")
  parser.add_argument('--min-val-fc-par2' , dest='min_val_fc_par2' , type=float, default=0        , help="For param 2: The first FC point is at this value")
  parser.add_argument('--max-val-fc-par2' , dest='max_val_fc_par2' , type=float, default=0        , help="For param 2: The last FC point is at this value")
  parser.add_argument("--pts-fc-par2"     , dest="pts_fc_par2"     , nargs="+" , default=None     , help="For param 2: If the FC grid points aren't evenly spaced, specify them all here instead of using --min... and --max...")

  # Marginalisation grid points
  parser.add_argument('--npt-marg-par1'     , dest='npt_marg_par1'     , type=int  , default=101      , help="For param 1: The number of marginalisation grid points to calculate Delta chisq values at")
  parser.add_argument('--start-pt-marg-par1', dest='start_pt_marg_par1', type=int  , default=0        , help="For param 1: Start at this marginalisation point number")
  parser.add_argument('--min-val-marg-par1' , dest='min_val_marg_par1' , type=float, default=-conv.kpi, help="For param 1: The first marginalisation grid point is at this value")
  parser.add_argument('--max-val-marg-par1' , dest='max_val_marg_par1' , type=float, default= conv.kpi, help="For param 1: The last marginalisation grid point is at this value")
  parser.add_argument('--npt-marg-par2'     , dest='npt_marg_par2'     , type=int  , default=1        , help="For param 2: The number of marginalisation grid points to calculate Delta chisq values at")
  parser.add_argument('--start-pt-marg-par2', dest='start_pt_marg_par2', type=int  , default=0        , help="For param 2: Start at this marginalisation point number")
  parser.add_argument('--min-val-marg-par2' , dest='min_val_marg_par2' , type=float, default=0        , help="For param 2: The first marginalisation grid point is at this value")
  parser.add_argument('--max-val-marg-par2' , dest='max_val_marg_par2' , type=float, default=0        , help="For param 2: The last marginalisation grid point is at this value")

  # Other
  parser.add_argument('--do-brazilian', dest='do_brazilian'   , default=False          , action="store_true"    , help="Specify this flag if you want to calculate dcisq for brazilian plot")
  parser.add_argument('--ntoys-fc'    , dest='ntoys_fc'       , type=int, default=10000, help="The number of FC toys at each grid point")
  parser.add_argument('--fc-2d'       , dest='fc_2d'          , action="store_true"    , help="Specify this flag if you are doing 2D FC")
  parser.add_argument('--fc-hierarchy', dest='fc_hier'        , default="normal"       , help="The mass hierarchy that the FC toys were generated at")
  
 
  o = parser.parse_args()

  if not(o.fc_hier == "normal" or o.fc_hier == "inverted"):
    print "ERROR: The --fc-hierarchy option must be either \'normal\' or \'inverted\'"
    sys.exit()

  main(o.output_tag, o.in_tag,
       o.npt_fc_par1, o.start_pt_fc_par1, o.min_val_fc_par1, o.max_val_fc_par1, o.pts_fc_par1,
       o.npt_fc_par2, o.start_pt_fc_par2, o.min_val_fc_par2, o.max_val_fc_par2, o.pts_fc_par2,
       o.ntoys_fc, o.fc_2d, o.fc_hier, o.do_brazilian,
       o.npt_marg_par1, o.start_pt_marg_par1, o.min_val_marg_par1, o.max_val_marg_par1,
       o.npt_marg_par2, o.start_pt_marg_par2, o.min_val_marg_par2, o.max_val_marg_par2)