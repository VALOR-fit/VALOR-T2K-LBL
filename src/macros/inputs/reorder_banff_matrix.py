#!/usr/bin/env python

# Re-order the BANFF matrix supplied by the BANFF group, using an ordered list of parameters taken from an XML file.
#
# NOTE: Don't use this for the SKDetFSI matrix. This macro assumes that the non-VALOR format input file
# contains an ordered list of parameter names, but unhelpfully, the SDKDetFSI file does not contain this.
#
# Francis Bench
# Sept. 2019


import sys
import ROOT
import re
import collections
import numpy

# Don't let ROOT ruin --help
ROOT.PyConfig.IgnoreCommandLineOptions = True

# load valor
env_var = "VALOR_EXPT"
try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $VALOR_EXPT/root/valor_expt.C")
  
ROOT.gROOT.SetBatch(1)
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()
rnd  =  valor.Rnd.Instance()

python_undefined = 2L*sys.maxint+1
python_undefined_short = 2**32 -1

# List strange untuned central values. Will assume 1 for all others.
# Want to grab these from the TVector of untuned values instead
abnormal_untuned_central_values = {
    "MAQE"                     : 1.21,
    "MARES"                    : 0.95,
    "CA5"                      : 1.01,
    "2p2h_shape_O"             : 0.0,
    "FEFQE"                    : 1.069,
    "FEFQEH"                   : 1.824,
    "FEFINEL"                  : 1.002,
    "FEFABS"                   : 1.404,
    "FEFCX"                    : 0.697,
    "ISO_BKG"                  : 1.3,
    "ISO_BKG_LowPPi"           : 1.3,
    "EB_dial_O_nu"             : 0.0,
    "EB_dial_O_nubar"          : 0.0,
    "CC_BY_DIS"                : 0.0,
    "CC_BY_MPi"                : 0.0,
    "CC_AGKY_Mult"             : 0.0,
    "Q2_norm_0"                : 0.0,
    "Q2_norm_1"                : 0.0,
    "Q2_norm_2"                : 0.0,
    "Q2_norm_3"                : 0.0,
    "Q2_norm_4"                : 0.0
}


def main(infile_name, parameter_config_filename, parameter_list_filename, del_not_in_list):

  # First read the parameter config XML file so we know how to categorise the parameters
  valor.ParameterPool.MutableInstance().ReadXML(parameter_config_filename)

  # Read the parameter list XML file to determine the order of parameters in the output matrix
  ordered_param_name_dict = get_ordered_param_names(parameter_list_filename)

  # Set up files
  assert(infile_name)

  for beam_mode_name in ["fhc","rhc","all"]:

    print "\n\n===================================="
    print "Beam mode:", beam_mode_name
    print "====================================\n"

    outfile_name = infile_name.replace(".root","_valor_order_"+beam_mode_name+".root")
    print "Will write to file:", outfile_name

    infile  = ROOT.TFile(infile_name)
    outfile = ROOT.TFile(outfile_name,"RECREATE")

    # Get the input parameters and count their types
    input_param_list = [str(s) for s in infile.Get("param_list")]
    n_banff_params = len(input_param_list)

    n_banff_flux_params = 0
    for s in input_param_list:
      if "SKNuMode" in s or "SKANuMode" in s:
        n_banff_flux_params += 1

    n_banff_xsec_params = n_banff_params - n_banff_flux_params

    check_param_list(ordered_param_name_dict, input_param_list, n_banff_params, n_banff_flux_params, n_banff_xsec_params, del_not_in_list)

    params_to_delete = [s for s in numpy.setdiff1d(input_param_list, ordered_param_name_dict.values(), assume_unique=True)]
    if not del_not_in_list and len(params_to_delete) > 0:
      print "ERROR: Did not specify that we should delete params, but the list of params to delete is non-empty"
      sys.exit()

    print "Original parameter order:\n\t", input_param_list
    print "New parameter order:\n\t", ordered_param_name_dict.values()
    print "Deleting the following params from the input BANFF matrix:\n\t", params_to_delete

    # The index of param_order_list is the NEW position of each parameter, the entry is the OLD position
    param_order_list = [input_param_list.index(ordered_param_name_dict[valor_name]) for new_i, valor_name in enumerate(ordered_param_name_dict)]

    # Save the reordered VALOR format names to the output list
    output_names = ROOT.TObjArray()
    for valor_name in ordered_param_name_dict.keys():
      output_names.Add(ROOT.TObjString(valor_name))
    output_names.Write("param_list",ROOT.TObject.kSingleKey)

    n_elements = len(param_order_list)
    print "Elements in new matrix", n_elements

    reordered_covs_d = dict()
    reordered_params_d = dict()

    # Reorder the parameters and the covariance
    for fit in ["prefit","postfit"]:
      orig_cov    = infile.Get(fit+"_cov")
      orig_params = infile.Get(fit+"_params")

      reordered_cov    = ROOT.TMatrixD(n_elements,n_elements)
      reordered_params = ROOT.TVectorD(n_elements)
      for i,orig_i in enumerate(param_order_list):
        reordered_params[i] = orig_params[orig_i]
        for j,orig_j in enumerate(param_order_list):
          reordered_cov[i][j] = orig_cov[orig_i][orig_j]

      outfile.cd()
      reordered_cov   .Write(fit+"_cov")
      reordered_params.Write(fit+"_params")

      reordered_covs_d  [fit] = reordered_cov
      reordered_params_d[fit] = reordered_params


    # Produce untuned central values - all
    vec = ROOT.TVectorD(len(ordered_param_name_dict))

    for i in xrange(vec.GetNrows()):
      vec[i] = 1.
    for name, val in abnormal_untuned_central_values.iteritems():
      try:
        ind = list(output_names).index(name)
        vec[ind] = val
      except ValueError:
        print "Making central values - failed to set value of ",name

    vec.Write("banff_untuned_central_values")


    # Produce untuned central values - xsec only
    xsec_vec = ROOT.TVectorD(n_banff_xsec_params)
    xsec_output_names = ROOT.TObjArray(output_names)
    while xsec_output_names.GetEntries() > n_banff_xsec_params:
      xsec_output_names.RemoveFirst()
      xsec_output_names.Compress()
    xsec_output_names.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)

    for i in xrange(xsec_vec.GetNrows()):
      #xsec_vec[i] =
      xsec_vec[i] = 1.
    for name, val in abnormal_untuned_central_values.iteritems():
      try:
        ind = list(xsec_output_names).index(name)
        xsec_vec[ind] = val
      except ValueError:
        print "Making central values - failed to set value of ",name

    xsec_vec.Write("banff_xsec_untuned_central_values")

    infile.Close()
    outfile.Close()


#_______________________________________________________________________
def get_ordered_param_names(parameter_list_filename):
  ordered_param_name_dict = collections.OrderedDict()

  # Parse syst list file and put parameters into the list
  xml_syst = ROOT.TXMLEngine()
  doc_syst = xml_syst.ParseFile(util.CompleteConfigPath(parameter_list_filename))
  root_node_syst  = xml_syst.DocGetRootElement(doc_syst)
  curr_child_node = xml_syst.GetChild(root_node_syst)

  while curr_child_node:
    content = util.RemoveAllWhitespace(xml_syst.GetNodeContent(curr_child_node))
    p = util.NameRegistry.Instance().PhysParamFromString(content)
    if p == python_undefined or p == python_undefined_short:
      print "ERROR: Could not recognise parameter with name", content, "in the input parameter list", parameter_list_filename
      sys.exit()
    param = valor.ParameterPool.Instance().GetParam(p)
    valor_name = param.GetName()
    banff_name = get_banff_name(param)
    ordered_param_name_dict[valor_name] = banff_name
    curr_child_node = xml_syst.GetNext(curr_child_node)

  xml_syst.FreeDoc(doc_syst)
  del xml_syst

  if len(ordered_param_name_dict.keys()) != len(list(set(ordered_param_name_dict.keys()))):
    print "ERROR: Duplicates found in the input parameter list", parameter_list_filename
    sys.exit()

  return ordered_param_name_dict


#_______________________________________________________________________
# Convert between VALOR param names and BANFF param names
def get_banff_name(param):
  valor_name = param.GetName()
  banff_name = valor_name

  if param.SearchCategories("flux"):
    banff_name = "SK"

    if "_rhc" in valor_name:
      banff_name += "A"
    banff_name += "NuMode"

    if "numubar" in valor_name:
      banff_name += "Numub"
    elif "numu" in valor_name:
      banff_name += "Numu"
    elif "nuebar" in valor_name:
      banff_name += "Nueb"
    elif "nue" in valor_name:
      banff_name += "Nue"
    else:
      print "ERROR: Could not parse VALOR flux param name"
      sys.exit()

    param_number = re.findall(r'\d+', valor_name)
    if len(param_number) != 1:
      print "ERROR: Could not parse VALOR flux param name"
      sys.exit()
    banff_name += str(int(param_number[0]))

  return banff_name


#_______________________________________________________________________
def check_param_list(ordered_param_name_dict, input_banff_param_list, n_banff_params, n_banff_flux_params, n_banff_xsec_params, del_not_in_list):
  # Does our input XML list have the right total number of params?
  if not del_not_in_list and n_banff_params != len(ordered_param_name_dict):
    print "ERROR: The number of parameters in the input XML param list (", len(ordered_param_name_dict), ") " \
          "does not match the number of parameters in the input matrix (", n_banff_params, ")"
    print "input XML params:\n", ordered_param_name_dict.keys()
    print "input matrix params:\n", input_banff_param_list
    sys.exit()

  # Does each converted param name from the XML list exist in the param list from the BANFF matrix file?
  for param_banff_name in ordered_param_name_dict.values():
    if param_banff_name not in input_banff_param_list:
      print "ERROR: Converted BANFF param param name", param_banff_name, "does not exist in the banff param list"
      sys.exit()

  # Do they all exist in the ParameterPool? Also check that we have the right number of parameters of each type
  n_banff_params_pool = 0
  n_banff_flux_params_pool = 0
  n_banff_xsec_params_pool = 0
  for param_valor_name in ordered_param_name_dict.keys():
    p = util.NameRegistry.Instance().PhysParamFromString(param_valor_name)
    if p == python_undefined or p == python_undefined_short:
      print "ERROR: Could not recognise parameter with name", param_valor_name
      sys.exit()
    if not valor.ParameterPool.Instance().IsParamInPool(p):
      print "ERROR: The parameter ", param_valor_name, " exists in the input list but not in the ParameterPool"
      sys.exit()
    param = valor.ParameterPool.Instance().GetParam(p)
    if param.SearchCategories("BANFF"):
      n_banff_params_pool += 1
      if param.SearchCategories("flux"): n_banff_flux_params_pool += 1
      if param.SearchCategories("int" ): n_banff_xsec_params_pool += 1

  if not del_not_in_list and n_banff_params_pool != n_banff_params:
    print "ERROR: Found", n_banff_params_pool, "BANFF params in the parameter pool, expected", n_banff_params
    sys.exit()
  if not del_not_in_list and n_banff_flux_params_pool != n_banff_flux_params:
    print "ERROR: Found", n_banff_flux_params_pool, "BANFF flux params in the parameter pool, expected", n_banff_flux_params
    sys.exit()
  if not del_not_in_list and n_banff_xsec_params_pool != n_banff_xsec_params:
    print "ERROR: Found", n_banff_xsec_params_pool, "BANFF xsec params in the parameter pool, expected", n_banff_xsec_params
    sys.exit()

  print "Num of raw BANFF params =", n_banff_params, "; Num of VALOR BANFF params =", n_banff_params_pool
  print "Num of raw BANFF flux params =", n_banff_flux_params, "; Num of VALOR BANFF flux params =", n_banff_flux_params_pool
  print "Num of raw BANFF xsec params =", n_banff_xsec_params, "; Num of VALOR BANFF xsec params =", n_banff_xsec_params_pool

  if n_banff_flux_params % 2 != 0:
    print "ERROR: Expected an even number of raw BANFF flux params"
    sys.exit()

  if n_banff_flux_params_pool % 2 != 0:
    print "ERROR: Expected an even number of VALOR BANFF flux params"
    sys.exit()


#_______________________________________________________________________
if __name__=="__main__":

  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Reformat and optionally reduce input SK MC trees")

  parser.add_option('-i', '--in-file'           , dest='infile_name'          , default=""                        , help="Input file name.")
  parser.add_option('-p', '--param-config'      , dest='param_config_filename', default=""                        , help="Config file specifying the systematic parameters and their categories. All BANFF parameters must be present.")
  parser.add_option('-l', '--param-list'        , dest='param_list_filename'  , default=""                        , help="Config file specifying the list of systematic parameters. The order of params in this list determines the order of the output matrix.")
  parser.add_option('-d', '--delete-not-in-list', dest='del_not_in_list'      , default=False, action="store_true", help="Delete params from the matrix (and associated objects) that do not appear in the input --param-list")

  o,other = parser.parse_args()

  main(o.infile_name, o.param_config_filename, o.param_list_filename, o.del_not_in_list)

  if other:
    print "Warning - ignored arguments: ", other
