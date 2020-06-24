#!/usr/bin/python

import ROOT
import os
import sys
from collections import namedtuple
import re

fmwk_var = "VALOR_FMWK"
sys.path.append(os.path.abspath(os.environ[fmwk_var] + "/src/valor/Utilities/"))
from UtlPy import c_to_pyroot_array_hack

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

modelist_names_true = [
    "rightsign_numu+oscnue_modes",
    "wrongsign_numu+oscnue_modes",
    "rightsign_nue_modes",
    "wrongsign_nue_modes",
]

modelist_names_reco = [
    "mulike",
    "elike"
]


#______________________________________________________________________________________________________________________________

# TODO: Do we have more than one combination of kinevars for xsec binning?

def main(tag, overwrite):
    #nbins_etrue_all = []
    #nbins_theta_all = []
    edges_etrue_all = []
    edges_ereco_all = []
    edges_theta_all = []
    modes_with_same_binning_reco = []
    modes_with_same_binning_true = []

    # Store the xsec binning schemes
    for sample_index in range(0, conv.kns):
        s = util.GetSKSample(sample_index)
        if s == conv.kSKsFCFVMultiRingMuLike or s == conv.kSKsUndefined or s == 2**64 -1:  # kSKsUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
            # Don't yet have multi-ring mu-like samples
            continue
        nbins_theta = util.NumXsecQ2Bins(s)
        edges_theta = c_to_pyroot_array_hack(util.XsecQ2BinEdges(s), nbins_theta+1) # Need to set the size of the read-write buffer, otherwise pyroot thinks it has 2**31 -1 entries
        nbins_ereco = conv.kNXsErecoBins[sample_index]
        edges_ereco = c_to_pyroot_array_hack(util.XsecErecoBinEdges(s), nbins_ereco+1)
        for mode_index in range(0, conv.knm):
            m = util.GetMode(sample_index, mode_index)
            if not (m == conv.kSKmUndefined or m == 2**64 -1):  # kSKmUndefined = -1, but pyroot converts this to an unsigned int, so is 2^61 - 1
                nbins_etrue = util.NumVariableXsEtrueBins(s,m)
                edges_etrue = c_to_pyroot_array_hack(util.VariableXsEtrueBinEdges(s,m), nbins_etrue+1)

                print "For sample", util.SampleAsString(s), "and mode", util.ModeAsString(m), \
                    ": \n  nbins theta =", nbins_theta, ", nbins etrue =", nbins_etrue, ", nbins ereco =", nbins_ereco, \
                    "\n  edges theta =", edges_theta, \
                    "\n  edges etrue =", edges_etrue, \
                    "\n  edges ereo =" , edges_ereco

                binning_index_reco = get_common_binning_index(edges_theta_all, edges_theta, edges_ereco_all, edges_ereco)
                binning_index_true = get_common_binning_index(edges_etrue_all, edges_etrue)
                if binning_index_reco == -1:
                    # We haven't saved this binning scheme yet, so save it and say what mode it corresponds to
                    print "Recording this reco binning scheme"
                    edges_theta_all.append(edges_theta)
                    edges_ereco_all.append(edges_ereco)
                    modes_with_same_binning_reco.append([[s, m]])
                else:
                    # We have already recorded this binning scheme, so just record that this mode uses this binning.
                    # To save the mode in the right position, need both the current etrue and xsec binning to match a
                    # previously added binning scheme
                    print "Reco binning scheme already exists, adding mode to list"
                    modes_with_same_binning_reco[binning_index_reco].append([s,m])
                if binning_index_true == -1:
                    print "Recording this true binning scheme"
                    edges_etrue_all.append(edges_etrue)
                    modes_with_same_binning_true.append([[s, m]])
                else:
                    print "True binning scheme already exists, adding mode to list"
                    modes_with_same_binning_true[binning_index_true].append([s,m])
                print "\n\n"

    if not ( len(modes_with_same_binning_reco) == len(edges_theta_all)
            and len(edges_ereco_all) == len(edges_theta_all) ):
        print "ERROR: Number of binning schemes recorded for each reco kinevar / number of mode lists. Exiting."
        util.Die(-1)
    if not ( len(modes_with_same_binning_true) == len(edges_etrue_all) ):
        print "ERROR: Number of binning schemes recorded for each true kinevar / number of mode lists. Exiting."
        util.Die(-1)


    # Write the Xsec binning schemes to XML files
    for ischeme in range (0, len(edges_theta_all)):
        mode_list_str = get_mode_sample_list_str(modes_with_same_binning_reco[ischeme])
        fname_reco_e2e = tag + "_e2e_Ereco_Thetareco_" + modelist_names_reco[ischeme] + ".xml"
        write_binning_xml_content(edges_ereco_all[ischeme], edges_theta_all[ischeme],
                                  mode_list_str, fname_reco_e2e, overwrite, False, True)
        fname_reco_mondrian = tag + "_mondrian_Ereco_Thetareco_" + modelist_names_reco[ischeme] + ".xml"
        write_binning_xml_content(edges_ereco_all[ischeme], edges_theta_all[ischeme],
                                  mode_list_str, fname_reco_mondrian, overwrite, False, False)
    for ischeme in range (0, len(edges_etrue_all)):
        mode_list_str = get_mode_sample_list_str(modes_with_same_binning_true[ischeme])
        fname_true = tag + "_e2e_Etrue_" + modelist_names_true[ischeme] + ".xml"
        write_binning_xml_content(edges_etrue_all[ischeme], [],
                                  mode_list_str, fname_true, overwrite, True, True)

    print "Done!"
    return



def get_mode_sample_list_str(sample_mode_pair_list):
    mode_list_str = ""
    for pair in sample_mode_pair_list:
        mode_list_str += util.SampleAsString(pair[0]) + "\t\t" + util.ModeAsString(pair[1]) + "\n"
    return mode_list_str



def get_common_binning_index(edges_all1, edges1, edges_all2 = [], edges2 = []):
    indicies1 = get_indicies_of_matching_list(edges_all1, edges1)
    indicies2 = get_indicies_of_matching_list(edges_all2, edges2)
    if len(indicies1) == 0 or (len(indicies2) == 0 and len(edges2) != 0):
        return -1   # This binning scheme has not yet been recorded
    common_index = -1
    found = False
    if len(edges2) != 0:
        for index1 in indicies1:
            for index2 in indicies2:
                if index1 == index2:
                    if not found:
                        common_index = index1
                        found = True
                    else:
                        print "ERROR: Should not be more than one common index. Exiting."
                        util.Die(-1)
    elif len(indicies1) > 1:
        print "ERROR: Should not be more than one common index. Exiting."
        util.Die(-1)
    else: common_index = indicies1[0]

    return common_index



def get_indicies_of_matching_list(edges_all, edges_search):
    matching_indicies = []
    for i in range(0, len(edges_all)):
        if edges_all[i] == edges_search:
            matching_indicies.append(i)
    return matching_indicies



def write_binning_xml_content(edges_energy, edges_theta, mode_list_str, filename, overwrite, is_true, is_e2e):
    print "Writing binning to XML file", filename, ", which applied to the following modes:\n", mode_list_str
    exists = os.path.isfile(filename)
    if exists and not overwrite:
        print "ERROR: Output file", filename, " already exists and the overwrite option was not specified. Exiting."
        util.Die(-1)
    f = open(filename, "w")

    f.write( "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n\n")
    #f.write( "<!--\nTODO: Make sure this xsec binning file is only applied to the SDK equivalent of the following legacy T2K modes:\n")
    f#.write( mode_list_str + "-->\n\n")
    f.write( "<binning>\n\n")
    if is_e2e:
        f.write( "\t<type>\n\t\tedge-to-edge\n\t</type>\n\n")
        if is_true:
            f.write( "\t<space>\n\t\tEtrue\n\t</space>\n\n")
            write_single_kinevar_e2e_binning("Etrue", edges_energy, f)
        else:
            f.write( "\t<space>\n\t\tEreco,Thetareco\n\t</space>\n\n")
            write_single_kinevar_e2e_binning("Ereco", edges_energy, f)
            write_single_kinevar_e2e_binning("Thetareco", edges_theta, f)
    else:
        f.write( "\t<type>\n\t\tmondrian\n\t</type>\n\n")
        nbins = ( len(edges_energy) -1 ) * ( len(edges_theta) -1 )
        if nbins < 1:
            print "ERROR: < 1 Mondrian bins. Exiting."
            util.Die(-1)
        f.write( "\t<nbins>\n\t\t"+str(nbins)+"\n\t</nbins>\n\n")

        # Treat mu-like and e-like binning separately
        if len(edges_theta) == 2:
            # mu-like
            for i in range(0, len(edges_energy) -1):
                f.write( "\t<bin_" + str(i) + ">  " + "{:.2f}".format(edges_energy[i]) + ",  " + \
                         "{:.2f}".format(edges_energy[i+1]) + "\t</bin_" + str(i) + ">\n" )
        else:
            # e-like
            for j in range(0, len(edges_theta) -1):
                for i in range(0, len(edges_energy) -1):
                    bindex = j * (len(edges_energy) -1) + i
                    f.write( "\t<bin_" + str(bindex) + ">  " + "{:.2f}".format(edges_energy[i]) + ",  " + \
                          "{:.2f}".format(edges_energy[i+1]) + ",  " + "{:.0f}".format(edges_theta[j]) + ",  " + \
                          "{:.0f}".format(edges_theta[j+1]) + "\t</bin_" + str(bindex) + ">\n" )

    f.write( "</binning>\n")



def write_single_kinevar_e2e_binning(name, edges, f):
    print_line_length = 5
    f.write( "\t<"+name+"_binning>\n")
    f.write( "\t\t<bin_size> variable </bin_size>\n")
    f.write( "\t\t<nbins> " + str(len(edges)-1) + " </nbins>\n")
    f.write( "\t\t<bin_edges>\n")
    edges_str = "\t\t\t"
    for i in range(0, len(edges)):
        edges_str += str(edges[i])
        if i != len(edges) -1:
            edges_str += ",\t"
            if i != 0 and i % print_line_length == 0:
                edges_str += "\n\t\t\t"
    f.write( edges_str + "\n")
    f.write( "\t\t</bin_edges>\n")
    f.write( "\t</"+name+"_binning>\n\n")



if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Create an SDK binning XML file for T2K xsec binning schemes.")
  parser.add_option('-t', '--tag', dest='tag', default="t2k_xsec_binning", help="Tag used in naming the output XML files")
  parser.add_option('-o', '--overwrite', dest='overwrite', action="store_true", help="If this option is specified, any existing file with the same output name will be deleted and overwritten.")

  o,other = parser.parse_args()
  main(o.tag, o.overwrite)
  
  if other:
    print "Warning - ignored arguments: ",other
