#!/usr/bin/env python

# Macro to disable specified parameters from the input BANFF matrix

import ROOT
import collections

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor = ROOT.valor
util  = valor.util
conv  = valor.conv


# To disable a parmeter, add it's name to the floowing list
disable_params = [
    "ebo"
]

# What the input file should contain:

#KEY: TObjArray	param_list;1	An array of objects
#  KEY: TMatrixT<double>	prefit_cov;1	Template of General Matrix class
#  KEY: TVectorT<double>	prefit_params;1	Template of Vector class
#  KEY: TMatrixT<double>	postfit_cov;1	Template of General Matrix class
#  KEY: TVectorT<double>	postfit_params;1	Template of Vector class
#  KEY: TVectorT<double>	banff_untuned_central_values;1	Template of Vector class
#  KEY: TObjArray	param_list_banff_xsec;1	An array of objects
#  KEY: TVectorT<double>	banff_xsec_untuned_central_values;1	Template of Vector class


def main(infilename, outfilename):
    banff_file = ROOT.TFile(infilename,"read")
    banff = banff_file.Get("postfit_cov")
    param_names = banff_file.Get("param_list")
    prefit_cov = banff_file.Get("prefit_cov")     # TODO:test
    prefit_params = banff_file.Get("prefit_params")     # TODO:test
    postfit_params = banff_file.Get("postfit_params")     # TODO:test
    banff_untuned_central_values = banff_file.Get("banff_untuned_central_values")     # TODO: test
    param_list_banff_xsec = banff_file.Get("param_list_banff_xsec")     # TODO: test
    banff_xsec_untuned_central_values = banff_file.Get("banff_xsec_untuned_central_values")     # TODO: test

    nrows = banff.GetNrows()
    ncols = banff.GetNcols()

    output_param_names = ROOT.TObjArray(param_names)
    output_postfit_cov = ROOT.TMatrixD(banff)
    output_prefit_cov = ROOT.TMatrixD(prefit_cov)
    output_prefit_params = ROOT.TVectorD(prefit_params)
    output_postfit_params = ROOT.TVectorD(postfit_params)
    output_param_list_banff_xsec = ROOT.TObjArray(param_list_banff_xsec)
    output_banff_untuned_central_values = ROOT.TVectorD(banff_untuned_central_values)
    output_banff_xsec_untuned_central_values = ROOT.TVectorD(banff_xsec_untuned_central_values)

    if ncols != nrows:
        print "Input matricies must be square! Exiting"
        util.Die(1)

    
    # param loop
    for par_name in disable_params:
        # Get the position of the parameter of interest
        par_index = -1
        for i in range(0, param_names.GetEntriesFast()):
            param_name = param_names.At(i).String().Data()
            #print "*** test: par_name =", par_name, "compare =", param_name
            if param_name == par_name:
                par_index = i
                break
        
        # Disable param in postfit and prefit covariance matricies
        for i in range(0, nrows):
            for j in range(0, nrows):
                if i == par_index or j == par_index:
                    if i == j:
                        output_postfit_cov[i][j] = 1.0
                        output_prefit_cov[i][j] = 1.0
                    else:
                        output_postfit_cov[i][j] = 0.0
                        output_prefit_cov[i][j] = 0.0
        
        # Disable param in prefit and postfit params vectors
        output_prefit_params[par_index] = 0.0
        output_postfit_params[par_index] = 0.0

        # Disable param in xsec vector
        par_index_xsec = -1
        for i in range(0, param_list_banff_xsec.GetEntriesFast()):
            param_name = param_list_banff_xsec.At(i).String().Data()
            if param_name == par_name:
                par_index_xsec = i
                break
        output_banff_xsec_untuned_central_values[par_index_xsec] = 0.0
                    
    banff_file.Close()
    outfile = ROOT.TFile(outfilename,"recreate")
    outfile.cd()
    output_param_names.Write("param_list",ROOT.TObject.kSingleKey)
    output_postfit_cov.Write("postfit_cov")
    output_prefit_cov.Write("prefit_cov")
    output_prefit_params.Write("prefit_params")
    output_postfit_params.Write("postfit_params")
    output_param_list_banff_xsec.Write("param_list_banff_xsec",ROOT.TObject.kSingleKey)
    output_banff_untuned_central_values.Write("banff_untuned_central_values")
    output_banff_xsec_untuned_central_values.Write("banff_xsec_untuned_central_values")
    outfile.Close()
    print "Done!"
    


if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Disable parameters in the input BANFF matrix (based on param names specified in the macro)")
    parser.add_option("--input", dest="inp", default="", help='Input file containing a BANFF matrix')
    parser.add_option("--output" , dest="out", default="", help='Output file name')

    v = optparse.Values()
    parser.parse_args(values = v)   # v contains only the values that were specified, not defaults
    o, other = parser.parse_args()

    main(o.inp, o.out)

    if other:
        print "Warning - ignored arguments: ", other
