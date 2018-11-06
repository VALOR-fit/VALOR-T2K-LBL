#!/usr/bin/env python

# Macro to compare different BANFF matricies and re-order the first matrix according to the positions of the parameters
# in the second matrix

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


# To append new parameters to the end of the new matrix, add them to the dictionary below with the value set to a boolean where True means that it's an xsec param. 
#At the moment, all covariance values for these new parameters is set to 0.
append_params = collections.OrderedDict([
    ("CC_norm_nu"    , True),
    ("CC_norm_nubar" , True),
    ("FEFABS"        , True),
    ("FEFQE"         , True),
    ("FEFCX"         , True),
    ("FEFQEH"        , True),
    ("FEFINEL"       , True)
])

# If a parameter needs to be renamed in the new matrix, add them to the dictionary below in the format
# <old name> : <new name>
rename = {
    "test_entry" : "test_entry"
}

#KEY: TObjArray	param_list;1	An array of objects
#  KEY: TMatrixT<double>	prefit_cov;1	Template of General Matrix class
#  KEY: TVectorT<double>	prefit_params;1	Template of Vector class
#  KEY: TMatrixT<double>	postfit_cov;1	Template of General Matrix class
#  KEY: TVectorT<double>	postfit_params;1	Template of Vector class
#  KEY: TVectorT<double>	banff_untuned_central_values;1	Template of Vector class
#  KEY: TObjArray	param_list_banff_xsec;1	An array of objects
#  KEY: TVectorT<double>	banff_xsec_untuned_central_values;1	Template of Vector class


def main(infilename1, infilename2, outfilename):
    banff_file1 = ROOT.TFile(infilename1,"read")
    banff_file2 = ROOT.TFile(infilename2,"read")
    banff1 = banff_file1.Get("postfit_cov")
    banff2 = banff_file2.Get("postfit_cov")
    param_names1 = banff_file1.Get("param_list")
    param_names2 = banff_file2.Get("param_list")
    prefit_cov1 = banff_file1.Get("prefit_cov")     # TODO:test
    prefit_params1 = banff_file1.Get("prefit_params")     # TODO:test
    postfit_params1 = banff_file1.Get("postfit_params")     # TODO:test
    banff_untuned_central_values1 = banff_file1.Get("banff_untuned_central_values")     # TODO: test
    param_list_banff_xsec1 = banff_file1.Get("param_list_banff_xsec")     # TODO: test
    param_list_banff_xsec2 = banff_file2.Get("param_list_banff_xsec")
    banff_xsec_untuned_central_values1 = banff_file1.Get("banff_xsec_untuned_central_values")     # TODO: test

    nrows1 = banff1.GetNrows()
    nrows2 = banff2.GetNrows()
    ncols1 = banff1.GetNcols()
    ncols2 = banff2.GetNcols()

    reorder_dict = {}
    reorder_dict_xsec = {}
    num_new_params = nrows1+len(append_params)
    num_new_xsec_params = param_list_banff_xsec1.GetEntriesFast()
    for param in append_params:
        if append_params[param] == True:
            num_new_xsec_params += 1
    output_param_names = ROOT.TObjArray()
    output_postfit_cov = ROOT.TMatrixD(num_new_params, num_new_params)
    output_prefit_cov = ROOT.TMatrixD(num_new_params, num_new_params)
    output_prefit_params = ROOT.TVectorD(num_new_params)
    output_postfit_params = ROOT.TVectorD(num_new_params)
    output_param_list_banff_xsec = ROOT.TObjArray()
    output_banff_untuned_central_values = ROOT.TVectorD(num_new_params)
    output_banff_xsec_untuned_central_values = ROOT.TVectorD(num_new_xsec_params)

    if ncols1 != nrows1 or ncols2 != nrows2:
        print "Input matricies must be square! Exiting"
        util.Die(1)

    if nrows1 != nrows2:
        print "Different number of rows: nrows matrix 1 =", nrows1, "; nrows matrix 2 =", nrows2
        print "Will re-order the common parameters anyway, but the resulting matrix may or may not be sensible " \
              "(will be okay if this difference is due to parameters that only appear at the end of matrix 2 and " \
              "do not appear in matrix 1)"

    # Find the new ordering
    for iold in range(0,nrows1):
        param_name = param_names1.At(iold).String().Data()
        inew = find_index(param_names2,param_name)
        if inew == -1:
            print "param", param_name, "isn't present in the 2nd matrix. Exiting"
            util.Die(1)
        reorder_dict.update({inew:iold})
    for iold in range(0,param_list_banff_xsec1.GetEntriesFast()):
        param_name = param_list_banff_xsec1.At(iold).String().Data()
        inew = find_index(param_list_banff_xsec2,param_name)
        if inew == -1:
            print "param", param_name, "isn't present in the 2nd xsec matrix. Exiting"
            util.Die(1)
        reorder_dict_xsec.update({inew:iold})

    # Re-order existing parameters
    for inew in range(0,nrows1):
        output_param_names.Add(ROOT.TObjString(param_names1.At(reorder_dict[inew]).String().Data()))
        output_prefit_params[inew] = prefit_params1[reorder_dict[inew]]
        output_postfit_params[inew] = postfit_params1[reorder_dict[inew]]
        output_banff_untuned_central_values[inew] = banff_untuned_central_values1[reorder_dict[inew]]
        for jnew in range(0,nrows1):
            output_postfit_cov[inew][jnew] = banff1(reorder_dict[inew],reorder_dict[jnew])
            output_prefit_cov[inew][jnew] = prefit_cov1(reorder_dict[inew],reorder_dict[jnew])
    for inew in range(0,param_list_banff_xsec1.GetEntriesFast()):        
        output_param_list_banff_xsec.Add(ROOT.TObjString(param_list_banff_xsec1.At(reorder_dict_xsec[inew]).String().Data()))
        output_banff_xsec_untuned_central_values[inew] = banff_xsec_untuned_central_values1[reorder_dict_xsec[inew]]

    # Append new parameters
    for new_param in append_params:
        output_param_names.Add(ROOT.TObjString(new_param))
        if append_params[new_param] == True:
            output_param_list_banff_xsec.Add(ROOT.TObjString(new_param))        
    for inew in range(nrows1,num_new_params):
        output_prefit_params[inew] = 0
        output_postfit_params[inew] = 0
        output_banff_untuned_central_values[inew] = 0
        for jnew in range(nrows1,num_new_params):
            if inew == jnew:
                output_postfit_cov[inew][jnew] = 1      # TODO: Was previously set to 0. Is 1 okay? The new params need to be turned off, but the matrix must be positive-definite
                output_prefit_cov[inew][jnew] = 1       # TODO: Was previously set to 0. Is 1 okay? The new params need to be turned off, but the matrix must be positive-definite
            else:
                output_postfit_cov[inew][jnew] = 0
                output_prefit_cov[inew][jnew] = 0
    for inew in range(param_list_banff_xsec1.GetEntriesFast(),num_new_xsec_params):
        output_banff_xsec_untuned_central_values[inew] = 0    

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
    


def find_index(names,param):
    for i in range(0,names.GetSize()):
        if names.At(i).String().Data() == param or (param in rename and names.At(i).String().Data() == rename[param] and not skip_renamed):
            return i
    return -1



if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Re-order the first banff matrix according to the order of the parameters in the second.")
    parser.add_option("--input-1", dest="in1", default="", help='Input file containing a BANFF matrix')
    parser.add_option("--input-2", dest="in2", default="", help='Input file containing a different version of the BANFF matrix')
    parser.add_option("--output" , dest="out", default="", help='Output file name')

    v = optparse.Values()
    parser.parse_args(values = v)   # v contains only the values that were specified, not defaults
    o, other = parser.parse_args()

    main(o.in1, o.in2, o.out)

    if other:
        print "Warning - ignored arguments: ", other
