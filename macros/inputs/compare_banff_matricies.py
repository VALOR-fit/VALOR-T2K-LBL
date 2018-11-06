#!/usr/bin/env python

# Macro to compare different versions of the BANFF matrix (does it have the same parameters in the same positions?
# Are their values the same? Same dimensions?)

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True

try:
    ROOT.valor
except:
    ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");

ROOT.gROOT.SetBatch(1)

valor = ROOT.valor
util  = valor.util
conv  = valor.conv


# If a parameter has been renamed between BANFF matrix versions, put the name as it appears in the 1st input matrix as
# the key and the name from the 2nd input matrix as the value in the following dictionary
rename = {
    "CC_nu"           : "CC_norm_nu"   ,
    "CC_nubar"        : "CC_norm_nubar",
    "SKFSI_Absorption": "FEFABS"       ,
    "SKFSI_QEl"       : "FEFQE"        ,
    "SKFSI_CXl"       : "FEFCX"        ,
    "SKFSI_QEh"       : "FEFQEH"       ,
    "SKFSI_Inel"      : "FEFINEL"
}

# Use this option to not count the renamed parameters when counting the number of mismatched elements
# (e.g. In one matrix I have the renamed parameters switched off in one version, but not the other and I only want to
# search for discrepancies in the other parameters)
skip_renamed = True



def main(infilename1, infilename2):
    identical = True
    banff_file1 = ROOT.TFile(infilename1,"read")
    banff_file2 = ROOT.TFile(infilename2,"read")
    banff1 = banff_file1.Get("postfit_cov")
    banff2 = banff_file2.Get("postfit_cov")
    param_names1 = banff_file1.Get("param_list")
    param_names2 = banff_file2.Get("param_list")

    nrows1 = banff1.GetNrows()
    nrows2 = banff2.GetNrows()
    ncols1 = banff1.GetNcols()
    ncols2 = banff2.GetNcols()

    if ncols1 != nrows1 or ncols2 != nrows2:
        print "Input matricies must be square!"
        util.Die(1)

    nrows       = 0
    matrix1     = 0
    matrix2     = 0
    names1      = 0
    names2      = 0
    matrix1_tag = ""
    matrix2_tag = ""
    if nrows1 != nrows2:
        print "Different number of rows: nrows matrix 1 =", nrows1, "; nrows matrix 2 =", nrows2
        print "Will compare the common parameters anyway"
        identical = False
        if nrows1 > nrows2:
            nrows   = nrows1
            matrix1 = banff1
            matrix2 = banff2
            names1  = param_names1
            names2  = param_names2
            matrix1_tag = "BANFF1"
            matrix2_tag = "BANFF2"
        else:
            nrows = nrows2
            matrix1 = banff2
            matrix2 = banff1
            names1  = param_names2
            names2  = param_names1
            matrix1_tag = "BANFF2"
            matrix2_tag = "BANFF1"
            # TODO: construct reverse dictionary too
    else:
        nrows = nrows1
        matrix1 = banff1
        matrix2 = banff2
        names1  = param_names1
        names2  = param_names2
        matrix1_tag = "BANFF1"
        matrix2_tag = "BANFF2"

    # Check for mismatches between covariance values between the two input BANFF matrices
    # Also check for params which appear in matrix1 but not in matrix2
    checked_params = []
    params_in_matrix1_not_in_matrix2 = []
    n_elems_mismatched = 0
    params_moved = False
    for row1 in range(0,nrows):
        # The parameter for this element of matrix1 may not be in the same position in matrix2, so search for this parameter in matrix2
        param_row = names1.At(row1).String().Data()
        if param_row in rename and skip_renamed:
            continue    # TODO: reversed map as well (in case the matrices have been swapped round)
        checked_params.append(param_row)
        row2 = find_index(names2,param_row)
        if row2 == -1:
            params_in_matrix1_not_in_matrix2.append(param_row)
            identical = False
            continue
        if row1 != row2 and not params_moved:
            params_moved = True
        for col1 in range(0,nrows):
            param_col = names1.At(col1).String().Data()
            col2 = find_index(names2,param_col)
            if col2 == -1:
                identical = False
                continue
            if matrix1(row1,col1) != matrix2(row2,col2):
                n_elems_mismatched += 1
                identical = False
                print matrix1_tag, "-", matrix2_tag, "for (", param_row, ",", param_col, ") =", matrix1(row1,col1) - matrix2(row2,col2)

    if n_elems_mismatched == 0:
        print "No elements are mismatched"
    else:
        print n_elems_mismatched, "elements are different between the two matrices"

    print "Parameters that appear in", matrix1_tag, "but not", matrix2_tag, ":"
    for param in params_in_matrix1_not_in_matrix2:
        print " ", param

    print "Parameters that appear in", matrix2_tag, "but not", matrix1_tag, ":"
    for i in range(0,names2.GetSize()):
        if names2.At(i).String().Data() not in checked_params:
            match = False
            for rname in rename:
                if rename[rname] == names2.At(i).String().Data():
                    match = True
                    break
            if not match:
                print " ", names2.At(i).String().Data()
                identical = False

    if identical:
        str1 = "The content of the input matrices are identical"
        str2 = ""
        if params_moved:
            str2 = ", other than the positions of some parameters having been moved"
        print str1 + str2



def find_index(names,param):
    for i in range(0,names.GetSize()):
        if names.At(i).String().Data() == param or (param in rename and names.At(i).String().Data() == rename[param] and not skip_renamed):
            return i
    return -1



if __name__=="__main__":
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..", description="Plot spectra")
    parser.add_option("--input-1", dest="in1", default="", help='Input file containing a BANFF matrix')
    parser.add_option("--input-2", dest="in2", default="", help='Input file containing a different version of the BANFF matrix')

    v = optparse.Values()
    parser.parse_args(values = v)   # v contains only the values that were specified, not defaults
    o, other = parser.parse_args()

    main(o.in1, o.in2)

    if other:
        print "Warning - ignored arguments: ", other
