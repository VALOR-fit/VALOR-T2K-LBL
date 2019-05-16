#!/usr/bin/env python

from datetime import date
import collections
import sys
#import numpy as np
from decimal import Decimal
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


#_______________________________________________________________________________________________________________________
# User defined variables
_slides_title = "$\\Delta m^2$ and $\\sin^{2}(\\theta_{23})$ PAC Marginalisation Stability Study"
_slides_title_small = "PAC Marginalisation Study"
_data_type  = collections.OrderedDict()
_allocation = collections.OrderedDict()
_hierarchy  = collections.OrderedDict()
_POT        = collections.OrderedDict()

# key = name in file, value = name in Tex file
_data_type["pure_AsimovA"] = "Asimov"
_data_type["hybrid"]       = "hybrid"

_allocation["FHC"]     = "all FHC"
_allocation["FHC_RHC"] = "equal"

_hierarchy["normal"]   = "NH"
_hierarchy["inverted"] = "IH"

_POT["7.8E21"] = "7.8"
_POT["10E21"]  = "10"
_POT["15E21"]  = "15"

#_______________________________________________________________________________________________________________________
def main(file_pattern, outfilename, nbatch, critval_str, var):

    crit_dchisq = get_crit_dchisq(critval_str)

    tex_output = open(outfilename,"w")
    tex_output.write("""\\documentclass{beamer}
\\usepackage{pgfpages}
\\usepackage{graphicx}
\\usepackage{amsmath}    %For easy SI units
\\usepackage{siunitx}    %For easy SI units
\\usepackage{booktabs}
\\usepackage{multirow}
\\usetheme{Madrid}
\\pgfpagesuselayout{resize to}[a4paper,landscape,border shrink=5mm]
\\beamertemplatenavigationsymbolsempty   % get rid of navigation bar
\\begin{document}\n\n
    """)

    tex_output.write("\\title[" + _slides_title_small + "] { \\begin{center} {" + _slides_title + "} \\end{center} }\n")
    tex_output.write("\\author[F. Bench] {Francis Bench, Maria Antonova}\n")
    tex_output.write("\\date{" + date.today().strftime("%B %d, %Y") + " \\ \\vspace{2cm}}\n")
    tex_output.write(""" \\titlegraphic{
  \\includegraphics[height=25px]{./images/logo/liverpool.png}
  \\hspace{3px}
  \\includegraphics[height=30px]{./images/logo/t2k_logo_large_white.png}
  \\hspace{3px}
  \\includegraphics[height=30px]{./images/logo/valor_nu_small_transparent.png}
}\n\n
\\begin{frame}[plain]
  \\titlepage
\\end{frame}\n\n
    """)

    for hier in _hierarchy:
        hier_tex = _hierarchy[hier]
        write_hier = True

        tex_output.write("\\begin{frame}{" + hier_tex + "; " + critval_str + " CL}\n")
        #tex_output.write("\t\\begin{tabular}{|c|c|c|c|c|c|} \\toprule\n")
        tex_output.write("\t\\begin{tabular}{|c|c|c|c|c|}\n\\hline\n")
        range_title = "Range"
        if var == "dmsq23":
            range_title = "$\\Delta m^2$ Range (1E-3)"
        elif var == "mixing23":
            range_title = "$\\sin^{2}(\\theta_{23})$ Range"
        #tex_output.write("\t\tMH & Data & Allocation & POT (1E21) & Mean " + range_title + " & Extreme " +
        #                 range_title + " \\\\ \n \\hline \n")
        tex_output.write("\t\tMH & Data & Allocation & POT (1E21) & Extreme " + range_title + " \\\\ \n \\hline \n")
        
        for type in _data_type:
            type_tex = _data_type[type]
            write_type = True
            
            for alloc in _allocation:
                alloc_tex = _allocation[alloc]
                write_alloc = True
                
                for pot in _POT:
                    pot_tex = _POT[pot]
                    
                    row_tex = "\t\t\t"
                    if write_hier:
                        nrows_hier = len(_data_type) + len(_allocation) + len(_POT)
                        row_tex += "\\multirow{" + str(nrows_hier) + "}{*}{" + hier_tex + "}"
                        write_hier = False
                    row_tex += " & "
                    
                    if write_type:
                        nrows_type = len(_allocation) + len(_POT)
                        row_tex += "\\multirow{" + str(nrows_type) + "}{*}{" + type_tex + "}"
                        write_type = False
                    row_tex += " & "

                    if write_alloc:
                        nrows_alloc = len(_POT)
                        row_tex += "\\multirow{" + str(nrows_alloc) + "}{*}{" + alloc_tex + "}"
                        write_alloc = False
                    row_tex += " & "

                    row_tex += pot_tex + " & "

                    range_low  = []
                    range_high = []
                    for ibatch in range(0, nbatch):
                        file_in = file_pattern.replace("[BATCH]"     , str(ibatch))
                        file_in = file_in.replace     ("[DATATYPE]"  , type)
                        file_in = file_in.replace     ("[ALLOCATION]", alloc)
                        file_in = file_in.replace     ("[MH]"        , hier)
                        file_in = file_in.replace     ("[POT]"       , pot)
                        file_in = file_in.replace     ("[VAR]"       , var)
                        print "Reading from file", file_in
                        f = ROOT.TFile(file_in, "read")
                        hchisq = f.Get("h_deltachisq")
                        
                        # Search from left:
                        ibin_high = -1
                        ibin_low  = -1
                        for ibin in range(1, hchisq.GetNbinsX()+1):
                            if hchisq.GetBinContent(ibin) > crit_dchisq:
                                ibin_high = ibin
                            elif hchisq.GetBinContent(ibin) <= crit_dchisq:
                                ibin_low = ibin
                                if ibin_high != -1:
                                    # Interpolate between bins
                                    range_low.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_high, ibin_low))
                                    break
                                else:
                                    print "ERROR: No bins > critval, when searching from left. Exiting."
                                    sys.exit()
                                    
                        # Search from right:
                        ibin_high = -1
                        ibin_low  = -1
                        for ibin in range(hchisq.GetNbinsX(), 0, -1):
                            if hchisq.GetBinContent(ibin) > crit_dchisq:
                                ibin_high = ibin
                            elif hchisq.GetBinContent(ibin) <= crit_dchisq:
                                ibin_low = ibin
                                if ibin_high != -1:
                                    # Interpolate between bins
                                    range_high.append(interpolate_crit_x(crit_dchisq, hchisq, ibin_low, ibin_high))
                                    break
                                else:
                                    print "ERROR: No bins > critval, when searching from right. Exiting."
                                    sys.exit()

                    # Print formatting
                    if var == "dmsq23":
                        range_low  = [1000.0 * i for i in range_low]
                        range_high = [1000.0 * i for i in range_high]


                    # Now that we have the x range for each fit, calculate the mean and error
                    # xrange_low_mean    = np.mean(range_low)
                    # xrange_low_stderr  = np.std(range_low, ddof=len(range_low)-1) / np.sqrt(len(range_low))
                    # xrange_high_mean   = np.mean(range_high)
                    # xrange_high_stderr = np.std(range_high, ddof=len(range_high)-1) / np.sqrt(len(range_high))
                    # mean_range_str = "[$" + "{0:.3g}".format(Decimal(xrange_low_mean)) + " \pm " + \
                    #                  "{0:.2g}".format(Decimal(xrange_low_stderr)) + ", " + \
                    #                  "{0:.3g}".format(Decimal(xrange_high_mean)) + \
                    #                  " \pm " + "{0:.2g}".format(Decimal(xrange_high_stderr))  + "$]"
                    extreme_range_str =  "[$" + "{0:.3g}".format(Decimal(min(range_low))) + ", " + \
                                     "{0:.3g}".format(max(range_high)) + "$]"

                    #row_tex += mean_range_str + " & " + extreme_range_str + "\\\\ \n"
                    row_tex += extreme_range_str + "\\\\ \n"
                    tex_output.write(row_tex)

                    # Lines between rows
                    cline_start = 1
                    cline_end   = 5
                    if pot != _POT.keys()[-1]:     # If we're not at the final iteration in the ordered dictionary
                        cline_start = 4
                    else:
                        if alloc != _allocation.keys()[-1]:
                            cline_start = 3
                        else:
                            cline_start = 2
                    tex_output.write("\t\t\t\\cline{" +str(cline_start) + "-" + str(cline_end) + "}\n")

        tex_output.write("\t\t\\hline\n\t\\end{tabular}\n")
        tex_output.write("\\end{frame}\n")
    tex_output.write("\\end{document}\n")

    print "Output written to ", outfilename


#_______________________________________________________________________________________________________________________
def interpolate_crit_x(crit_dchisq, hchisq, ibin1, ibin2):
    y1 = hchisq.GetBinContent(ibin1)
    y2 = hchisq.GetBinContent(ibin2)
    x1 = hchisq.GetBinLowEdge(ibin1)
    x2 = hchisq.GetBinLowEdge(ibin2) + hchisq.GetBinWidth(ibin2)
    gradient = (y2 - y1) / (x2 - x1)
    intercept = y1 - gradient * x1
    return (crit_dchisq - intercept) / gradient


#_______________________________________________________________________________________________________________________
def get_crit_dchisq(critval_str_in):
    # Assume the input can either be "Xsigma", "x%" or "x". Treat the latter as a decimal.
    critval = -1
    is_sigma = False
    if "sigma" in critval_str_in:
        critval = float(critval_str_in.replace("sigma", ""))
        is_sigma = True
    elif "%" in critval_str_in:
        critval = float(critval_str_in.replace("%", ""))/100.
    else:
        critval = float(critval_str_in)

    if is_sigma:
        return critval**2
    else:
        return sigma(critval)**2


#_______________________________________________________________________________________________________________________
def sigma(frac):
    return ROOT.TMath.Sqrt(2) * ROOT.TMath.ErfInverse(frac)


#_______________________________________________________________________________________________________________________
if __name__=="__main__":
  
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")
    parser.add_option('-p', '--file-pattern', dest='file_pattern', default="batch[BATCH]/fit_oa2018_[MH]_reactor_[VAR]_nexp_1_0.003_-999.000-job0_[ALLOCATION]_[POT]_[DATATYPE]-_0_0_contour.root", help="File naming pattern for input ROOT files.")
    parser.add_option('-o', '--out-file'    , dest='output_file' , default="PAC_conf_interval_comparison.tex"                                                                                      , help="Name for output TeX file.")
    parser.add_option('-n', '--n-batches'   , dest='nbatch'      , default=3,   type=int                                                                                                           , help="The number of batches of fits to compare.")
    parser.add_option('-c', '--crit-val'    , dest='critval'     , default="4sigma"                                                                                                                , help="Critical value (in sigma or percent) to evaluate the range at.")
    parser.add_option('-v', '--var'         , dest='var'         , default="dmsq23"                                                                                                                , help="Oscillation param used in fit. \'dmsq23\' or \'mixing23\'")

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o,other = parser.parse_args()

    main(o.file_pattern, o.output_file, o.nbatch, o.critval, o.var)

    if other:
        print "Warning - ignored arguments: ",other
