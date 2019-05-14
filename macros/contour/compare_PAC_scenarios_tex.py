#!/usr/bin/env python

from datetime import date
import collections
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True


#_______________________________________________________________________________________________________________________
# User defined variables
_slides_title = "$\\Delta m^2$ and $\\sin^{2}(\\theta_{23})$ PAC Marginalisation Stability Study"
_slides_title_small = "PAC Marginalisation Study"
_vars       = collections.OrderedDict()   # One section per variable
_data_type  = collections.OrderedDict()   # -|
_allocation = collections.OrderedDict()   #  | A single page for each combination of values within these three lists
_hierarchy  = collections.OrderedDict()   # -|
_POT        = collections.OrderedDict()   # Along the columns
_ntoys      = ["40k", "100k", "300k"]     # Along the rows

# key = name in file, value = name in Tex file
_vars["dmsq23"]  = "$\\Delta m_{23}^{2}$"
_vars["theta23"] = "$\\sin^{2}(\\theta_{23})$"

_data_type["pure_AsimovA"] = "Pure Asimov A"
_data_type["hybrid"]       = "hybrid"

_allocation["FHC"]     = "100\% FHC"
_allocation["FHC_RHC"] = "50/50 FHC/RHC"

_hierarchy["normal"]   = "NH"
_hierarchy["inverted"] = "IH"

_POT["7.8E21"] = "$7.8 \\times 10^{21}$"
_POT["10E21"]  = "$10 \\times 10^{21}$"
_POT["15E21"]  = "$15 \\times 10^{21}$"

# If a combination of a variable an ntoys appears here, then skip it.
_skip_var_ntoys = dict()
_skip_var_ntoys["theta23"] = ["100k", "300k"]

#_______________________________________________________________________________________________________________________
def main(file_pattern, outfilename):
    ncols = len(_POT)
    nrows = len(_ntoys)

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

    for var in _vars:
        # Write a page showing only the name of the variable that the following contours correspond to
        var_tex = _vars[var]
        tex_output.write("\\begin{frame}{} \n\t \\begin{center} \n\t\t \\huge{" + var_tex + "} \n\t \\end{center} \n\\end{frame}\n\n")

        for type in _data_type:
            type_tex = _data_type[type]

            for alloc in _allocation:
                alloc_tex = _allocation[alloc]

                for hier in _hierarchy:
                    hier_tex = _hierarchy[hier]
                    slide_title = var_tex + "; " + type_tex + "; Future POT " + alloc_tex + "; " + hier_tex
                    # Make a new page
                    tex_output.write("\\begin{frame}{" + slide_title + "}\n")
                    tex_output.write("\t\\begin{tabular}{c|" + "c"*(ncols) + "} \\toprule\n")

                    # Make the table title rows
                    table_title_row1 = "\t\t\\multirow{2}{*}{\\# Toys} & \\multicolumn{" + str(ncols) + "}{c}{Exposure (POT)}\\\\ \n"
                    table_title_row2 = "\t\t"
                    for pot in _POT:
                        pot_tex = _POT[pot]
                        table_title_row2 += " & " + pot_tex
                    table_title_row2 += "\\\\ \n"
                    tex_output.write(table_title_row1)
                    tex_output.write("\t\t\\cline{2-" + str(ncols+1) + "}\n")
                    tex_output.write(table_title_row2)

                    # Write each row of graphs
                    for n in _ntoys:
                        if var in _skip_var_ntoys:
                            if n in _skip_var_ntoys[var]:
                                continue
                        table_row_tex = "\t\t\t" + n
                        for pot in _POT:
                            pdf_file_in = file_pattern.replace("[VARS]"      , var)
                            pdf_file_in = pdf_file_in.replace ("[DATATYPE]"  , type)
                            pdf_file_in = pdf_file_in.replace ("[ALLOCATION]", alloc)
                            pdf_file_in = pdf_file_in.replace ("[MH]"        , hier)
                            pdf_file_in = pdf_file_in.replace ("[NTOYS]"     , n)
                            pdf_file_in = pdf_file_in.replace ("[POT]"       , pot)
                            pdf_file_in = pdf_file_in.replace (".pdf"        , "")
                            pdf_file_in = "{" + pdf_file_in + "}.pdf"   # This is to stop latex throwing an error when there are dots in the file name
                            fig_width = ( 1.0/float(ncols) ) - 0.05
                            table_row_tex += "\n & \\includegraphics[width={:.2f}".format(fig_width) + "\\textwidth]{" + pdf_file_in + "}"
                        table_row_tex += "\\\\\n"
                        tex_output.write(table_row_tex)
                    tex_output.write("\t\\bottomrule\n\\end{tabular}\n")
                    tex_output.write("\\end{frame}\n\n")
    tex_output.write("\\end{document}\n")

    print "Output written to ", outfilename


#_______________________________________________________________________________________________________________________
if __name__=="__main__":
  
    import optparse
    parser = optparse.OptionParser(usage="%prog [options]..",description="Make tables of validation spectra")

    parser.add_option('-p', '--file-pattern', dest='file_pattern', default="images/T2K/PAC_Apr2019/[VARS]/[NTOYS]_nominal_range/contour_[VARS]_[MH]_[ALLOCATION]_[POT]_[DATATYPE].pdf", help="File naming pattern for input pdf files.")
    parser.add_option('-o', '--out-file'    , dest='output_file' , default="PAC_contours_comparison.tex"                                                                              , help="Name for output TeX file.")

    v = optparse.Values()
    parser.parse_args(values=v) # v contains only the values that were specified, not defaults

    o,other = parser.parse_args()

    main(o.file_pattern, o.output_file)

    if other:
        print "Warning - ignored arguments: ",other
