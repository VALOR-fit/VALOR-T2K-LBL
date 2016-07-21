#!/usr/bin/env python
import sys
import glob


def main(images, outfile, ncols):
  
  figure_width  = str(0.99/float(ncols))+"\\textwidth"
  
  # Start writing tex
  
  tex_output = open(outfile,"w")
  
  tex_output.write("""\documentclass[10pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage{amsmath}
\\usepackage{amsfonts}
\\usepackage{amssymb}
\\usepackage{fullpage}
\\usepackage{graphicx}
\\usepackage{float}
\\usepackage{setspace}

\setlength\parindent{0pt}

\\begin{document}

""")
  
  tex_output.write("\\begin{doublespacing}\n")
  
  for image in images:
    
    label =  "fig:"+image.split("/")[-1].rsplit(".",1)[0]
    tex_output.write("\n\\begin{minipage}[t]{"+figure_width+"}\n")
    tex_output.write(" \\centering\n \\includegraphics[width=0.95\\textwidth]{"+image+"}\n")
    tex_output.write(" \\label{"+label+"}\n")
    tex_output.write("\\end{minipage}")
  
  tex_output.write("\n\\end{doublespacing}\n")
  
  tex_output.write("\n\n\\end{document}\n")
  
  #~ print "LaTeX written to ",outfile
  

if __name__=="__main__":
  import optparse
  parser = optparse.OptionParser(usage="%prog [options]..",description="Make contours")
  
  parser.add_option('-x' , "--n-plots-per-row",dest='ncols' , type=int, default=2, help="Number of plots per row")
  parser.add_option('-f','--output-file' , dest='outfile', default="figures.tex", help="Output TeX file name")
  
  o,images = parser.parse_args()
  
  
  
  main(images, o.outfile,o.ncols)
  
