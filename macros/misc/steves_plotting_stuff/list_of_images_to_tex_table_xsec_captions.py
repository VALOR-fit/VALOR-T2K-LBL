#!/usr/bin/python
import sys
import glob
import ROOT


try:
  ROOT.valor
except:
  ROOT.gROOT.ProcessLine(".x $T2KOSCANA/root/root2k.C");
  
valor = ROOT.valor
util =  valor.util
conv =  valor.conv
opt  =  valor.Opt.Instance()
env  =  valor.Env.Instance()


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
\\usepackage{morefloats}

\setlength\parindent{0pt}

\\begin{document}

""")
  
  tex_output.write("\\begin{doublespacing}\n")
  
  curr_col = 0
  
  for image in images:
    
    label =  "fig:"+image.split("/")[-1].rsplit(".",1)[0]
    param_name = "_".join(label.split("_")[3:-1])
    param = util.PhysParamIdFromString(param_name)
    
    if curr_col == 0:
      tex_output.write("\n\\begin{figure*}[htb]")
    curr_col+=1
      
    tex_output.write("\n \\begin{minipage}[t]{"+figure_width+"}\n")
    tex_output.write("  \\centering\n  \\includegraphics[width=0.95\\textwidth]{"+image+"}\n")
    tex_output.write("  \\label{"+label+"}\n")
    tex_output.write("  \\caption{+1, +3, -1 and -3 $\sigma$ variations of the "+util.ParamAsStringPaperTable(param)+ " systematic.}\n")
    tex_output.write(" \\end{minipage}")
    
    if curr_col == ncols:
      tex_output.write("\n\\end{figure*}")
      curr_col = 0
    else:
      tex_output.write(" \\hspace{0.06\\linewidth}")
  
  if curr_col:
    tex_output.write("\n\\end{figure*}")
  
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
  
