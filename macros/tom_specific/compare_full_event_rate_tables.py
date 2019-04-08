#!/bin/bash python

from glob import glob


def make_table(binning, oscname):
    s = '\\begin{frame}\n' \
        '  \\frametitle{' + binning.replace('_',' ') + ' ' + oscname + '}\n' \
        '  \\scriptsize\n' \
        ' \\begin{center}\n' \
        '  \\begin{tabular}{|l|c|c|c|c|c|}\n' \
        '    \\hline\n' \
        '    & \\multicolumn{3}{c|}{FHC} & \\multicolumn{2}{c|}{RHC} \\\\\n' \
        '    & 1R$\\mu$ & 1R$e$ & 1R$e$+decay & 1R$\\mu$ & 1R$e$ \\\\\n' \
        '    \\hline\n'
    for fname in sorted(glob(binning + '/' + '*' + oscname + '*.tex')):
        print fname
        temp = fname.split('/')[-1].split('_',3)[-1].rsplit('_',1)[0]
        s += '    ' + temp.replace('_',' ')
        temp = []
        temposc = []
        with open(fname) as f:
            for l in f.readlines():
                if l.startswith('Total:'):
                    temp.append(float(l.split()[1]))
                elif l.startswith('\\item[]'):
                    x = l.split(None, 1)[1].rsplit(None,2)
                    temposc.append([x[0],x[-1]])
                elif l.startswith('Hierarchy:'):
                    temposc.append(l.split()[-1])
        for i in [0,1,4,2,3]:
            s += ' & %.3f' % temp[i]
        s += ' \\\\\n'
    s += '  \\hline\n' \
         '  \\end{tabular}\n\n' \
         '  \\vspace{1cm}\n' \
         '  \\begin{tabular}{l l l l}\n'
    for i in xrange(3):
        s += '    ' + temposc[i]  [0] + ' & ' + temposc[i]  [1]
        s += ' & '  + temposc[i+3][0] + ' & ' + temposc[i+3][1] + '\\\\\n'
    s += '    Hierarchy & ' + temposc[-1] + '\\\\\n'
    s += '  \\end{tabular}\n' \
         '  \\end{center}\n' \
         '\\end{frame}\n'
    return s

t = '\\documentclass[table,12pt]{beamer}\n' \
    '\\errorcontextlines 10000\n' \
    '\\usepackage[utf8]{inputenc}\n' \
    '\\usepackage{amsmath}\n' \
    '\\usetheme{Boadilla}\n' \
    '\\usecolortheme{whale}\n' \
    '\\date{\\today}\n' \
    '\\title{VALOR event rates}\n' \
    '\\author{Tom Dealtry}\n' \
    '\\institute{Lancaster University}\n' \
    '\\begin{document}\n'

for o in ['hkdrosc', 't2kasimovAosc', 'noosc']:
    for b in ['t2k_binning', 'osc3pp_binning', 'osc3ppextended_binning']:
        print b, o
        t += make_table(b,o)
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions}\n' \
     '  \\begin{center}\n' \
     '  \\begin{tabular}{l|l|l|l}\n' \
     '    \\hline\n' \
     '    & $E_{\\textrm{true}}$ & $\\mu$-like $E_{\\textrm{reco}}$ & $e$-like $E_{\\textrm{reco}}$ \\\\\n' \
     '    \\hline\n' \
     '    T2K & VALOR (84) & VALOR (73) & VALOR (25) \\\\\n' \
     '    osc3pp & Osc3++ (98) & VALOR (73) & VALOR (25) \\\\\n' \
     '    osc3ppextended & Osc3++ (98) & Osc3++X (87) & VALOR (25) \\\\\n' \
     '    \\hline\n' \
     '  \\end{tabular}\n' \
     '  \\end{center}\n' \
     '\\end{frame}\n'
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions: $E_{\\textrm{true}}$}\n' \
     '  VALOR 84 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 50-MeV bins from 0.0  to 0.3  GeV (6  bins)\n' \
     '  \\item  25-MeV bins from 0.3  to 1.0  GeV (28 bins)\n' \
     '  \\item  50-MeV bins from 1.0  to 3.0  GeV (40 bins)\n' \
     '  \\item 100-MeV bins from 3.0  to 3.5  GeV (5 bins)\n' \
     '  \\item 500-MeV bin from 3.5 to 4.0 GeV (1 bin)\n' \
     '  \\item 1000-MeV bin from 4.0 to 5.0 GeV (1 bin)\n' \
     '  \\item 2000-MeV bin from 5.0 to 7.0 GeV (1 bin)\n' \
     '  \\item 3000-MeV bin from 7.0 to 10.0 GeV (1 bin)\n' \
     '  \\item 20000-MeV bin from 10.0 to 30.0 GeV (1 bin)\n' \
     '  \\end{itemize}\n' \
     '\\end{frame}\n'
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions: $E_{\\textrm{true}}$}\n' \
     '  Osc3++ 98 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 50-MeV bins from 0.0  to 0.2  GeV (4  bins)\n' \
     '  \\item  25-MeV bins from 0.2  to 1.0  GeV (32 bins)\n' \
     '  \\item  50-MeV bins from 1.0  to 2.5  GeV (30 bins)\n' \
     '  \\item 100-MeV bins from 2.5  to 4.0  GeV (15 bins)\n' \
     '  \\item 200-MeV bins from 4.0  to 6.0  GeV (10 bins)\n' \
     '  \\item 500-MeV bins from 6.0  to 8.0  GeV (4 bins)\n' \
     '  \\item 1000-MeV bins from 8.0  to 10.0  GeV (2 bins)\n' \
     '  \\item 20000-MeV bin from 10.0 to 30.0 GeV (1 bin)\n' \
     '  \\end{itemize}\n' \
     '  \\footnotesize Note: This is an extension of the default Osc3++ binning, with extra bins at 0.00--0.05 and 10.0--30.0 GeV' \
     ' Requried due to spline binning\n' \
     '\\end{frame}\n'
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions: $E_{\\textrm{reco}}$}\n' \
     '  VALOR $\\nu_\\mu$ 73 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 50-MeV bins from 0.0  to 3.0  GeV (60  bins)\n' \
     '  \\item  250-MeV bins from 3.0  to 4.0  GeV (4 bins)\n' \
     '  \\item  500-MeV bins from 4.0  to 6.0  GeV (4 bins)\n' \
     '  \\item 1000-MeV bins from 6.0  to 10.0  GeV (4 bins)\n' \
     '  \\item 20000-MeV bin from 10.0 to 30.0 GeV (1 bin)\n' \
     '  \\end{itemize}\n' \
     '  VALOR $\\nu_e$ 25 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 50-MeV bins from 0.0  to 1.25  GeV (25  bins)\n' \
     '  \\end{itemize}\n' \
     '  VALOR $\\nu_e$ $\\theta_{\\textrm{reco}}$ 15 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 10-degree bins from 0  to 140 degrees (14  bins)\n' \
     '  \\item 40-degree bin from 140  to 180 degrees (1  bin)\n' \
     '  \\end{itemize}\n' \
     '\\end{frame}\n'
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions: $E_{\\textrm{reco}}$}\n' \
     '  Osc3++X $\\nu_\\mu$ 87 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 50-MeV bins from 0.0  to 3.0  GeV (60  bins)\n' \
     '  \\item  100-MeV bins from 3.0  to 4.0  GeV (10 bins)\n' \
     '  \\item  200-MeV bins from 4.0  to 6.0  GeV (10 bins)\n' \
     '  \\item  500-MeV bins from 6.0  to 8.0  GeV (4 bins)\n' \
     '  \\item 1000-MeV bins from 8.0  to 10.0  GeV (2 bins)\n' \
     '  \\item 20000-MeV bin from 10.0 to 30.0 GeV (1 bin)\n' \
     '  \\end{itemize}\n' \
     '  \\footnotesize Note: Osc3++ binning is by default different.' \
     ' VALOR can\'t use it without regenerating splines. This combines the Osc3++ binning with the spline binning\n' \
     '\\end{frame}\n'
t += '\\begin{frame}\n' \
     '  \\frametitle{Binning definitions: $E_{\\textrm{reco}}$}\n' \
     '  Osc3++ 56 bins\n' \
     '  \\begin{itemize}\n' \
     '  \\item 100-MeV bins from 0.0  to 4.0  GeV (40  bins)\n' \
     '  \\item 200-MeV bins from 4.0  to 6.0  GeV (10 bins)\n' \
     '  \\item 500-MeV bins from 6.0  to 8.0  GeV (4 bins)\n' \
     '  \\item 1000-MeV bins from 8.0  to 10.0  GeV (25 bins)\n' \
     '  \\end{itemize}\n' \
     '  \\footnotesize Note: This is the default Osc3++ binning for both $\\mu$-like and $e$-like.' \
     ' It\'s not used in the slides. Just here for reference\n' \
     '\\end{frame}\n'
t += '\\end{document}\n'
with open('osctables.tex','w') as f:
    f.write(t)
print t
