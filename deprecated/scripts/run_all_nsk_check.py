#!/usr/bin/env python

# Runs nsk_check with and without banff for 3f, 2f and no osc hypotheses
# Steve Dennis 27-03-12

import os
#~ import itertools

def main():
  macroname = '$T2KOSCANA/macros/physin/nsk_check.C'

  hypotheses = ['3f_osc', '2f_numu_disap', '' ]
  boolstrs = ['true' , 'false' ]

  #~ for hyp, banff in itertools.product(hypotheses, boolstrs):
  for hyp in hypotheses:
    for banff in boolstrs:
      if os.system("root2k -q -b '"+macroname+"("+banff+",\""+hyp+"\")'"):
        print "Something went wrong"
        return 1;
        
if __name__ == "__main__":
    main()
