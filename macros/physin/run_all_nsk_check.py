#!/usr/bin/env python

# Runs nsk_check with and without banff for 3f, 2f and no_osc hypotheses
# Steve Dennis 27-03-12

import os
import itertools

macroname = '$T2KOSCANA/macros/physin/nsk_check.C'

hypotheses = ['t2kana::conv::kHypOsc3F', 't2kana::conv::kHypNuMuDisapp2F' , 't2kana::conv::kHypUndefined' ]
boolstrs = ['true' , 'false' ]

for hyp, banff in itertools.product(hypotheses, boolstrs):
  if os.system("root2k -q -b '"+macroname+"("+banff+","+hyp+")'"):
    print "Something went wrong"
    break
  
  
  break
