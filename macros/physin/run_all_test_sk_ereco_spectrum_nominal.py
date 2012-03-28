#!/usr/bin/env python

# Runs test_sk_ereco_spectrum_nominal.C  with and without banff for 3f, 2f and no osc hypotheses
# Steve Dennis 28-03-12

import os
#~ import itertools

def main():
  macroname = '$T2KOSCANA/macros/physin/test_sk_ereco_spectrum_nominal.C'
  
  sample = "\\\"1Rmu\\\""
  theta_12 = 0.86
  theta_13 = 0.1
  theta_23 = 0.0
  delta_cp = 0.0
  dm12sq   = 7.58E-5
  dm13sq   = 2.32E-3

  hypotheses = ["\\\"3f_osc\\\"", "\\\"2f_numu_disap\\\"" ,"\\\"\\\"" ]
  boolstrs = ["true" , "false" ]
  #hypotheses = ["\\\"3f_osc\\\"" ]
  #boolstrs = ["true"]

  print sample

  #~ for hyp, banff in itertools.product(hypotheses, boolstrs):
  for hyp in hypotheses:
    for banff in boolstrs:
      
      command = "root2k -q -b \""+macroname+"("
      command += sample + ","

      command += hyp + ","

      command += `theta_12` + ","
      command += `theta_13` + ","
      command += `theta_23` + ","
      command += `delta_cp` + ","
      command += `dm12sq`   + ","
      command += `dm13sq`   + ","

      command += banff
      
      command += ")"
      command += "\""

      #print command

      if os.system(command):
        print "Something went wrong"
        return 1;
        
if __name__ == "__main__":
    main()
