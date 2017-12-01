#!/usr/bin/env python
import argparse

samples =["1Rmu",
          "1Re",
          "RHC1Rmu",
          "RHC1Re",
          "MultiRe",          
]

modes=["CCQE",
       "CC1$\pi$",
       "CC coherent",
       "CCn$\pi$",
       "CC other",
       "NC$\pi^{0}$",
       "NC$\pi^{+/-}$",
       "NC coherent",
       "NC other",
       "2p-2h",
       "NC 1$\gamma$"]

sheet_order = ["CCQE",
               "CC1$\pi$",
               "CCcoherent",
               "2p-2h",
               "CCother",
               "NC1pi",
               "NCcoherent",
               "NCother",
               "NC1$\gamma$"]

sum_modes = { "CCQE"         :["CCQE"        ],
              "CC1$\pi$"     :["CC1$\pi$"    ],
              "CCcoherent"  :["CC coherent" ],
              "2p-2h"        :["2p-2h"       ],
              "CCother"     :["CC other"    ],
              "NC1pi"       :["NC 1pi"      ],
              "NCcoherent"  :["NC coherent" ],
              "NCother"     :["NC other"    ],
              "NC1$\gamma$" :["NC 1$\gamma$"],
              "NC1pi"       :["NC$pi^{0}$", "NC$pi^{+/-}$"],
              "CCother"     :["CCn$pi$", "CC other"]} 


def main(filename):

  list_mode(filename)
    
  return

def list_mode(filename):
  sum_dict = {}
  raw_dict = {}
  nu_mu    = {}
  nu_e     = {}
  
  nu_mu_bar = {}
  nu_e_bar  = {}

  osc_nu_e     = {}
  osc_nu_e_bar = {}

  total = {}
  
  sample_str = ""
  sheet_sum = {}
  with open(filename) as file:
    for line in file:
      if "\begin{tabular}" in line:
        sum_dict = {}

      for sample in samples:
        if sample in line:
          #sample_str = sample.replace("\\","")
          print ""
          print sample_str
          print "-------o-------"

      for mode in modes:
        if mode in line:
          entries = [x.strip() for x in (line.replace("\\","")).split('&')]

          raw_dict[entries[0].replace("\\","")] = tuple(map(float, entries[1:]))
          nu_mu[entries[0].replace("\\","")] = float(entries[1])
          nu_e [entries[0].replace("\\","")] = float(entries[2])

          nu_mu_bar[entries[0].replace("\\","")] = float(entries[3])
          nu_e_bar [entries[0].replace("\\","")] = float(entries[4])

          osc_nu_e     [entries[0].replace("\\","")] = float(entries[5])
          osc_nu_e_bar [entries[0].replace("\\","")] = float(entries[6])

          total[entries[0].replace("\\","")] = float(entries[1]) + float(entries[2]) + float(entries[3]) + float(entries[4]) + float(entries[5]) + float(entries[6])

          
          sum_dict = dict((k, sum(i for i in v)) for k, v in raw_dict.items())

      if "\end{tabular}" in line:
       
        #sheet_sum = {k: sum(sum_dict[i.replace("\\","")] for i in v if isinstance(sum_dict[i.replace("\\","")], float)) for k, v in sum_modes.items()}
        
        # for key in modes:
        #   for keys in sum_dict:
        #     if key.replace("\\","") == keys:
        #       print ("{0:14s} {1:.4f}".format(keys, sum_dict[keys]))

        # print "\nComparison Sheet Order\n"

        # for key in sheet_order:
        #   for s_key in sum_modes:
        #     if key == s_key:
        #       mode_total = 0
        #       for values in sum_modes[s_key]:
        #         t_key = values.replace("\\","")
        #         #print "------------------>{}".format(t_key)                
        #         if t_key in sum_dict:
        #           mode_total += sum_dict[t_key]
        #           #print "The value {0:20s}: {1:.4f}".format(values, sum_dict[t_key])
        #       print "{0:14s} {1:.4f}".format(key, mode_total)

        print "\nComparison Sheet Order Nu_mu\n"


        sum_print(total)
        # sum_print(nu_mu)
        # sum_print(nu_e)
        # sum_print(nu_mu_bar)
        # sum_print(nu_e_bar)
        # sum_print(osc_nu_e)
        # sum_print(osc_nu_e_bar)
        print ""
        print "{0:14s} {1:.4f}".format("nu_mu"        ,sum(nu_mu.itervalues()))
        print "{0:14s} {1:.4f}".format("nu_e"         ,sum(nu_e.itervalues()))
        print "{0:14s} {1:.4f}".format("nu_mu_bar"    ,sum(nu_mu_bar.itervalues()))
        print "{0:14s} {1:.4f}".format("nu_e_bar"     ,sum(nu_e_bar.itervalues()))        
        print "{0:14s} {1:.4f}".format("osc_nu_e"     ,sum(osc_nu_e.itervalues()))
        print "{0:14s} {1:.4f}".format("osc_nu_e_bar" ,sum(osc_nu_e_bar.itervalues()))        
              
        # print sum_dict
      
  return sum_dict

def sum_print(dictionary):
  for key in sheet_order:
    for s_key in sum_modes:
      if key == s_key:
        mode_total = 0
        for values in sum_modes[s_key]:
          t_key = values.replace("\\","")
          #print "------------------>{}".format(t_key)                
          if t_key in dictionary:
            mode_total += float(dictionary[t_key])
            #print "The value {0:20s}: {1:.4f}".format(values, dictionary[t_key])
        print "{0:14s} {1:.4f}".format(key, mode_total)


if __name__=="__main__":
  parser = argparse.ArgumentParser(description="Get input file")
  parser.add_argument('-f', '--latex-file', dest='filename', help="LaTeX file containing event rates")
  
  args = parser.parse_args()
  
  main(args.filename)
