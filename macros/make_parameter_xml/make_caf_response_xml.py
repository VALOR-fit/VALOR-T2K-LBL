#!/usr/bin/env python

import ROOT
import re
import sys

parameter_regex = re.compile(r"(\w+)_cvwgt")

file_form = """<param_list name="from_caf">
[PARAMS]</param_list>
"""

param_form = """  <param>
    <name>
      <simple> [SNAME] </simple>
    </name>
    
    <type>
      response
    </type>
    
    <categories>
      <category> int </category>
      <category> caf </category>
    </categories>
      
    <apply_to>
      <detector>  all    </detector>
      <beam>      all    </beam>
      <selection> all    </selection>
      <reac_type> all    </reac_type>
      <flav>      all    </flav>
      <target>    all    </target>
    </apply_to>
  </param>

"""

def main(caf_filename):
  
  caf_file = ROOT.TFile(caf_filename)
  caf = caf_file.Get("caf")
  assert(caf)
  
  param_list = []
  
  
  for branch in caf.GetListOfBranches():
    branch_name = branch.GetName()
    
    result = parameter_regex.match(branch_name)
    
    if result is not None:
      # ~ print result.group()
      sname = result.group(1)
      # ~ print "sname: ",sname
      param_list.append(sname)
  # end branch list
  
  param_list_string = ""
  for param_sname in param_list:
    param_list_string += param_form.replace("[SNAME]",param_sname)
    print param_list_string
    # ~ print param_sname
    
  xml_output = open("caf_params.xml","w")
  
  # ~ print param_list_string
  
  xml_output.write(file_form.replace("[PARAMS]",param_list_string))
  
  xml_output.close()
  
if __name__ == "__main__":
  main(sys.argv[1])
  print "Fin!"
  
