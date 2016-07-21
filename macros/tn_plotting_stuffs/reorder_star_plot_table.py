#!/usr/bin/env python

# Sorts a tex file containing a star plot table
# according to it's column's numerical values
# Not very robust
# Does not preserve order for lines with identical values

import sys

def main(infile):
  list_of_columns_to_sort_by = [1,2]
  descending_order = True

  f = open(infile)

  lines = [l for l in f]
  f.close()

  for column in list_of_columns_to_sort_by:
    sortable = []
    for line in lines:
      split_line = line.split("     &")
      try:
        sort_by = float(split_line[column].rstrip("\\ \n"))
      except:
        continue
      
      sortable.append( (sort_by,line) )
    
    for j,line in sorted(sortable, reverse=descending_order):
      print line.rstrip()
      
    print "\n\n_________\n\n"

if __name__=="__main__":
  main(sys.argv[1])
