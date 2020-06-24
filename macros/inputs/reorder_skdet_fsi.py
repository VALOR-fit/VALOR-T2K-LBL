#!/usr/bin/env python

# Re-order the input SKDetFSI matrix, using a hard-coded order.
#
# Andy Chappell
# Francis Bench
# Sept. 2019 (run 1-9d analysis compatible)

import ROOT
  
ROOT.gROOT.SetBatch(1)
#______________________________________________________________________________________________________________________________

def main(skdet_file, output_file="skdetsipn.root"):
  print "Reading", skdet_file 
  skdet_tf  = ROOT.TFile(skdet_file)

  sk_matrix = skdet_tf.Get("SKJointError_Erec_Total")

  n = sk_matrix.GetNrows()
  for i in range(0,n):
    print "%02d:  %.8f" % (i, sk_matrix(i,i))
    

  out_tf = ROOT.TFile(output_file,"RECREATE")

  reordered_range = range(12,18)+range(0,12)+range(30,36)+range(18,30)+range(36,44)+range(44,45)
  sk_reordered = reorder_tmatrix(sk_matrix, reordered_range)

  sk_reordered.Write("skdetsipn")
  

  print "reordered matrix"
  for i in range(0,n):
    print "%02d: %.8f" % (i, sk_reordered[i][i])



  
  out_tf.Close()


    
def reorder_tmatrix(matrix,new_order):
  n = matrix.GetNrows()
  assert(n==matrix.GetNcols())
  assert(n==len(new_order))
  new_matrix = ROOT.TMatrixD(n,n)
  for inew,iold  in enumerate(new_order):
    for jnew,jold in enumerate(new_order):
      new_matrix[inew][jnew] = matrix[iold][jold]
  return new_matrix


if __name__=="__main__":
  
  import optparse
  parser = optparse.OptionParser(usage="%prog [options] skdet.root ",description="reorder SK matrix for VALOR")
  
  parser.add_option('-o', '--output_file', dest='output_file', default="skdetsipn.root", help="Current matrix")
  
  o,other = parser.parse_args()
  
  main(other[0],o.output_file)
  
  if other[3:]:
    print "Warning - ignored arguments: ",other[3:]
