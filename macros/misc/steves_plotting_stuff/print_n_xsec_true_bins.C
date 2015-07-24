
using namespace valor;
using namespace valor::conv;
using namespace valor::util;

int print_n_xsec_true_bins()
{
  for (int is = 0 ; is < kns ; is++) {
    SuperkSample_t s = kSKSamples[is];
    if (s==kSKsUndefined) continue;
    
    for (int im = 0 ; im < knm ; im++) {
      SuperkMode_t m = kSKReactionModes[is][im];
      if(m == kSKmUndefined) continue;
      
      int nbins = 0;
      for (int ibin = 1 ; ibin <= kXsNEtrueBins ; ibin ++) {
        if (kXsEtrueBinEdges[is][im][ibin] < 0) break;
        nbins++;
      }
      
      cout<<AsString(s)<<"\t"<<AsString(m)<<"\t"<<nbins<<endl;
    }
  }
  
  return 0;
}
