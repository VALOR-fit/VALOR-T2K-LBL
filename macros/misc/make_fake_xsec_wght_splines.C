//
// ---------------------------------------------------------------------------------------------------------
// Build "fake" cross-section weight splines to be used as analysis inputs for debugging/testing purposes.
//
// C.Andreopoulos
// ---------------------------------------------------------------------------------------------------------
//

#include <iostream>
#include <sstream>

using namespace std;

void make_fake_xsec_wght_splines()
{
  using namespace t2kana;
  using namespace t2kana::conv;

  const int n_knots = 100;

  const double twkd_min  = -5.0;
  const double twkd_max  = 5.0;
  const double twkd_step = (twkd_max - twkd_min) / (n_knots-1);

  double wght[n_knots];
  double twkd[n_knots];

  TFile f("splines.root","recreate");

  // loop over samples
  for(int is = 0; is < kns; is++) {
    SuperkSample_t s = kSKSamples[is];
    if(s == kSKsUndefined) continue;

    // loop over modes for current sample
    for(int im = 0; im < knm; im++) {
       SuperkMode_t m = kSKReactionModes[is][im];
       if(m == kSKmUndefined) continue;

       // loop over xsec systematic parameters
       for(int ip = 0; ip < kNumOfXSecSyst; ip++) {
          PhysParamId_t p = kListOfXSecSyst[ip];
          if(p == kPhypUndefined) continue;
          // splines expected by fitting code?
          bool need_splines = kXSecSystHasWghtSplines[ip];
          if(!need_splines) continue;
          // is this spline going to be applied to this sample/mode
          bool applied = kXsApplyXSecSystToSKMode[is][im][ip];
          if(!applied) continue;
          
          // loop over variable-size Ereco and Etrue bins for which xsec wght splines
          // are computed
          for(int ir = 0; ir < kXsNErecoBins; ir++) {
             double Ereco = kXsErecoBinEdges[ir] + kOneMeV;
             for(int it = 0; it < kXsNEtrueBins; it++) {
                double Etrue = kXsEtrueBinEdges[is][im][it] + kOneMeV;

                // Build spline name
                ostringstream splname;
                splname << "spl_"
                        << util::AsStringCompact(s) << "_"
                        << util::AsStringCompact(m) << "_"
                        << util::AsStringCompact(p) << "_"
                        << ir << "_"
                        << it;

                cout << " - Creating cross-section weight spline: " << splname.str() << endl;

                for(int i = 0; i < n_knots; i++) {
		   // twkd[i] is x-value of current point
		   // wght[i] is y-value of current point
                   twkd[i] = twkd_min + twkd_step * i;

                   if (Ereco > 0.400 && Ereco < 1.100 && m==kSKmNuMuCCQE  && p==kNuisBANFF22) {
                      cout << "Not null for QE" << endl;
                      wght[i] = (twkd[i] / twkd_max) + 1;
                   }
                   else
                   if (Ereco > 0.400 && Ereco < 1.100 && m==kSKmNuMuCC1pi && p==kNuisBANFF23) {
                      cout << "Not null for 1pi" << endl;
                      //wght[i] = 1./((twkd[i] / twkd_max) + 1 + conv::kASmallNum);
		      wght[i] = (twkd[i] / twkd_max) + 1;
                   }
                   else {
 		      wght[i] = (twkd[i] / twkd_max) + 1;
                   }
          
                }
                TSpline3 * spl = new TSpline3(splname.str().c_str(),twkd,wght,n_knots);
                f.cd();
                spl->Write(splname.str().c_str());

             }//true energy bins
          }//reco energy bins
       }//systematic params
    }//modes
  }//samples
  
  f.Close();
}
