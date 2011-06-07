//
// ------------------------------------------------------------------------------------------------
// Test Cholesky decomposition of SuperK CCQE ring-counting efficiency covariance matrix
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_cholesky_decomp_sk_ccqe_rc_covariance
(
  const int N               = 1000,  /* number of vectors generated with the Cholesky method */
  bool      save_cov_matrix = false, /* */
  bool      save_ntpl       = false  /* */
)
{
  using namespace std;
  using namespace t2kana;

  TMatrixD V(3,3);

  const double sig0 = 0.017; // Ereco: < 0.4 GeV
  const double sig1 = 0.035; // Ereco: 0.4 - 1.1 GeV
  const double sig2 = 0.093; // Ereco: > 1.1 GeV

  V(0,0) = sig0 * sig0;
  V(1,1) = sig1 * sig1;
  V(2,2) = sig2 * sig2;
  V(0,1) = 0.59 * sig0 * sig1;
  V(0,2) = 0.33 * sig0 * sig2;
  V(1,2) = 0.61 * sig1 * sig2;
  V(1,0) = V(0,1);
  V(2,0) = V(0,2);
  V(2,1) = V(1,2);

  cout << "Covariance matrix: " << endl;

  V.Print();
//V.Invert().Print();

  if(save_cov_matrix) {
     TFile f("./sk_ccqe_rc_err_covariance.root","recreate");
     V.Write("ccqe_rc_err_covariance");
     f.Close();
  }

  cout << "Running Cholesky decomposition: " << endl;

  TMatrixD L = util::CholeskyDecomposition(V);
  L.Print();

  TNtuple * ntfRC = new TNtuple("ntfRC","","twkRC0:twkRC1:twkRC2:wghtRC0:wghtRC1:wghtRC2");

  // mean tweaking dial
  TVectorD mean(3);
  mean(0)=0.;
  mean(1)=0.;
  mean(2)=0.;
 
  for(int i=0; i<N; i++) {
    // generate correlated tweaking dials
    TVectorD twkdial = util::CholeskyGenerateCorrelatedParams(L, mean);
    // print-out a few
    if(i%50==0) {
       cout << "CCQE RC efficiency tweaking vector " 
         << setw(6) << setfill(' ') << setprecision(3) << i << " : (";
       for(int ip=0; ip<3; ip++) {
          cout << setw(8) << setfill(' ') << setprecision(4) << twkdial(ip);
          if ( ip <  2) { cout << ", "; }
       }//ip
       cout << ")" << endl;
    }
    // calc weights
    TVector3 wght(3);
    wght(0) = 1 + twkdial(0);
    wght(1) = 1 + twkdial(1);
    wght(2) = 1 + twkdial(2);
    // store
    ntfRC->Fill(twkdial(0),twkdial(1),twkdial(2),wght(0),wght(1),wght(2));
  }//i

  if(save_ntpl) {
     TFile f("./sk_ccqe_rc_err_weights.root","recreate");
     ntfRC->Write();
     f.Close();
  }

}

