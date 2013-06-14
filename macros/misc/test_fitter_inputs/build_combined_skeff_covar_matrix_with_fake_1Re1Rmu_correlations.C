//
// ----------------------------------------------------------------------------------------------------
// Builds SK efficiency covariance matrices assuming correlations between 1Rmu and 1Re bins.
//
// Covariance matrix is 19x19
//
// The matrix is organized in the following 7 blocks:
// ----------------------------------------------
// Block id  | Comment         | Matrix elements
// ----------------------------------------------
//    a      | 1Rmu numu CC    |   0-4
//    b      | 1Rmu NC         |   5
//    c      | 1Rmu nue CC     |   6
// ...............................................
//    d      | 1Re  osc nue CC |   7-9
//    e      | 1Re  numu CC    |   10-12
//    f      | 1Re  nue CC     |   13-15
//    g      | 1Re  NC         |   16-18
// ----------------------------------------------
//
// The C_{1Rmu} (7x7) and C_{1Re} (12x12) blocks of the combined 19x19 matrix are the default ones.
// This macro works out the C_{1Rmu;1Re} (7x12) block that describes the 1Rmu-1Re correlations.
// 
// C_{1Rmu;1Re} is set to C^{klm} where k,l,m = +,- or 0.
//
// + stands-for highly correlated
// - stands-for highly anti-correlated
// 0 stands-for not correlated
// 
// Symbol k describes the correlation between the following blocks,
// describing similar modes in 1Rmu and 1Re samples:
//  >  a - e
//  >  b - g
//  >  c - d,f
// Symbol l describes the correlations between numu and nue CC across the 2 samples:
//  >  a - d,f
//  >  c - e
// Symbol m describes the correlations between CC and NC across the 2 samples:
//  > a,c - g
//  > b - d,e,f
//
// So C^{+-0} is a matrix where:
//  > (1Rmu numu CC - 1Re numu CC), 
//    (1Rmu NC      - 1Re NC), 
//    (1 Rmu nue CC - 1Re nue CC / osc nue CC) are highly correlated
//  > (1Rmu numu CC - 1Re nue CC / osc nue CC),
//    (1Re nue CC   - 1Re numu CC) are highly anti-correlated
//  > other correlations are ignored 
//
// So C^{--+} is a matrix where:
//  > (1Rmu numu CC - 1Re numu CC), 
//    (1Rmu NC      - 1Re NC), 
//    (1 Rmu nue CC - 1Re nue CC / osc nue CC) are highly anti-correlated
//  > (1Rmu numu CC - 1Re nue CC / osc nue CC),
//    (1Re nue CC   - 1Re numu CC) are highly anti-correlated
//  > (1Rmu numu/nue CC - 1Re NC),
//  > (1Rmu NC - 1Re numu/nue/osc nue CC) are highly correlated.
//
//
// When a + or - index is set for k,l,m the maximum amount of correlation or anti-correlation
// is set and it is checked whether the combined covariance matrix is invertible.
// If not, the amount of correlation or anti-correlation is relaxed till the matrix is invertible.
//
// Will generate matrices for the following 14 scenarios:
//   - C+00, C-00
//   -    C++0, C-+0
//   -       C+++, C-++
//   -       C++-, C-+-
//   -    C+-0, C--0
//   -       C+-+, C--+
//   -       C+--, C---
//
// C.Andreopoulos
// ----------------------------------------------------------------------------------------------------
//

#include <sstream>
#include <iomanip>

using namespace std;

void print_matrix(TMatrixD * m, int w, int prec);

void build_combined_skeff_covar_matrix_with_fake_1Re1Rmu_correlations(int k=0, int l=0, int m=0)
{
  gSystem->Load("libMatrix");

  //
  // Constants
  //

  const int n     = 19;  // num of elements in combined matrix
  const int n1Rmu =  7;  // num of elements in 1Rmu block
  const int n1Re  = 12;  // num of elements in 1Re block

  const int ngroups = 7; // num of groups identified above (a-g)

  // first/last index (inclusive) of each group in the combined covariance matrix
  const int idx_range[ngroups][2] = {
     /* a  */  { 0,  4}, 
     /* b  */  { 5,  5},
     /* c  */  { 6,  6},
     /* d  */  { 7,  9},
     /* e  */  {10, 12},
     /* f  */  {13, 15},
     /* g  */  {16, 18}
  };

  // the k inded controls correlations/anti-correlations between pairs of
  // groups as indicated below (pairs marked with 1)
  const int kcov[ngroups][ngroups] = {
               /*   a,   b,   c,   d,   e,   f,   g   */
     /* a  */  {    0,   0,   0,   0,   1,   0,   0   }, 
     /* b  */  {    0,   0,   0,   0,   0,   0,   1   },
     /* c  */  {    0,   0,   0,   1,   0,   1,   0   },
     /* d  */  {    0,   0,   1,   0,   0,   0,   0   },
     /* e  */  {    1,   0,   0,   0,   0,   0,   0   },
     /* f  */  {    0,   0,   1,   0,   0,   0,   0   },
     /* g  */  {    0,   1,   0,   0,   0,   0,   0   }
  };

  // the l index controls correlations/anti-correlations between pairs of
  // groups as indicated below (pairs marked with 1)
  const int lcov[ngroups][ngroups] = {
               /*   a,   b,   c,   d,   e,   f,   g   */
     /* a  */  {    0,   0,   0,   1,   0,   1,   0   }, 
     /* b  */  {    0,   0,   0,   0,   0,   0,   0   },
     /* c  */  {    0,   0,   0,   0,   1,   0,   0   },
     /* d  */  {    1,   0,   0,   0,   0,   0,   0   },
     /* e  */  {    0,   0,   1,   0,   0,   0,   0   },
     /* f  */  {    1,   0,   0,   0,   0,   0,   0   },
     /* g  */  {    0,   0,   0,   0,   0,   0,   0   }
  };

  // the m index controls correlations/anti-correlations between pairs of
  // groups as indicated below (pairs marked with 1)
  const int mcov[ngroups][ngroups] = {
               /*   a,   b,   c,   d,   e,   f,   g   */
     /* a  */  {    0,   0,   0,   0,   0,   0,   1   }, 
     /* b  */  {    0,   0,   0,   1,   1,   1,   0   },
     /* c  */  {    0,   0,   0,   0,   0,   0,   1   },
     /* d  */  {    0,   1,   0,   0,   0,   0,   0   },
     /* e  */  {    0,   1,   0,   0,   0,   0,   0   },
     /* f  */  {    0,   1,   0,   0,   0,   0,   0   },
     /* g  */  {    1,   0,   1,   0,   0,   0,   0   }
  };
   
  // Instantiante covar matrix
  TMatrixD * pcov = new TMatrixD(n,n);
  TMatrixD & cov = *pcov;

  // Reset all elements to 0
  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      cov(i,j) = 0.;
    }
  }

  //
  // Add standard 1Rmu, 1Re parts
  //

  // 1Rmu 
  cov(0,0) = TMath::Power(  1.7/100., 2.);
  cov(1,1) = TMath::Power(  3.5/100., 2.);
  cov(2,2) = TMath::Power(  9.3/100., 2.);
  cov(3,3) = TMath::Power(  1.8/100., 2.);
  cov(4,4) = TMath::Power( 20.0/100., 2.);
  cov(5,5) = TMath::Power(111.0/100., 2.);
  cov(6,6) = TMath::Power(100.0/100., 2.);
  cov(0,1) = 0.59 * ( 1.7/100.) * ( 3.5/100.);
  cov(0,2) = 0.33 * ( 1.7/100.) * ( 9.3/100.);
  cov(1,2) = 0.61 * ( 3.5/100.) * ( 9.3/100.);
  cov(1,0) = cov(0,1);
  cov(2,0) = cov(0,2);
  cov(2,1) = cov(1,2);

  // 1Re 
  TFile * f1Re = new TFile("./inputs/official/2012/superk/err/SKMay.root","read");
  TMatrixD * pcov1Re = (TMatrixD*)f1Re->Get("mat_cov_all");
  TMatrixD & cov1Re = *pcov1Re;
  for(int i=0; i<n1Re; i++) {
    for(int j=0; j<n1Re; j++) {
      cov(n1Rmu+i,n1Rmu+j) = cov1Re(i,j);
    }
  }

  //
  // Now work-out correlations between the 1Rmu and 1Re bins
  //

  // size of correlations (1 fully correlated, -1 fully anti-correlated)
  // Start with max (anti-)correlation, then tweak till matrix is invertible.
  double kcovsize = 0;
  double lcovsize = 0;
  double mcovsize = 0;
  if(k!=0) {
    kcovsize = 1. * k;
  }
  if(l!=0) {
    lcovsize = 1. * l;
  }
  if(m!=0) {
    mcovsize = 1. * m;
  }

  bool done = false;
  while(!done) {

     // work-out k-correlations
     cout << "Correlations decided by k index: " << endl;
     if(k!=0) {
       for(int ig=0; ig<ngroups; ig++) {
         for(int jg=0; jg<ngroups; jg++) {
            if(kcov[ig][jg]!= 1) continue;
            if(jg<ig) continue;
            for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
               for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
                  double sigi = TMath::Sqrt(cov[i][i]);
                  double sigj = TMath::Sqrt(cov[j][j]);
                  cov[i][j] = kcovsize * sigi * sigj;
                  cout << "Covar matrix element (" << i << ", " << j << ") set to " 
                       << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
               }//j
            }//i
         }//jg
       }//ig
     }//k!=0

     // work-out l-correlations
     cout << "Correlations decided by l index: " << endl;
     if(l!=0) {
       for(int ig=0; ig<ngroups; ig++) {
         for(int jg=0; jg<ngroups; jg++) {
            if(lcov[ig][jg]!= 1) continue;
            if(jg<ig) continue;
            for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
               for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
                  double sigi = TMath::Sqrt(cov[i][i]);
                  double sigj = TMath::Sqrt(cov[j][j]);
                  cov[i][j] = lcovsize * sigi * sigj;
                  cout << "Covar matrix element (" << i << ", " << j << ") set to " 
                       << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
               }//j
            }//i
         }//jg
       }//ig
     }//k!=0

     // work-out m-correlations
     cout << "Correlations decided by m index: " << endl;
     if(m!=0) {
       for(int ig=0; ig<ngroups; ig++) {
         for(int jg=0; jg<ngroups; jg++) {
            if(mcov[ig][jg]!= 1) continue;
            if(jg<ig) continue;
            for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
               for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
                  double sigi = TMath::Sqrt(cov[i][i]);
                  double sigj = TMath::Sqrt(cov[j][j]);
                  cov[i][j] = mcovsize * sigi * sigj;
                  cout << "Covar matrix element (" << i << ", " << j << ") set to " 
                       << mcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
               }//j
            }//i
         }//jg
       }//ig
     }//k!=0

     // make sure matrix is symmetric
     for(int i=0; i<n; i++) {
       for(int j=0; j<n; j++) {
          if(i<=j) continue;
          cov(i,j) = cov(j,i);
       }
     }
     assert(cov.IsSymmetric());    

     // try to invert
     TMatrixD invcov(cov);
     invcov.Invert();

     TMatrixD id = invcov*cov;
     const double epsilon = 1E-4;
     bool isid = true;
     for(int i=0; i<n; i++) {
       for(int j=0; j<n; j++) {
         if(i==j) { 
           if(TMath::Abs(1-id[i][i]) > epsilon) {
             isid = false;
             break;
           }
         } else { 
           if(TMath::Abs(id[i][j]) > epsilon) {
             isid = false;
             break;
           }
         }
       }//j
     }//i
     if(!isid) {
       print_matrix(id, 10, 6);
     }

     if(isid) { 
       done = true; 
       cout << "The inverse was sucessfully checked. Here is the prod = cov*covInv: " << endl;
       print_matrix(id, 10, 6);
     } else {
       // tweak size of correlations
       cerr << "The matrix is not positive definite - INVERSE CANNOT BE CALCULATED!!! " << endl;
       exit(1);  // not done
     }
  }

  // printout matrix
  cout << "Here is the covariance matrix:" << endl;
  print_matrix(cov, 10, 6);

  // save matrix
  ostringstream fname;
  fname << "cov";
  
  if(k==0) fname << "0";
  else {
    if(k>0) fname << k;//"+";
    if(k<0) fname << "m" << -k;//"-";
  }
  if(l==0) fname << "0";
  else {
    if(l>0) fname << l;//"+";
    if(l<0) fname << "m" << -l; //"-";
  }
  if(m==0) fname << "0";
  else {
    if(m>0) fname << m;//"+";
    if(m<0) fname << "m" << -m ; //"-";
  }
  
  // fname << ".root";
  //  TFile f(fname.str().c_str(),"recreate");
  TFile f("cov.root","update");
  cov.Write(fname.str().c_str());
  f.Close();
}

void print_matrix(TMatrixD & m, int w, int prec)
{
  assert(m.GetNcols() == m.GetNrows());
  int n =  m.GetNrows();
  cout << fixed;
  for(int i=0; i<n; i++) {
    if(i==0) {
       cout << setw(w) << setfill(' ') << "r\\c:" << " | ";
       for(int j=0; j<n; j++) {
         cout << setw(w) << setfill(' ') << setprecision(0) << j << " | ";
       }
       cout << endl;
       cout << setw((w+3)*(n+1)) << setfill('-') << " " << endl;
    }
    if(i==7) {
       cout << setw((w+3)*(n+1)) << setfill('-') << " " << endl;
    }
    cout << setw(w) << setfill(' ') << setprecision(0) << i << " | ";
    for(int j=0; j<n; j++) {
         cout << setw(w) << setfill(' ') << setprecision(prec)
              << m[i][j] << " | ";
    }//rows
    cout << endl;
  }//columns
  cout << scientific;

}
