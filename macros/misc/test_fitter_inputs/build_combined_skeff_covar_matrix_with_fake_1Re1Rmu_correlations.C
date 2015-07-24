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
  TFile * f1Re = new TFile("./SKMay.root","read");
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

  double k_values[15];
  double l_values[15];
  double m_values[15];

  for (int i = 0; i<15; i++){
    k_values[i] = 0.0;
    l_values[i] = 0.0;
    m_values[i] = 0.0;
  }

  // size of correlations (1 fully correlated, -1 fully anti-correlated)
  // Start with max (anti-)correlation, then tweak till matrix is invertible.
  // double kcovsize = 0;
  //double lcovsize = 0;
  //double mcovsize = 0;
 
  //we start in 0 for all
  /* if(k!=0) {
    kcovsize = 1. * k;
  }
  if(l!=0) {
    lcovsize = 1. * l;
  }
  if(m!=0) {
    mcovsize = 1. * m;
    }*/

  // FIRST "LOOP" - MOVING ONLY K
  double kmaxaux = -999;
  double kminaux =  999;

  int nsteps = 201;

  //while(!done) {
  for (int ik = 0; ik<201; ik++){
    double kcovsize = -1 + (ik)*0.01;
    cout << "kcovsize = " << kcovsize << endl;
    bool done = false;    
    // work-out k-correlations
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    //  }//k!=0
    
    // make sure matrix is symmetric
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
	if(i<=j) continue;
	cov(i,j) = cov(j,i);
      }
    }
    assert(cov.IsSymmetric());    
    // MATRIX IS BUILT
    
    //FIRST CHECK -> POSITIVE DEFINITE
    
    //Cholesky-decompose the matrix to check for positive-definiteness
    TMatrixD  SymCov(cov);
    TDecompChol choli(SymCov);
    Bool_t dek; 
    dek = choli.Decompose(); 
    if (!dek) {
      cerr << "matrix not positive definite" << endl;
      done = false;
      // exit(1);
    }
    else {
      done = true;
      cout << "matrix IS positive definite!!!! You may now continue and invert it" << endl;
    }
    //For reference, the standard method used for testing positive-definiteness, as used by ROOT:
    //see also http://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky.E2.80.93Banachiewicz_and_Cholesky.E2.80.93Crout_algorithms
    //for a theory exposition of the algorithm
    /*
      Int_t i,j,icol,irow;
      const Int_t     n  = cov.GetNrows();
      Double_t *pU = cov.GetMatrixArray();
      for (icol = 0; icol < n; icol++) {
      const Int_t rowOff = icol*n;
      
      //Compute fU(j,j) and test for non-positive-definiteness.
      Double_t ujj = pU[rowOff+icol];
      for (irow = 0; irow < icol; irow++) {
      const Int_t pos_ij = irow*n+icol;
      ujj -= pU[pos_ij]*pU[pos_ij];
      }
      if (ujj <= 0) { NOT positive definite}
      ...(actually do the decomp)
      } */
    
    // try to invert - ONLY IF IT IS POSITIVE DEFINITE
    if (done){
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
	 // print_matrix(id, 10, 6);
       } else {
	 // tweak size of correlations
	 cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	 // exit(1);  // not done
	 done = false;  // can be positve definite but not invertible...
       }
     }
    
    
    if(done){
      if(kcovsize!=0){
	if ((kcovsize>0) && (kcovsize>kmaxaux)) kmaxaux = kcovsize;
	if ((kcovsize<0) && (kcovsize<kminaux)) kminaux = kcovsize;
      } 
      else done = false; 
    }
    
  }//end of first "loop"

  k_values[1] = kmaxaux;
  k_values[2] = kminaux;

  //SECOND "LOOP" - 1
  double lmaxaux = -999;
  double lminaux =  999;
  
  double found = false;
  done = false;
  kcovsize = k_values[1];
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int il = 0; il<201; il++){
      double lcovsize = -1 + (il)*0.01;
      cout << "lcovsize = " << lcovsize << endl;
      bool done = false;    
      // work-out l-correlations
      cout << "Correlations decided by l index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(lcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = lcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!!  You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(lcovsize!=0){
	  if ((lcovsize>0) && (lcovsize>lmaxaux)) lmaxaux = lcovsize;
	  if ((lcovsize<0) && (lcovsize<lminaux)) lminaux = lcovsize;
      } 
	else done = false; 
      }
      
    }//end of second "loop"
    if ((lmaxaux!=-999) && (lminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of k" << endl;
      kcovsize -= 0.01;//if not done, change k	
      if (kcovsize < 0) {
	cout << "Not reached scenario with kmax != 0 and l definite positve!" << endl;
	continue;
      }
    }
  }
  k_values[3] = kcovsize;
  k_values[4] = kcovsize;
  l_values[3] = lmaxaux;
  l_values[4] = lminaux;


  //SECOND "LOOP" - 2
  lmaxaux = -999;
  lminaux =  999;
  
  found = false;
  done = false;
  kcovsize = k_values[2];
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int il = 0; il<201; il++){
      double lcovsize = -1 + (il)*0.01;
      cout << "lcovsize = " << lcovsize << endl;
      bool done = false;    
      // work-out l-correlations
      cout << "Correlations decided by l index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(lcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = lcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!! You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(lcovsize!=0){
	  if ((lcovsize>0) && (lcovsize>lmaxaux)) lmaxaux = lcovsize;
	  if ((lcovsize<0) && (lcovsize<lminaux)) lminaux = lcovsize;
	} 
	else done = false; 
      }
      
    }//end of second "loop"
    if ((lmaxaux!=-999) && (lminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of k" << endl;
      kcovsize += 0.01;//if not done, change k	
      if (kcovsize < 0) {
	cout << "Not reached scenario with kmax != 0 and l definite positve!" << endl;
	continue;
      }
    }
  }
  k_values[5] = kcovsize;
  k_values[6] = kcovsize;
  l_values[5] = lmaxaux;
  l_values[6] = lminaux;



  //THIRD "LOOP" - 1
  double mmaxaux = -999;
  double mminaux =  999;
  
  found = false;
  done = false;
  //PART 1 - kmax + lmax (latest)
  kcovsize = k_values[3];
  lcovsize = l_values[3];
  
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig


    //fixing l to lmax in previous
    cout << "Correlations decided by l index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(lcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = lcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int im = 0; im<201; im++){
      double mcovsize = -1 + (im)*0.01;
      cout << "mcovsize = " << mcovsize << endl;
      bool done = false;    
      // work-out m-correlations
      cout << "Correlations decided by m index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(mcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = mcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!! You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(mcovsize!=0){
	  if ((mcovsize>0) && (mcovsize>mmaxaux)) mmaxaux = mcovsize;
	  if ((mcovsize<0) && (mcovsize<mminaux)) mminaux = mcovsize;
	}
	else done = false;
      }
      
    }//end of second "loop"
    if ((mmaxaux!=-999) && (mminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of l" << endl;
      lcovsize -= 0.01;//if not done, change l	
      if (lcovsize < 0) {
	cout << "Changing value of k" << endl;
	lcovsize = 1.0;
	kcovsize -= 0.01;//if not done, change l
	if (kcovsize<0) continue;
      }
    }
  }
  k_values[7] = kcovsize;
  k_values[8] = kcovsize;
  l_values[7] = lcovsize;
  l_values[8] = lcovsize;
  m_values[7] = mmaxaux;
  m_values[8] = mminaux;


  //THIRD "LOOP" - 2
  mmaxaux = -999;
  mminaux =  999;
  
  found = false;
  done = false;
  //PART 2 - kmax + lmin (latest)
  kcovsize = k_values[3];
  lcovsize = l_values[6];
  
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig


    //fixing l to lmax in previous
    cout << "Correlations decided by l index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(lcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = lcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int im = 0; im<201; im++){
      double mcovsize = -1 + (im)*0.01;
      cout << "mcovsize = " << mcovsize << endl;
      bool done = false;    
      // work-out m-correlations
      cout << "Correlations decided by m index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(mcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = mcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!!  You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(mcovsize!=0){
	  if ((mcovsize>0) && (mcovsize>mmaxaux)) mmaxaux = mcovsize;
	  if ((mcovsize<0) && (mcovsize<mminaux)) mminaux = mcovsize;
	}
	else done = false;
      }
      
    }//end of second "loop"
    if ((mmaxaux!=-999) && (mminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of l" << endl;
      lcovsize += 0.01;//if not done, change l	
      if (lcovsize > 0) {
	cout << "Changing value of k" << endl;
	lcovsize = -1.0;
	kcovsize -= 0.01;//if not done, change l
	if (kcovsize<0) continue;
      }
    }
  }
  k_values[9] = kcovsize;
  k_values[10] = kcovsize;
  l_values[9] = lcovsize;
  l_values[10] = lcovsize;
  m_values[9] = mmaxaux;
  m_values[10] = mminaux;


  //THIRD "LOOP" - 3
  mmaxaux = -999;
  mminaux =  999;
  
  found = false;
  done = false;
  //PART 3 - kmax + lmin (latest)
  kcovsize = k_values[6];
  lcovsize = l_values[3];
  
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig


    //fixing l to lmax in previous
    cout << "Correlations decided by l index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(lcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = lcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int im = 0; im<201; im++){
      double mcovsize = -1 + (im)*0.01;
      cout << "mcovsize = " << mcovsize << endl;
      bool done = false;    
      // work-out m-correlations
      cout << "Correlations decided by m index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(mcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = mcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!! You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(mcovsize!=0){
	  if ((mcovsize>0) && (mcovsize>mmaxaux)) mmaxaux = mcovsize;
	  if ((mcovsize<0) && (mcovsize<mminaux)) mminaux = mcovsize;
	}
	else done = false;
      }
      
    }//end of second "loop"
    if ((mmaxaux!=-999) && (mminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of l" << endl;
      lcovsize -= 0.01;//if not done, change l	
      if (lcovsize < 0) {
	cout << "Changing value of k" << endl;
	lcovsize = 1.0;
	kcovsize += 0.01;//if not done, change l
	if (kcovsize>0) continue;
      }
    }
  }
  k_values[11] = kcovsize;
  k_values[12] = kcovsize;
  l_values[11] = lcovsize;
  l_values[12] = lcovsize;
  m_values[11] = mmaxaux;
  m_values[12] = mminaux;
  
  //THIRD "LOOP" - 4
  mmaxaux = -999;
  mminaux =  999;
  
  found = false;
  done = false;
  //PART 4 - kmin + lmin (latest)
  kcovsize = k_values[6];
  lcovsize = l_values[6];
  
  while(!found) {
    
    //fixing k to kmax in previous
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig


    //fixing l to lmax in previous
    cout << "Correlations decided by l index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(lcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = lcovsize * sigi * sigj;
	    //cout << "Covar matrix element (" << i << ", " << j << ") set to " 
	    //	 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    for (int im = 0; im<201; im++){
      double mcovsize = -1 + (im)*0.01;
      cout << "mcovsize = " << mcovsize << endl;
      bool done = false;    
      // work-out m-correlations
      cout << "Correlations decided by m index: " << endl;
      for(int ig=0; ig<ngroups; ig++) {
	for(int jg=0; jg<ngroups; jg++) {
	  if(mcov[ig][jg]!= 1) continue;
	  if(jg<ig) continue;
	  for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	    for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	      double sigi = TMath::Sqrt(cov[i][i]);
	      double sigj = TMath::Sqrt(cov[j][j]);
	      cov[i][j] = mcovsize * sigi * sigj;
	    }//j
	  }//i
	}//jg
      }//ig
      
      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    
      // MATRIX IS BUILT
      
      //FIRST CHECK -> POSITIVE DEFINITE
      
      //Cholesky-decompose the matrix to check for positive-definiteness
      TMatrixD  SymCov(cov);
      TDecompChol choli(SymCov);
      Bool_t dek; 
      dek = choli.Decompose(); 
      if (!dek) {
	cerr << "matrix not positive definite" << endl;
	done = false;
	// exit(1);
      }
      else {
      done = true;
      cout << "matrix IS positive definite!!!! You may now continue and invert it" << endl;
      }
      
      // try to invert - ONLY IF IT IS POSITIVE DEFINITE
      if (done){
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
	  // print_matrix(id, 10, 6);
	} else {
	  // tweak size of correlations
	  cerr << "INVERSE CANNOT BE CALCULATED!!! " << endl;
	  // exit(1);  // not done
	  done = false;  // can be positve definite but not invertible...
	}
      }
      if(done){
	if(mcovsize!=0){
	  if ((mcovsize>0) && (mcovsize>mmaxaux)) mmaxaux = mcovsize;
	  if ((mcovsize<0) && (mcovsize<mminaux)) mminaux = mcovsize;
	}
	else done = false;
      }
      
    }//end of second "loop"
    if ((mmaxaux!=-999) && (mminaux!=999)) found = true;
    if(!found){
      cout << "Changing value of l" << endl;
      lcovsize += 0.01;//if not done, change l	
      if (lcovsize > 0) {
	cout << "Changing value of k" << endl;
	lcovsize = -1.0;
	kcovsize += 0.01;//if not done, change l
	if (kcovsize>0) continue;
      }
    }
  }
  k_values[13] = kcovsize;
  k_values[14] = kcovsize;
  l_values[13] = lcovsize;
  l_values[14] = lcovsize;
  m_values[13] = mmaxaux;
  m_values[14] = mminaux;
  


  //Final print out 
  for (int i = 0; i<15; i++){
    cout << " k_values[" << i << "] = " << k_values[i] << " l_values[" << i<< "] = " << l_values[i] << " m_values[" << i<< "] = " << m_values[i] << endl;
  }

  //save matrices
  

  TFile f("cov_new.root","update");
  
  double kk = 0;
  double ll = 0;
  double mm = 0;

  for (int index = 0; index<15; index++){
    kk = k_values[index];
    ll = l_values[index];
    mm = m_values[index];
  
    cout << "kk = " << kk << " , ll = " << ll << " and mm = " << mm << endl;
    cout << "Correlations decided by k index: " << endl;
    //    if(k!=0) {
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(kcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = kk * sigi * sigj;
	    cout << "Covar matrix element (" << i << ", " << j << ") set to " <<  cov[i][j]  << endl;
	      //	    	 << kcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    

    //fixing l to lmax in previous
    cout << "Correlations decided by l index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(lcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = ll * sigi * sigj;
	    cout << "Covar matrix element (" << i << ", " << j << ") set to " <<  cov[i][j]  << endl;
	    //	    	 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    
    
    //fixing m to lmax in previous
    cout << "Correlations decided by m index: " << endl;
    for(int ig=0; ig<ngroups; ig++) {
      for(int jg=0; jg<ngroups; jg++) {
	if(mcov[ig][jg]!= 1) continue;
	if(jg<ig) continue;
	for(int i=idx_range[ig][0]; i<=idx_range[ig][1]; i++) {
	  for(int j=idx_range[jg][0]; j<=idx_range[jg][1]; j++) {
	    double sigi = TMath::Sqrt(cov[i][i]);
	    double sigj = TMath::Sqrt(cov[j][j]);
	    cov[i][j] = mm * sigi * sigj;
	    cout << "Covar matrix element (" << i << ", " << j << ") set to " <<  cov[i][j]  << endl;
	    //		 << lcovsize << " * " << 100*sigi << "% * " << 100*sigj << "%" << endl;
	  }//j
	}//i
      }//jg
    }//ig
    

      // make sure matrix is symmetric
      for(int i=0; i<n; i++) {
	for(int j=0; j<n; j++) {
	  if(i<=j) continue;
	  cov(i,j) = cov(j,i);
	}
      }
      assert(cov.IsSymmetric());    

    // printout matrix
    cout << "Here is the covariance matrix:" << endl;
    //   print_matrix(cov, 10, 6);
    
    // save matrix
    ostringstream fname;
    fname << "cov";
    
    if(kk==0) fname << "0";
    else {
      if(kk>0) fname << "1";//"+";
      if(kk<0) fname << "m1";//"m" << -k;//"-";
    }
    if(ll==0) fname << "0";
    else {
      if(ll>0) fname << "1";//l;//"+";
      if(ll<0) fname << "m1";// << -l; //"-";
    }
    if(mm==0) fname << "0";
    else {
      if(mm>0) fname << "1";//m;//"+";
      if(mm<0) fname << "m1" ;//<< -m ; //"-";
    }
    

    cout << "name = " << fname.str().c_str() << endl;
    // fname << ".root";
    //  TFile f(fname.str().c_str(),"recreate");
    
    cov.Write(fname.str().c_str());
  }
  //}
  //}	

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
