// Test the Cholesky method
//
// Generate a random default-param vector (np - dimensions) and a random covariance
// matrix (np X np - dimensional). Run the Cholesky method N times and calculate the 
// parameter covariance matrix using the N output tweaked parameter vectors. 
// Compare the calculated covariance matrix with the one that was used as input 
// for the Cholesky method.
//
// C.Andreopoulos

{
  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libMinuit");   
  gSystem->Load("libMathCore");   
  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  //
  // generate inputs for the Cholesky method 
  //

  const int np = 5;

  TVectorD def_vector (np);
  TMatrixD cov_matrix (np,np);

  TRandom3 random(0);

  for(int ip=0; ip<np; ip++) {
    def_vector[ip] = random.Gaus();
    for(int jp=0; jp<np; jp++) {
      cov_matrix[ip][jp] = -1 + 2.*random.Rndm();
    }//jp
  }//ip

  //
  // run the Cholesky method N times and store the tweaked output vectors
  //

  const int N = 1000;

  TMatrixD x(N,np);

  for(int i=0; i<N; i++) {
     cout << "Iteration: " << i << " --> ";
     TVectorD twk_vector = util::RunCholeskyThingy(cov_matrix, def_vector);
     for(int ip=0; ip<np; ip++) {
       cout << " p[" << ip << "] = " << twk_vector[ip];
       x[i][ip] = twk_vector[ip];
     }//ip  
     cout << endl;
  }//i

  //
  // calculate the covariance matrix
  //

  TMatrixD calc_cov_matrix (np,np);

  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
      double cov=0;
      for(int i=0; i<N; i++) {
         double dxip = x[i][ip] - def_vector[ip];
         double dxjp = x[i][jp] - def_vector[jp];
         cov += (dxip*dxjp);
      }//i
      cov /= N;
      calc_cov_matrix[ip][jp] = cov;
    }//jp
  }//ip

  //
  // compare the input covariance matrix and the one calculated 
  // from the generated vectors
  //

  cout << "\n\nComparing ovariance matrices" << endl;
  cout << "Printing out `calculated(input)' covariance matrix elements\n" << endl;
  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
         cout << calc_cov_matrix[ip][jp] << "(" << cov_matrix[ip][jp] << ")\t";
    }//jp
    cout << endl;
  }//ip

}
