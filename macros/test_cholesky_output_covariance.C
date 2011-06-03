//
// ------------------------------------------------------------------------------------------------
// Test the Cholesky method
//
// Generate a random default-param vector (np - dimensions) and a random covariance
// matrix (np X np - dimensional). Run the Cholesky method N times and calculate the 
// parameter covariance matrix using the N output tweaked parameter vectors. 
// Compare the calculated covariance matrix with the one that was used as input 
// for the Cholesky method.
//
// Costas Andreopoulos
// ------------------------------------------------------------------------------------------------
//

#include <iomanip.h>

void test_cholesky_output_covariance(
  const int np=6, const char * inpfile="inputs/test/cholesky_inputs_6dim_1.txt")
{
  using namespace std;
  using namespace t2kana;

  //
  // read inputs for the Cholesky method 
  //

  TVectorD input_mean_params (np);
  TMatrixD input_cov_matrix  (np,np);

  ifstream finp(inpfile);
  int iline=0;
  while(1) {
    char c = finp.peek();
    char comment[1024];
    if(c == '#' || c == '\n') {
       finp.getline(comment,1024);
     //cout << comment << endl;
       continue;
    }
    double val;
    if(iline==0) {
      for(int jp=0; jp<np; jp++) {
        finp >> val; 
        input_mean_params[jp] = val;
      //cout << "Reading input_mean_params(" << jp << ") = " << input_mean_params[jp] << endl;
      }
    } else {
      int ip = iline-1;
      for(int jp=0; jp<np; jp++) {
        finp >> val; 
        input_cov_matrix[ip][jp] = val;
      //cout << "Reading input_cov_matrix(" << ip << ", " << jp << ") = " << input_cov_matrix[ip][jp] << endl;
      //cout << val << endl;
      }
    }
    iline++;
    if(iline == np+1) break;
  }//eof

  cout << "Covariant matrix determinant = " << input_cov_matrix.Determinant() << endl;

  //
  // run the Cholesky decomposition method to obtain the triangular
  // matrix that can be used for Monte Carlo generation of correlated params
  //
  TMatrixD cholesky_triangular = util::CholeskyDecomposition(input_cov_matrix); 

  //
  // run the Cholesky method N times and store the tweaked output vectors
  //
  const int N = 200000;
  TMatrixD x(N,np);
  for(int i=0; i<N; i++) {
     if(i%100) {
       cout << "Number of vectors generated so far: " << i << endl;
     }
     TVectorD twk_vector = util::CholeskyGenerateCorrelatedParams(cholesky_triangular, input_mean_params);
     for(int ip=0; ip<np; ip++) {
       x[i][ip] = twk_vector[ip];
     }//ip  
  }//i

  //
  // calculate the expectation, correlation & covariance matrix
  //

  TVectorD calc_mean_params (np);
  TMatrixD calc_cov_matrix  (np,np);
  TMatrixD calc_cor_matrix  (np,np);

  for(int ip=0; ip<np; ip++) {
    calc_mean_params[ip] = 0.;
    for(int i=0; i<N; i++) {
      calc_mean_params[ip] += x[i][ip];
    }//i
    calc_mean_params[ip] /= N;
  }//ip

  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
      double cov=0;
      for(int i=0; i<N; i++) {
         double dxip = x[i][ip] - calc_mean_params[ip];
         double dxjp = x[i][jp] - calc_mean_params[jp];
         cov += (dxip*dxjp);
      }//i
      cov /= N;
      calc_cov_matrix[ip][jp] = cov;
    }//jp
  }//ip

  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
      double sigip =  TMath::Sqrt(calc_cov_matrix[ip][ip]);
      double sigjp =  TMath::Sqrt(calc_cov_matrix[jp][jp]);
      calc_cor_matrix[ip][jp] = calc_cov_matrix[ip][jp] / (sigip*sigjp);
    }//jp
  }//ip

  //
  // print-out so as to compare the covariance matrices
  //

  cout << "\n\nCalculated mean from generated vectors:" << endl;
  for(int ip=0; ip<np; ip++) {
    cout << setw(15) << setfill(' ') << setprecision(6) << calc_mean_params[ip] << " | ";
  }
  cout << "\n\nInput mean" << endl;
  for(int ip=0; ip<np; ip++) {
    cout << setw(15) << setfill(' ') << setprecision(6) << input_mean_params[ip] << " | ";
  }
  cout << "\n\nCorrelation matrix calculated from generated vectors:" << endl;
  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
         cout << setw(15) << setfill(' ') << setprecision(6) << calc_cor_matrix[ip][jp] << " | ";
    }//jp
    cout << endl;
  }//ip
  cout << "\n\nCovariance matrix calculated from generated vectors:" << endl;
  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
         cout << setw(15) << setfill(' ') << setprecision(6) << calc_cov_matrix[ip][jp] << " | ";
    }//jp
    cout << endl;
  }//ip
  cout << "\n\nInput covariance matrix:" << endl;
  for(int ip=0; ip<np; ip++) {
    for(int jp=0; jp<np; jp++) {
         cout << setw(15) << setfill(' ') << setprecision(6) << input_cov_matrix[ip][jp] << " | ";
    }//jp
    cout << endl;
  }//ip

  cout << "Done!" << endl;
}
