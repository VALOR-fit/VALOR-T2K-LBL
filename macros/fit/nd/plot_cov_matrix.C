/*
This macro reads covariance matrix/matrices from the input ROOT file and makes a TH2D(s) of the matrix/matrices 
and a TH1D(s) of the square root(s) of its/their diagonal elements. It writes these histograms to an output ROOT 
file. The input ROOT file is made by generating and fitting some toy experiments using run_fit.C.

There is an option to marginalise the input covariance matrix/matrices. This can only be done if the input matrix/
matrices have dimensions n_fit_params x n_fit_params. The marginalisation is done by deleting the rows and 
columns that correspond to the parameters to be marginalised. This is equivalent to integrating over the parameters 
to be marginalised if the parameters have a multivariate Gaussian distribution - please see, for example

section 3.1 in http://www.cs.uoi.gr/~arly/full/conferences/C69.pdf

and

section 3.1 in http://arxiv.org/pdf/0906.4123.pdf

The parameters to be marginalised are taken from the enum valor::conv::kNDFitParamMarginalize.

If there is only one entry in the input ROOT tree, the covariance matrix and the square roots of its diagonal 
elements are drawn.

Nick Grant  21 July 2014 
*/

void plot_cov_matrix(
TString inputfile = "/opt/ppd/t2k/users/nick/DUNE/macros/fit/nd/fit51params/toymc_fhcflux_plus10_4.root",
bool marginalise = true,
TString config = "analysis/dune_tdr.xml"
)
{
  using namespace valor;
  using namespace valor::conv;

  // Read the config file to populate the ParameterPool
  std::string config_path = util::CompleteConfigPath(config.Data());
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(config_path.c_str());
  const std::vector<Parameter> & params = ParameterPool::Instance()->GetParams();
  const int n_fit_params = ParameterPool::Instance()->GetNumSystParams();

  // Open input file
  TFile f(inputfile, "read");

  PhysParamList * true_param_array     = 0; // input true params
  FitParamList *  best_fit_param_array = 0; // best-fit params
  TObjArray *     best_fit_model_array = 0; // an array of various best-fit model predictions for the current toy-MC experiment 
  TObjArray *     data_array           = 0; // an array of all fitted data distributions for the current toy-MC experiment

  TTree * fitnt = (TTree*)f.Get("fitnt");
  fitnt->SetBranchAddress("true_params",        &true_param_array);
  fitnt->SetBranchAddress("best_fit_params",    &best_fit_param_array);
  fitnt->SetBranchAddress("best_fit_model",     &best_fit_model_array);
  fitnt->SetBranchAddress("data",               &data_array);

  int num_params, i, j;
  double element, x_low, x_high;
  int n_params_after_marg;

  if(marginalise)
    {
    n_params_after_marg = 0;
    //find number of parameters after marginalisation
    for(int i=0; i<n_fit_params; i++)
      {
      if(!kNDFitParamMarginalize[i])
        n_params_after_marg++;
      }

    num_params = n_params_after_marg;
    cout<<endl;
    cout<<"Number of parameters after marginalisation is "<<n_params_after_marg<<"."<<endl;
    cout<<endl;
    }
  else
    {
    fitnt->GetEntry(0);
    if(best_fit_param_array->GetCovMatrix()->GetNrows() != best_fit_param_array->GetCovMatrix()->GetNcols())
      {
      cerr<<"Covariance matrix in experiment 0 is not square, it has "<<best_fit_param_array->GetCovMatrix()->GetNrows()
         << " rows and "<<best_fit_param_array->GetCovMatrix()->GetNcols()<<" columns."<<endl;
      exit(1);
      }
    num_params = best_fit_param_array->GetCovMatrix()->GetNrows();
    }

  const int n_params = num_params;

  int num_entries = fitnt->GetEntries();
  const int n_entries = num_entries;

  TString histname;
  TH2D *hCovMatrix[n_entries] = {0};
  for(int i=0; i<n_entries; i++)
    {
    histname = Form("hCovMatrix_%d", i);
    hCovMatrix[i] = new TH2D(histname, "", n_params, 0, n_params, n_params, 0, n_params);
    //need this line if any of these histograms is drawn in this macro 
    hCovMatrix[i]->SetDirectory(gROOT);
    }

  TH1D *hCovMatDiagonal[n_entries] = {0};
  for(int i=0; i<n_entries; i++)
    {
    histname = Form("hCovMatDiagonal_%d", i);
    hCovMatDiagonal[i] = new TH1D(histname, "", n_params, 0, n_params);
    //need this line if any of these histograms is drawn in this macro
    hCovMatDiagonal[i]->SetDirectory(gROOT);
    }

  TPaveText *tCovMatDiagonal[n_params] = {0};
  TText *axislabel[n_params] = {0};

  TFile *covmatfile = new TFile("CovMatrix.root", "RECREATE");
  covmatfile->cd();

  // Loop over entries (toy-MC experiments) in the tree
  for(Long_t toy=0; toy < fitnt->GetEntries(); toy++) 
    {
    fitnt->GetEntry(toy);

    if(toy%100 == 0)
      {
      cout<<endl;
      cout<<"Doing toy experiment "<<toy<<endl;
      cout<<endl;
      }

    TMatrixDSym * cov = best_fit_param_array->GetCovMatrix();

    cov->Print();

    if(cov->GetNrows() != cov->GetNcols())
      {
      cerr<<"Covariance matrix in experiment "<<toy<<" is not square, it has "<<cov->GetNrows()<< " rows and "<<cov->GetNcols()<<" columns."<<endl;
      exit(1);
      }

    if(marginalise)
      {
      //check that input covariance matrix has correct dimensions
      if(cov->GetNrows() != n_fit_params || cov->GetNcols() != n_fit_params)
        {
        cerr<<"Input covariance matrix in experiment "<<toy<<" has wrong dimensions to be marginalised - it must be "<<n_fit_params<<" x "<<n_fit_params<<"."<<endl;
        exit(1);
        }
      }//end of if(marginalise)

    double *covmat_array = cov->GetMatrixArray();
    if(marginalise)
      double *end_covmat_array = covmat_array + n_fit_params * n_fit_params;
    else
      double *end_covmat_array = covmat_array + n_params * n_params;

    i = -1;
    j = -1;

    for(; covmat_array<end_covmat_array; covmat_array++)
      { 
      //get element from input covariance matrix
      element = *covmat_array;
      //need to find 2D indices to write element to marginalised covariance matrix
      if(marginalise)
	{
        if((j+1)%n_fit_params==0)
	  {
	  i++;
	  j = -1;
	  }
	}
      else
	{
	if((j+1)%n_params==0)
	  {
	  i++;
	  j = -1;
	  }
	}
      j++;
      //skip this element if it is in a row or column corresponding to a parameter to be marginalised
      if(marginalise && (kNDFitParamMarginalize[i] || kNDFitParamMarginalize[j]))
        continue;
      //write element to marginalised covariance matrix
      hCovMatrix[toy]->SetBinContent(i+1, j+1, element);
      if(i == j)
	{
	hCovMatDiagonal[toy]->SetBinContent(i+1, sqrt(element));       
	}
      }

    hCovMatrix[toy]->Write();
    hCovMatDiagonal[toy]->Write();

    }//end of for(Long_t toy; toy < fitnt->GetEntries(); toy++)

  covmatfile->Close();
  
  if(n_entries == 1)
    {
    TCanvas *CovMatrix = new TCanvas("CovMatrix", "CovMatrix", 800, 600);
    CovMatrix->Divide(1);
    CovMatrix->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.14);
    gPad->SetLogz();
    hCovMatrix[0]->GetXaxis()->SetTitle("Parameter");
    hCovMatrix[0]->GetYaxis()->SetTitle("Parameter");
    hCovMatrix[0]->GetXaxis()->SetTitleOffset(1.15);
    hCovMatrix[0]->GetYaxis()->SetTitleOffset(1.2);
    hCovMatrix[0]->GetXaxis()->CenterTitle();
    hCovMatrix[0]->GetYaxis()->CenterTitle();
    hCovMatrix[0]->Draw("colz");

    TCanvas *CovMatDiagonal = new TCanvas("CovMatDiagonal", "CovMatDiagonal", 800, 600);
    CovMatDiagonal->Divide(1);
    CovMatDiagonal->cd(1);
    gPad->SetLeftMargin(0.18);
    gPad->SetBottomMargin(0.35);
    hCovMatDiagonal[0]->GetXaxis()->SetNdivisions(000);
    hCovMatDiagonal[0]->GetYaxis()->SetRangeUser(0.0, 1.1 * hCovMatDiagonal[0]->GetMaximum());
    hCovMatDiagonal[0]->GetYaxis()->SetTitle("#sqrt{Diagonal of cov matrix}");
    hCovMatDiagonal[0]->GetYaxis()->SetTitleOffset(1.2);
    hCovMatDiagonal[0]->GetYaxis()->CenterTitle();
    hCovMatDiagonal[0]->SetLineColor(2);
    hCovMatDiagonal[0]->SetLineWidth(2);
    hCovMatDiagonal[0]->Draw();

    for(Parameter param : params)
      {
      PhysParamId_t p = param.GetID();
      if(!util::IsSystematic(p)) continue;
      //const char * label = util::AsStringTeX(p,true).c_str();
      const char * label = util::NameRegistry::Instance()->ParamAsString(p).c_str();
      x_low = param;
      x_high = x_low + 1.0;
      tCovMatDiagonal[param] = new TPaveText(x_low,-0.04,x_high,-0.001);
      axislabel[param] = tCovMatDiagonal[param]->AddText(label);
      axislabel[param]->SetTextAngle(90);
      //SetTextAlign(11) to left justify, 21 to centre and 31 to right justify
      axislabel[param]->SetTextAlign(31);
      axislabel[param]->SetTextSize(0.03);
      tCovMatDiagonal[param]->SetFillColor(0);
      tCovMatDiagonal[param]->SetBorderSize(0);
      tCovMatDiagonal[param]->Draw();
      }

   }
  
}//end of plot_cov_matrix
