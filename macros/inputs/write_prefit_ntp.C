
void write_prefit_ntp(
  TString config = "fit/dune_ndfit.xml",
  TString outfilename = "prefit_ntp.root",
  TString hypname = "3f_osc",
  TString fit_param_opt = "no-osc-params"
)
{
  using namespace std;
  using namespace valor;
  util::PrintBanner(Env::Instance()->Banner(),500);
  
  valor::StdFitNtp ntp;
  ntp.Create(outfilename.Data(), "valor_dune_nd_ntp", "VALOR DUNE / ND-Fit prior");
  
  string configpath = util::CompleteConfigPath(config.Data());
  // Get ND model & configure
  cout<< "Getting and configuring the fit elements"<<endl;
  CompoundFitElements * elements = new CompoundFitElements();
  elements->Read(configpath.c_str(),"nomc nosplines");
  FitModelABC * model = elements->FitModel();
  
  // Reset output tree variables
  ntp.Reset();
  
  // Get the systematics object
  CorrelatedSystematics * syst = model->GetSystematics();
  // Get the specified fit parameter list (specifies which parameters are floated/fixed, 
  // parameter limits, initial values etc).
  FitParamList * fp = util::GetFitParameterList(hypname.Data(), fit_param_opt);
  PhysParamList true_params = fp->AsPhysParamList("current");
  ntp.TrueParamArr()->Append(true_params);
  
  // Copy bestfit
  fp->CopyCurrToBestFitValues();
  
  assert(syst->Size()==fp->Size());
  TMatrixDSym prior_cov = util::MatrixToSymMatrix(syst->CovMatrix());
  fp->SetCovMatrix(prior_cov);

  // Copy best-fit params
  ntp.BestFitParamArr()->Copy(*fp);
  ntp.Fill();
  
  ntp.Save();
  
  cout<<"Done! Written to "<<outfilename<<endl;
}
  
