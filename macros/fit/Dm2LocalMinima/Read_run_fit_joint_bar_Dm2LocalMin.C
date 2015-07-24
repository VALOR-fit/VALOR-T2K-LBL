
using namespace std;
using namespace valor;
using namespace valor::conv;

void Read_run_fit_joint_bar_Dm2LocalMin
(
 //const char * infile = "$T2KOSCANA/outputs/test/toyexptfits_Dm2Loc_10-run_0-3f+3fbar.root",
 const char * infile = "$T2KOSCANA/outputs/sensi/RHCmode_FixNu/Dm2localminima_0p170e21/toyexptfits_Dm2Loc_1kToys-run_0-3f+3fbar.root"
 )
{
  const double Dm2Thr = 0.005; // Above Dm2=Dm2Thr we consider it a second minimum
  const double eps = 0.001;
  const int NminMy = 9;
  const int Dm2True_Idx = 0;

  TH1F *hDm2Loc = new TH1F("hDm2Loc","hDm2Loc",100,0.001,0.04);

  int NOtherMin = 0;

  TFile *fin = new TFile(infile,"READ");
  TTree *fitnt = (TTree*)fin->Get("fitnt");
  
  // Create tree to get outputs

  int Nmin;
  double Dm2LocalMin[NminMy];
  double Dm2LocalMinValVec[NminMy];
  double Dm2LocalMinChisqVec[NminMy];
  PhysParamList * input_params; //only used to initialise branch

  fitnt->SetBranchAddress("Nmin",&Nmin);
  fitnt->SetBranchAddress("Dm2LocalMin",&Dm2LocalMin);
  fitnt->SetBranchAddress("Dm2LocalMinValVec",&Dm2LocalMinValVec);
  fitnt->SetBranchAddress("Dm2LocalMinChisqVec",&Dm2LocalMinChisqVec);
  fitnt->SetBranchAddress("input_params",&input_params);
  
  int NEntries = fitnt->GetEntries();  
  for(int itoy=0; itoy<NEntries; itoy++) {
    fitnt->GetEntry(itoy);

    if(Nmin != NminMy){
      cerr << endl;
      cerr << "Nmin != NminMy" << endl;
      cerr << "Nmin = " << Nmin << endl;
      cerr << "NminMy = " << NminMy << endl;
      cerr << endl;
      exit(1);
    }
  
    int ChisqBF_Idx = -999;
    double ChisqBF = 9999.;
    
    int Dm2BF_Idx = -999;
    double Dm2BF = 9999.;
     
    for(int imin=0;imin<Nmin;imin++){
      
      // search for the bf chisq
      
      if( Dm2LocalMinChisqVec[imin] < ( ChisqBF-ChisqBF*eps) ){
	ChisqBF = Dm2LocalMinChisqVec[imin];
	ChisqBF_Idx = imin;
	//cout << "(ChisqBF-ChisqBF*eps)/ChisqBF = " << (ChisqBF-ChisqBF*eps)/ChisqBF << endl;
      }
      
      // search for the lowest Dm2 
      if( Dm2LocalMinValVec[imin] < ( Dm2BF-eps ) ){
	Dm2BF = Dm2LocalMinValVec[imin];
	Dm2BF_Idx = imin;
	//cout << "(Dm2BF-Dm2BF*eps)/Dm2BF = " << (Dm2BF-Dm2BF*eps)/Dm2BF << endl;
      }

    } // loop dm2 minima

	
    //if(Dm2BF_Idx != ChisqBF_Idx){
    if( Dm2LocalMinValVec[ChisqBF_Idx] > Dm2Thr ){

      cout << endl;
      cout << endl;
      cout << "****************************" << endl;
      cout << "Toy " << itoy << endl;
      cout << "ChisqBF_Idx != Dm2BF_Idx" << endl;
      cout << "True Dm2 = " << input_params->Value(conv::kPhypDm23sqBar) << endl;
      cout << "First Min Dm2 " << Dm2LocalMinValVec[0] << endl;
      cout << "First ChisqBF " << Dm2LocalMinChisqVec[0] << endl;
      cout << "Min Dm2 = " << Dm2LocalMinValVec[ChisqBF_Idx] << endl;
      cout << "Min Chisq = " << Dm2LocalMinChisqVec[ChisqBF_Idx] << endl;
      cout << "Dm2BF_Idx = " << Dm2BF_Idx << endl;
      cout << "ChisqBF_Idx = " << ChisqBF_Idx << endl;
      cout << endl;

      hDm2Loc->Fill(Dm2LocalMinValVec[ChisqBF_Idx]);

      for(int imin=0;imin<Nmin;imin++){
	cout << Dm2LocalMin[imin] << " --> " 
	     << Dm2LocalMinChisqVec[imin] << " --> " 
	     << Dm2LocalMinValVec[imin] 
	     << endl;
      }
      
      NOtherMin++;
    }

  } // loop toys
  
  TCanvas *cDm2Loc = new TCanvas("cDm2Loc","cDm2Loc");
  hDm2Loc->Draw("histo");

  cout << endl;
  cout << endl;
  cout << "Prob. other minima = " << (double)NOtherMin / (double)NEntries * 100 << " %" << endl;
  cout << endl;
  
  
  return;
}
