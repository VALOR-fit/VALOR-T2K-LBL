{
//TString opt = "SuperK+Beam+NeutrinoGenerator";
  TString opt = "SuperK";

  gSystem->Load("libGpad");  
  gSystem->Load("libPhysics");  
  gSystem->Load("libMatrix");   
  gSystem->Load("libHist");     
  gSystem->Load("libGraf");     
  gSystem->Load("libTree");     
  gSystem->Load("libRIO");   
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");   

  gSystem->Load("$T2KANA/lib/libT2KAnaYr1NuMuOscAnalysis.so");

  using namespace std;
  using namespace t2kana;
  using namespace t2kana::yr1numu;

  SuperkExpectation * skexp = new SuperkExpectation;
  skexp->Read("inputs/analysis_configuration.xml");
  
  if(opt.Contains("Beam")) {
    skexp->GetBeamInputs()->Draw(new TCanvas(), "SKFlux");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "ND280Flux");
    skexp->GetBeamInputs()->Draw(new TCanvas(), "OffAxisRatio");
  }
  if(opt.Contains("SuperK")) {
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "ReconstructionEfficiency");
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "SelectionEfficiency/OneMuRingTight");
    skexp->GetSuperkInputs()->Draw(new TCanvas(), "EnergyResolutionMatrix/OneMuRingTight");
  }
  if(opt.Contains("NeutrinoGenerator")) {
    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKSigTotal");
    skexp->GetGeneratorInputs()->Draw(new TCanvas(), "SKFracMode");
  }
}

