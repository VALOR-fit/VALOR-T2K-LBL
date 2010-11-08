{
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

  SuperkExpectation skexp;
  skexp.GetBeamInputs()->Read("inputs/analysis_configuration.xml");

  TCanvas * cskflux = new TCanvas("cskflux","",50,50,400,400);

  skexp.GetBeamInputs()->Draw(cskflux, "SKFlux");
}
