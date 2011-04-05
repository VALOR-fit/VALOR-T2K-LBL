{
  TString libs0 = gSystem->GetDynamicPath();
  TString libs  = libs0 + ":/usr/lib:/usr/local/lib:/opt/lib:/opt/local/lib";
  gSystem->SetDynamicPath(libs.Data());

  gSystem->Load("libGpad");
  gSystem->Load("libPhysics");
  gSystem->Load("libMatrix");
  gSystem->Load("libHist");
  gSystem->Load("libGraf");
  gSystem->Load("libTree");
  gSystem->Load("libRIO");
  gSystem->Load("libXMLIO");
  gSystem->Load("libMinuit");

  gSystem->Load("libT2KNuOsc.so");
////  gSystem->Load("$T2KOSCANA/fitter/libT2KNuOsc.so");
}

