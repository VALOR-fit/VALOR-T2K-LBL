//
// Take a lot of chisq surfaces, average them, and plot the contours. Also plot KIT contour
// Its a bit of a hack at the moment (not general), but I'll clean it up later
//
// Takes runs from e.g.
//
// $T2KOSCANA/production/submit_toyexpt_fitting_job_fitsyst_minuit2 --resource RAL --softw-dir $T2KOSCANA --job-dir $PWD --cpu-time 48:00:00 --tweak-systematics 1 --nruns 100 --nexpts 1 --sinsq2theta12 -99 --sinsq2theta13 -99 --sinsq2theta23 0.957 --deltacp -99 --dm12sq -99 --dm23sq 0.00239 --fit-param-list nuosc3f_numudisapp_fix4+nuisanceselection --config-file analysis_configuration.xml --tag sensitivity --stat-fluctuations --store-chisq-withsyst-surf
//
// Tom Dealtry
//

void sensitivity_from_many_chisq_surfaces
(
 const char * fileprefix = "sensitivity-run",               // prefix of the root file where the surfaces are located
 const char * filesuffix = "systfluc_All-statfluc_yes-fitparam_nuosc3f_numudisapp_fix4+nuisanceselection-sinsq2theta[12_0.870,13_0.100,23_0.957],dcp_0.000,dmsq[12_0.00007600,23_0.002390]eV2.root",
 const char * surfacename   = "chisq_surface_syst_1Rmu_0",  // name of the dchisq surfaces
 int          nfiles        = 100,                          // number of files to access
 bool         ischisq       = true                          // Input is chisq if true; dchisq if false
)
{
  //need to initialise avdchisq outside the loop, or segfaults ensue
  const char * filename1 = Form("%s_%d-%s", fileprefix, 0, filesuffix);
  TFile f1(filename1);
  TH2D * dchisq1 = 0;
  f1.GetObject(surfacename,dchisq1);
  if(!dchisq1) { f1.Close(); return; }
  TH2D * avdchisq = new TH2D(dchisq1->GetTitle(), dchisq1->GetName(), dchisq1->GetNbinsX(), dchisq1->GetXaxis()->GetXmin(), dchisq1->GetXaxis()->GetXmax(), dchisq1->GetNbinsY(), dchisq1->GetYaxis()->GetXmin(), dchisq1->GetYaxis()->GetXmax());

  //first add everything up
  int nfilesexist = 0;
  for(int ifile = 0; ifile < nfiles; ifile++) {
    const char * filename = Form("%s_%d-%s", fileprefix, ifile, filesuffix);
    TFile f(filename);
    TH2D * dchisq = 0;
    f.GetObject(surfacename,dchisq);
    if(!dchisq) { f.Close(); continue;}
    nfilesexist++;
    for(int ix = 1; ix <= dchisq->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= dchisq->GetNbinsY(); iy++) {
	double current = avdchisq->GetBinContent(ix,iy);
	double extra   =   dchisq->GetBinContent(ix,iy);
	avdchisq->SetBinContent(ix,iy,current+extra);	
      }//iy
    }//ix
    f.Close();
  }//ifile

  //now do the average
  for(int ix = 1; ix <= avdchisq->GetNbinsX(); ix++) {
    for(int iy = 1; iy <= avdchisq->GetNbinsY(); iy++) {
      double current = avdchisq->GetBinContent(ix,iy);
      avdchisq->SetBinContent(ix,iy,current/(double)nfilesexist);
    }//iy
  }//ix

  //find the mimumum of the average
  double minchisq = avdchisq->GetMinimum();
  int minbin = avdchisq->GetMinimumBin();
  int minx, miny, minz;
  avdchisq->GetBinXYZ(minbin, minx, miny, minz);
  cout << "average chisq " << minchisq << " at " << avdchisq->GetXaxis()->GetBinCenter(minx) << " " << avdchisq->GetYaxis()->GetBinCenter(miny) << endl;

  //subtract the average (if required)
  if(ischisq) {
    for(int ix = 1; ix <= avdchisq->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= avdchisq->GetNbinsY(); iy++) {
	double current = avdchisq->GetBinContent(ix,iy);
	avdchisq->SetBinContent(ix,iy,current - minchisq);
      }//iy
    }//ix
  }//ischisq

  //and plot
  double constchisq[2] = {2.30, 4.61};
  avdchisq->SetContour(2, constchisq);
  TCanvas *c1 = new TCanvas(); avdchisq->Draw("colz");  c1->SaveAs("sensitivitycolz.pdf");
  TCanvas *c2 = new TCanvas();
  avdchisq->Draw("cont1");
  TFile fkit("of_constchi2_v3.root");
  TH2D *kit68 = 0, *kit90=0;
  fkit.GetObject("h_cont_68cl",kit68);
  fkit.GetObject("h_cont_90cl",kit90);
  kit68->SetLineStyle(2);
  kit68->SetLineColor(kBlack);
  kit90->SetLineStyle(2);
  kit68->Draw("cont1same");
  kit90->Draw("cont1same");
  c2->SaveAs("sensitivitycont.pdf");

  TFile fout("sensitivity.root","RECREATE");
  avdchisq->SetDirectory(0);
  avdchisq->Write();
  fout.Close();
}
