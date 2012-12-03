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
  cout << endl << endl << "average chisq " << minchisq << " at " << avdchisq->GetXaxis()->GetBinCenter(minx) << " " << avdchisq->GetYaxis()->GetBinCenter(miny) << endl;

  //subtract the average (if required)
  if(ischisq) {
    for(int ix = 1; ix <= avdchisq->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= avdchisq->GetNbinsY(); iy++) {
	double current = avdchisq->GetBinContent(ix,iy);
	avdchisq->SetBinContent(ix,iy,current - minchisq);
      }//iy
    }//ix
  }//ischisq

  // Get external contours
  std::string line;
  int ip = 0;
  TGraph ext1_contour(100);
  std::ifstream ext1file("/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/inputs/external/contours/minos2012_10p71E20POT.data");
  while (std::getline(ext1file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    ext1_contour.SetPoint(ip, a, b);
    ip++;
  }
  ext1_contour.Set(ip);
  ip = 0;
  TGraph ext2_contour(100);
  std::ifstream ext2file("/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/inputs/external/contours/minos2012_10p71E20POT+37p9ktyr.data");
  while (std::getline(ext2file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    ext2_contour.SetPoint(ip, a, b);
    ip++;
  }
  ext2_contour.Set(ip);
  ip = 0;
  TGraph ext3_contour(100);
  std::ifstream ext3file("/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/inputs/external/contours/sk2012_runs1to4_LE2nu.data");
  while (std::getline(ext3file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    ext3_contour.SetPoint(ip, a, b);
    ip++;
  }
  ext3_contour.Set(ip);
  ip = 0;
  TGraph ext4_contour(100);
  std::ifstream ext4file("/Users/thomasdealtry/T2K/ralrepo/t2k/t2kana/nuosc/trunk/inputs/external/contours/sk2012_runs1to4_zenith3nu.data");
  while (std::getline(ext4file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    ext4_contour.SetPoint(ip, a, b);
    ip++;
  }
  ext4_contour.Set(ip);

  //and plot
  gStyle->SetOptStat(0);
  avdchisq->SetTitle("deltachisq surface");
  TCanvas *c1 = new TCanvas(); avdchisq->Draw("colz");  c1->SaveAs("sensitivitycolz.pdf");
  double constchisq[2] = {2.30, 4.61};
  avdchisq->SetContour(2, constchisq); // this also reduces the colz palette to only 2 colors!
  TCanvas *c2 = new TCanvas();
  avdchisq->SetTitle("");
  avdchisq->Draw("cont1");
  TFile fkit("of_constchi2_v3.root");
  TH2D *kitdchisq = 0;
  fkit.GetObject("h_dchi2",kitdchisq);
  kitdchisq->SetLineStyle(2);
  kitdchisq->SetContour(2, constchisq);
  kitdchisq->Draw("cont1same");
  ext1_contour.SetLineColor(kBlack);
  ext1_contour.SetLineStyle(1);
  ext2_contour.SetLineColor(kBlack);
  ext2_contour.SetLineStyle(3);
  ext3_contour.SetLineColor(kBlue);
  ext3_contour.SetLineStyle(1);
  ext4_contour.SetLineColor(kBlue);
  ext4_contour.SetLineStyle(3);
  ext1_contour.SetName("extcont1");
  ext2_contour.SetName("extcont2");
  ext3_contour.SetName("extcont3");
  ext4_contour.SetName("extcont4");
  ext1_contour.Draw("Csame");
  ext2_contour.Draw("Csame");
  ext3_contour.Draw("Csame");
  ext4_contour.Draw("Csame");
  TLegend *leg = new TLegend(0.0,0.8,0.25,1.0);
  leg->AddEntry(avdchisq,"VaLOR");
  leg->AddEntry(kitdchisq,"KIT");
  leg->AddEntry("ext1cont","MINOS beam 90%","l");
  leg->AddEntry("ext2cont","MINOS beam + atmospherics 90%","l");
  leg->AddEntry("ext3cont","SK atmospherics L/E 90%","l");
  leg->AddEntry("ext3cont","SK atmospherics zenith angle 90%","l");
  leg->Draw();
  c2->SaveAs("sensitivitycont.pdf");

  TFile fout("sensitivity.root","RECREATE");
  avdchisq->SetDirectory(0);
  avdchisq->Write();
  fout.Close();
}
