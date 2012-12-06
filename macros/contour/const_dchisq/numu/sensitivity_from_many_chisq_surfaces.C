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

#include <TFile.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>

void sensitivity_from_many_chisq_surfaces
(
 const char * fileprefix         = "minos",                      // prefix of the root file where the surfaces are located
 const char * kitfilename        = "of_constchi2_v3.root",       // location of a kit file, to compare with
 const char * kitsurfacename     = "h_dchi2",                    // name of the surface in the KIT file
 const char * kitlegendname      = "KIT",                        // name of the 'kit' contour in the legend (i.e. it can be non KIT!)
 const char * surfacenameprefix  = "chisqwithsyst",              // name of the dchisq surfaces
 int          nfiles             = 1000,                         // number of files to access
 bool         ischisq            = true                          // Input is chisq if true; dchisq if false
)
{
  //need to initialise avdchisq outside the loop, or segfaults ensue
  const char * filename1 = Form("%s-run_%d-expts_%d_%d.root", fileprefix, 0, 0, 0);
  TFile f1(filename1);
  TH2D * dchisq1 = 0;
  f1.GetObject(Form("%s_%d",surfacenameprefix,0),dchisq1);
  if(!dchisq1) { f1.Close(); return; }
  TH2D * avdchisq = new TH2D(dchisq1->GetTitle(), dchisq1->GetName(), dchisq1->GetNbinsX(), dchisq1->GetXaxis()->GetXmin(), dchisq1->GetXaxis()->GetXmax(), dchisq1->GetNbinsY(), dchisq1->GetYaxis()->GetXmin(), dchisq1->GetYaxis()->GetXmax());

  //first add everything up
  int nfilesexist = 0;
  for(int ifile = 0; ifile < nfiles; ifile++) {
    const char * filename = Form("%s-run_%d-expts_%d_%d.root", fileprefix, ifile, ifile, ifile);
    TFile f(filename);
    if(f.IsZombie()) continue;
    TH2D * dchisq = 0;
    f.GetObject(Form("%s_%d",surfacenameprefix,ifile),dchisq);
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
  cout << "Average of " << nfilesexist << " files" << endl;

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
  std::ifstream ext1file(Form("%s/inputs/contours/minos2012_10p71E20POT_90CL.data",gSystem->Getenv("T2KOSCANA")));
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
  std::ifstream ext2file(Form("%s/inputs/contours/minos2012_10p71E20POT+37p9ktyr_90CL.data",gSystem->Getenv("T2KOSCANA")));
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
  std::ifstream ext3file(Form("%s/inputs/contours/sk2012_runs1to4_LE2nu_90CL.data",gSystem->Getenv("T2KOSCANA")));
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
  std::ifstream ext4file(Form("%s/inputs/contours/sk2012_runs1to4_zenith3nu_90CL.data",gSystem->Getenv("T2KOSCANA")));
  while (std::getline(ext4file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    ext4_contour.SetPoint(ip, a, b);
    ip++;
  }
  ext4_contour.Set(ip);
  ip = 0;
  TGraph t2k_I_68_contour(100);
  std::ifstream t2k_I_68file(Form("%s/inputs/contours/t2k2010_run1_68CL_statsyst.data",gSystem->Getenv("T2KOSCANA")));
  while (std::getline(t2k_I_68file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    t2k_I_68_contour.SetPoint(ip, a, b);
    ip++;
  }
  t2k_I_68_contour.Set(ip);
  ip = 0;
  TGraph t2k_I_90_contour(100);
  std::ifstream t2k_I_90file(Form("%s/inputs/contours/t2k2010_run1_90CL_statsyst.data",gSystem->Getenv("T2KOSCANA")));
  while (std::getline(t2k_I_90file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    t2k_I_90_contour.SetPoint(ip, a, b);
    ip++;
  }
  t2k_I_90_contour.Set(ip);
  ip = 0;
  TGraph t2k_I_II_68_contour(100);
  std::ifstream t2k_I_II_68file(Form("%s/inputs/contours/t2k2010_run1+2_68CL_statsyst.data",gSystem->Getenv("T2KOSCANA")));
  while (std::getline(t2k_I_II_68file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    t2k_I_II_68_contour.SetPoint(ip, a, b);
    ip++;
  }
  t2k_I_II_68_contour.Set(ip);
  ip = 0;
  TGraph t2k_I_II_90_contour(100);
  std::ifstream t2k_I_II_90file(Form("%s/inputs/contours/t2k2010_run1+2_90CL_statsyst.data",gSystem->Getenv("T2KOSCANA")));
  while (std::getline(t2k_I_II_90file, line)) {
    std::istringstream iss(line);
    double a, b;
    if (!(iss >> a >> b)) { continue; } // header
    t2k_I_II_90_contour.SetPoint(ip, a, b);
    ip++;
  }
  t2k_I_II_90_contour.Set(ip);

  //set up external/old contour stydles
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
  t2k_I_68_contour.SetLineColor(kOrange+3);
  t2k_I_68_contour.SetLineStyle(1);
  t2k_I_90_contour.SetLineColor(kOrange+3);
  t2k_I_90_contour.SetLineStyle(3);
  t2k_I_II_68_contour.SetLineColor(kOrange+7);
  t2k_I_II_68_contour.SetLineStyle(1);
  t2k_I_II_90_contour.SetLineColor(kOrange+7);
  t2k_I_II_90_contour.SetLineStyle(3);
  t2k_I_68_contour.SetName("t2k_I_68cont");
  t2k_I_90_contour.SetName("t2k_I_90cont");
  t2k_I_II_68_contour.SetName("t2k_I_II_68cont");
  t2k_I_II_90_contour.SetName("t2k_I_II_90cont");

  //and plot
  gStyle->SetOptStat(0);
  avdchisq->SetTitle("deltachisq surface");
  //dchisq surface
  TCanvas *c1 = new TCanvas(); avdchisq->Draw("colz");  c1->SaveAs("sensitivitycolz.pdf");
  TFile fout("sensitivity.root","RECREATE");
  avdchisq->SetDirectory(0);
  avdchisq->Write("dchisq_for_colz");
  double constchisq90[1] = {4.61};
  avdchisq->SetContour(1, constchisq90);
  avdchisq->Write("dchisq_for_cont1_90");
  double constchisq68[1] = {2.30};
  avdchisq->SetContour(1, constchisq68);
  avdchisq->Write("dchisq_for_cont1_68");
  double constchisq[2] = {2.30, 4.61};
  avdchisq->SetContour(2, constchisq); // this also reduces the colz palette to only 2 colors!
  avdchisq->Write("dchisq_for_cont1");
  //contour
  TCanvas *c2 = new TCanvas();
  avdchisq->SetTitle("");
  avdchisq->Draw("cont1");
  TLegend *leg = new TLegend(0.0,0.8,0.25,1.0);
  leg->SetFillColor(0);
  leg->AddEntry(avdchisq,"VaLOR","l");
  TH2D *kitdchisq = 0;
  TFile fkit(kitfilename);
  if(!fkit.IsZombie()) {
    fkit.GetObject(kitsurfacename,kitdchisq);
    kitdchisq->SetLineStyle(2);
    kitdchisq->SetContour(2, constchisq);
    kitdchisq->Draw("cont1same");
    leg->AddEntry(kitdchisq,kitlegendname,"l");
    leg->Draw();
    c2->SaveAs("sensitivitycont_comparisononly.pdf");
  }
  //add externals (minos/sk)
  ext1_contour.Draw("Csame");
  ext2_contour.Draw("Csame");
  ext3_contour.Draw("Csame");
  ext4_contour.Draw("Csame");
  leg->AddEntry("extcont1","MINOS beam 90%","l");
  leg->AddEntry("extcont2","MINOS beam + atmospherics 90%","l");
  leg->AddEntry("extcont3","SK atmospherics L/E 90%","l");
  leg->AddEntry("extcont4","SK atmospherics zenith angle 90%","l");
  leg->Draw();
  c2->SaveAs("sensitivitycont.pdf");
  leg->Clear();
  //another, without externals but with old t2k result
  avdchisq->Draw("cont1");
  leg->AddEntry(avdchisq,"VaLOR","l");
  t2k_I_II_68_contour.Draw("Csame");
  t2k_I_II_90_contour.Draw("Csame");
  leg->AddEntry("t2k_I_II_68cont","T2K run I-II 68%","l");
  leg->AddEntry("t2k_I_II_90cont","T2K run I-II 90%","l");
  leg->Draw();
  c2->SaveAs("sensitivitycont_witholdt2k.pdf");
  //now add t2k run I
  t2k_I_90_contour.SetTitle("");
  t2k_I_90_contour.Draw("LA");
  t2k_I_68_contour.Draw("Lsame");
  t2k_I_II_68_contour.Draw("Csame");
  t2k_I_II_90_contour.Draw("Csame");
  avdchisq->Draw("cont1same");
  leg->AddEntry("t2k_I_68cont","T2K run I 68%","l");
  leg->AddEntry("t2k_I_90cont","T2K run I 90%","l");
  leg->Draw();
  c2->SaveAs("sensitivitycont_witht2kI.pdf");

  fout.Close();
}
