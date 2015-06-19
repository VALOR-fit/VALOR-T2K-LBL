void dothework(int ican, TCanvas** c) {
 c[0]->cd(ican); bestfit_params->Draw("chisq","status>2.99");
 c[1]->cd(ican); bestfit_params->Draw("chisqp","status>2.99");
 c[2]->cd(ican); bestfit_params->Draw("status","status<2.99");
 c[3]->cd(ican); pullalt_a1->Draw();
 c[4]->cd(ican); pullalt_a2->Draw();
 c[5]->cd(ican); pullalt_a3->Draw();
 c[6]->cd(ican); pullalt_a4->Draw();
 c[7]->cd(ican); pull_mean->Draw();
 c[8]->cd(ican); pull_dip->Draw();
 c[9]->cd(ican); meansum->Draw();
 c[10]->cd(ican); sigmasum->Draw();
 c[11]->cd(ican); bestfit_params->Draw("m","status>2.99");
 c[12]->cd(ican); bestfit_params->Draw("d","status>2.99");
 c[13]->cd(ican); bestfit_params->Draw("a1","status>2.99");
 c[14]->cd(ican); bestfit_params->Draw("a2","status>2.99");
 c[15]->cd(ican); bestfit_params->Draw("a3","status>2.99");
 c[16]->cd(ican); bestfit_params->Draw("a4","status>2.99");
 c[17]->cd(ican); bestfit_params->Draw("a1err","status>2.99");
 c[18]->cd(ican); bestfit_params->Draw("a2err","status>2.99");
 c[19]->cd(ican); bestfit_params->Draw("a3err","status>2.99");
 c[20]->cd(ican); bestfit_params->Draw("a4err","status>2.99");
 c[21]->cd(ican); bestfit_params->Draw("merr","status>2.99");
 c[22]->cd(ican); bestfit_params->Draw("derr","status>2.99");
}

int plot() {
  TFile *file[6];
  file[0] = TFile::Open("fit5_chisqcor_fixnotfit/out.root");
  file[1] = TFile::Open("fit5_chisqcor_movenotfit/out.root");
  file[2] = TFile::Open("fit5_chisqun_fixnotfit/out.root");
  file[3] = TFile::Open("fit5_chisqun_movenotfit/out.root");
  file[4] = TFile::Open("fit6_chisqcor/out.root");
  file[5] = TFile::Open("fit6_chisqun/out.root");

  const int nplots = 23;
  TCanvas *c[nplots];
  for(int i = 0; i < nplots; i++) {
    c[i] = new TCanvas();
    c[i]->Divide(3,2);
  }

  for(int i = 0; i < 6; i++) {
    file[i]->cd();
    dothework(i+1,c);
    //dothework(i+1,c,cp,cs);
  }
}
