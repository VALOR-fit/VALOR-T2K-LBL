

const double conf68 = 2.30;
const double conf90 = 4.61;

void draw_marginalization_contour_2D(const char * input_file_pattern = "data.root")
{
  gStyle->SetPalette(56);
  
  gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(0.1);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.13);
  
  
  TFile * tf = new TFile(input_file_pattern);
  
  TH2D * h_chisq = new TH2D(*((TH2D*)tf->Get("hchisq_marg_2D")));
  h_chisq->SetTitle("");
  h_chisq->SetContour(100);
  
  TH2D * h_delta = new TH2D(*h_chisq);
  
  double minimum = 1E100;
  double xmin = 1E100;
  double ymin = 1E100;
  int min_bin_x = -1;
  int min_bin_y = -1;
  
  for (int ix = 1 ; ix <= h_chisq->GetXaxis()->GetNbins();ix++) {
    for (int iy = 1 ; iy <= h_chisq->GetYaxis()->GetNbins();iy++) {
      double content = h_chisq->GetBinContent(ix,iy);
      if (content < minimum) {
        min_bin_x = ix;
        min_bin_y = iy;
        xmin = h_chisq->GetXaxis()->GetBinCenter(ix);
        ymin = h_chisq->GetYaxis()->GetBinCenter(iy);
        minimum = content;
      }
    }
  }
  
  for (int ix = 1 ; ix <= h_chisq->GetXaxis()->GetNbins();ix++) {
    for (int iy = 1 ; iy <= h_chisq->GetYaxis()->GetNbins();iy++) {
      h_delta->SetBinContent( ix, iy, h_chisq->GetBinContent(ix,iy)-minimum );
    }
  }
  
  TMarker * bf = new TMarker(xmin,ymin,kFullStar);
  bf->SetMarkerColor(kViolet);
  bf->SetMarkerSize(1.5);
  
  TCanvas * c_chisq = new TCanvas();
  c_chisq->cd();
  
  h_chisq->Draw("colz");
  bf->Draw("same");
  c_chisq->Update();
  
  TCanvas * c_delta = new TCanvas();
  c_delta->cd();
  
  h_delta->Draw("colz");
  bf->Draw("same");
  c_delta->Update();
  
  TH2D * h_contour = new TH2D(*h_delta);
  h_contour->SetContour(2);
  h_contour->SetContourLevel(0,conf68);
  h_contour->SetContourLevel(1,conf90);
  TCanvas * c_contour = new TCanvas();
  c_contour->cd();
  
  h_contour->Draw("cont1");
  bf->Draw("same");
  
}
