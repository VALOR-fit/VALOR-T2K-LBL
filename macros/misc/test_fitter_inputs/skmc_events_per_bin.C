//take the "ori" templates (raw events per bin, not even POT scaling) & plot number of events per bin
//for each template, and for each MC sample
void skmc_events_per_bin
(
 const char * filename = "2013results/SPECTRA/mcstaterror/numu/Template_mcstaterr_84.root"
)
{
  using namespace t2kana;
  TFile f(filename);
  TCanvas * c = new TCanvas();
  c->SaveAs("events_per_template_bin.pdf[");
  const int ntemplates = 32;
  TH1I * h[ntemplates];
  TH2D * hin = 0;
  //loop over templates
  for(int it = 0; it < ntemplates; it++) {
    f.GetObject(Form("nsk_ori_0_%d",it),hin);
    TString test = hin->GetTitle();
    if(test.Contains("oscillated numu")) continue;
    h[it] = new TH1I(Form("h%d",it),"",10000,0,10000);
    double meanme = 0;
    int counter = 0;
    //loop over bins
    for(int ix = 1; ix <= hin->GetNbinsX(); ix++) {
      for(int iy = 1; iy <= hin->GetNbinsY(); iy++) {
	counter++;
	//if(hin->GetXaxis()->GetBinCenter(ix) > (3.0 + conv::kOneMeV)) continue;
	//if(hin->GetYaxis()->GetBinCenter(iy) > (3.0 + conv::kOneMeV)) continue;
	//if(hin->GetXaxis()->GetBinCenter(ix) < (0.3 - conv::kOneMeV)) continue;
	//if(hin->GetYaxis()->GetBinCenter(iy) < (0.3 - conv::kOneMeV)) continue;
	double content = hin->GetBinContent(ix,iy);
	if(content != 0) {
	  //cout << content << endl;
	  meanme += content;
	}
	h[it]->Fill(content);
      }//iy 
    }//ix 
    double ntotal = h[it]->Integral();
    double nempty = h[it]->GetBinContent(1);
    double emptyfrac = 100 * nempty / ntotal;
    double mean   = h[it]->GetMean();
    double rms    = h[it]->GetRMS();
    meanme /= (ntotal - nempty);
    h[it]->GetXaxis()->SetRange(2,50);
    h[it]->Draw();
    TLegend * l = new TLegend(0.3,0.4,0.9,0.9);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->AddEntry((TObject*)0,hin->GetTitle(),"");
    l->AddEntry((TObject*)0,Form("%f bins in range",ntotal),"");
    l->AddEntry((TObject*)0,Form("%d bins in range",counter),"");
    l->AddEntry((TObject*)0,Form("%f empty bins",nempty),"");
    l->AddEntry((TObject*)0,Form("(%f\%)",emptyfrac),"");
    l->AddEntry((TObject*)0,Form("Mean: %f",mean),"");
    l->AddEntry((TObject*)0,Form("RMS:  %f",rms),"");
    l->AddEntry((TObject*)0,Form("Mean (neglect 0 bin): %f",meanme),"");
    l->Draw();
    c->SaveAs("events_per_template_bin.pdf");
  }//it
  c->SaveAs("events_per_template_bin.pdf]");
}

