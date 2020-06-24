{
  using namespace valor;

  int zone_color[] = { kRed-9, kBlue-10, kGreen-10, kMagenta-7, kYellow+1 };

  // sample
  dune::conv::Sample_t sample = dune::conv::ksNDNuMuCC1Trk;

  // overall range
  const double kemin = 0;
  const double kemax = 50;
  const double kymin = 0;
  const double kymax = 1;

  TCanvas * c = new TCanvas;
  TH2F * hframe = (TH2F*) c->DrawFrame(kemin,kymin,kemax,kymax);
  hframe->Draw();

  TLine * l = 0;
  TBox * b = 0;

  // loop over binning zones for the current sample
  for(int izone = 0; izone < dune::conv::kNMaxBinZones; izone++) {

     // for current binning zone, get the energy and y bins (variable size)
     const int      ne = dune::util::NumErecoBins  (sample, izone);
     const double * eb = dune::util::ErecoBinEdges (sample, izone);
     const int      ny = dune::util::NumYrecoBins  (sample, izone);
     const double * yb = dune::util::YrecoBinEdges (sample, izone);         
     
     if(yb==0 || eb==0) continue;

     // work out the energy and y range of the current zone
     double emin = eb[TMath::LocMin(ne+1,eb)];
     double emax = eb[TMath::LocMax(ne+1,eb)];
     double ymin = yb[TMath::LocMin(ny+1,yb)];
     double ymax = yb[TMath::LocMax(ny+1,yb)];

     cout << "** zone = " << izone << endl;
     cout << "num of e bins: " << ne << endl;
     cout << "num of y bins: " << ny << endl;
     cout << "ranges: e = [" << emin << ", " << emax << "], y = [" << ymin << ", " << ymax << "]" << endl;
 
     b = new TBox(emin,ymin,emax,ymax);
     b->SetFillColor(zone_color[izone]);
     b->Draw();

     // draw horizontal and vertical lines at the bin boundaries
     for(int i = 0; i < ne+1; i++) {
       l = new  TLine(eb[i],ymin,eb[i],ymax);
       l->Draw();
     }
     for(int i = 0; i < ny+1; i++) {
       l = new  TLine(emin,yb[i],emax,yb[i]);
       l->Draw();
     }
  }//zone

}
