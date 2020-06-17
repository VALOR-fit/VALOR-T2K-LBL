/*
// ---------------------------------------------------------------------------------------------
// Test Mondrian binning specified in the input XML file.
// An example binning file is provided.
//
// Run as:
// % valor_dune 'test_binning_mondrian.C("example_binning-Ereco,yreco-mondrian_01.xml")'
//
// C.Andreopoulos
// ---------------------------------------------------------------------------------------------
*/

void test_binning_mondrian(const char * binning_file)
{
  using namespace std;
  using namespace valor;

  KBCubHcmbMondrian binning("test","A test binning scheme");
  binning.Read(binning_file);

  conv::KineSpace_t ksp = binning.KineSpace();
  cout << "> Kine space: " << util::AsString(ksp) << endl;

  unsigned int ndim = binning.NDimensions();
  cout << "> Number of dimensions: " << ndim << endl;

  unsigned int nbin = binning.NBins();
  cout << "> Number of bins: " << nbin << endl;

  for(unsigned int idim = 0; idim < ndim; idim++) {
     conv::KineVar_t kv = binning.KineVar(idim);
     cout << "\t Dimension: #" << idim << " -> variable: " << util::AsString(kv) << endl;

     double min = binning.Minimum(idim);
     double max = binning.Maximum(idim);
     cout << "\t\t Range: [" << min << ", " << max << "]" << endl;
  }

  cout << "> Bin edges: " << endl;

  for(unsigned int ibin = 0; ibin < nbin; ibin++) {
     cout << "\t Bin: " << ibin << endl;
     for(unsigned int idim = 0; idim < ndim; idim++) {
       conv::KineVar_t kv = binning.KineVar(idim);
       if(idim == 0) cout << "\t\t Edges: ";
       cout << util::AsString(kv) << " = [" 
            << binning.BinMinimum(ibin,idim) << ", " 
            << binning.BinMaximum(ibin,idim) << "] ";
       if(idim < ndim-1) cout << " - ";
     }
     cout << endl;
  }

  cout << "Testing bin finding: " << endl;
  const unsigned int np = 20;
  TRandom3 rand(0);
  for(unsigned int ip = 0; ip < np; ip++) {
     cout << "\t Random point = " << ip << endl;
     vector<double> xarr(ndim);
     for(unsigned int idim = 0; idim < ndim; idim++) {
        conv::KineVar_t kv = binning.KineVar(idim);
        double min = binning.Minimum(idim);
        double max = binning.Maximum(idim);
        double x = min + (max - min) * rand.Uniform();
        xarr[idim] = x;
        if(idim == 0) cout << "\t\t (";
        cout << util::AsString(kv) << " = " << x;
        if(idim < ndim-1) cout << ", ";        
     }
     cout << ") --> bin = ";
     unsigned int bin = binning.FindBin(xarr);
     cout << bin << endl;
  }

  // visualize binning in some usual casesø
  if(ndim == 2) {
     double xmin = binning.Minimum(0);
     double xmax = binning.Maximum(0);
     double ymin = binning.Minimum(1);
     double ymax = binning.Maximum(1);
     TCanvas * c = new TCanvas;     
     TH2F * h = (TH2F*)c->DrawFrame(xmin,ymin,xmax,ymax);
     h->GetXaxis()->SetTitle( util::AsString(binning.KineVar(0)).c_str() );
     h->GetYaxis()->SetTitle( util::AsString(binning.KineVar(1)).c_str() );
     h->Draw();
     for(unsigned int ibin = 0; ibin < nbin; ibin++) {
        double x0 = binning.BinMinimum(ibin,0);
        double x1 = binning.BinMaximum(ibin,0);
        double y0 = binning.BinMinimum(ibin,1);
        double y1 = binning.BinMaximum(ibin,1);
        TLine * lxl = new TLine(x0,y0,x1,y0);
        TLine * lxh = new TLine(x0,y1,x1,y1);
        TLine * lyl = new TLine(x0,y0,x0,y1);
        TLine * lyh = new TLine(x1,y0,x1,y1);
        lxl->Draw();
        lxh->Draw();
        lyl->Draw();
        lyh->Draw();
        TLatex * label = new TLatex(x0,y0,Form("%d",ibin));
        label->SetTextSize(0.05);
        label->SetTextColor(kRed);
        label->Draw();
     }
  }
}
