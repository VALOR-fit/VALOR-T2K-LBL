/*
// ---------------------------------------------------------------------------------------------
// Test edge-to-edge binning specified in the input XML file.
// Several example binning files are provided. 
//
// Run as:
// % valor_dune 'test_binning_e2e.C("example_binning-Ereco,yreco-e2e_01.xml")'
// % valor_dune 'test_binning_e2e.C("example_binning-Etrue,Q2true-e2e_01.xml")'
// % valor_dune 'test_binning_e2e.C("example_binning-Etrue,Q2true-e2e_02.xml")'
//
// C.Andreopoulos
// ---------------------------------------------------------------------------------------------
*/

void test_binning_e2e(const char * binning_file)
{
  using namespace std;
  using namespace valor;

  KBCubHcmbEdge2Edge binning("test","A test binning scheme");
  binning.Read(binning_file);

  conv::KineSpace_t ksp = binning.KineSpace();
  cout << "> Kine space: " << util::AsString(ksp) << endl;

  unsigned int ndim = binning.NDimensions();
  cout << "> Number of dimensions: " << ndim << endl;

  unsigned int nbin = binning.NBins();
  cout << "> Number of multi-dimensional bins: " << nbin << endl;

  for(unsigned int idim = 0; idim < ndim; idim++) {
     conv::KineVar_t kv = binning.KineVar(idim);
     cout << "\t Dimension: #" << idim << " -> variable: " << util::AsString(kv) << endl;

     unsigned int nbin_dim = binning.NBins(idim);
     cout << "\t\t NBins: " << nbin_dim << endl;

     bool fixed_width = binning.FixedWidth(idim);
     cout << "\t\t Fixed width? " << ((fixed_width) ? "true" : "false") << endl;

     double min = binning.Minimum(idim);
     double max = binning.Maximum(idim);
     cout << "\t\t Range: [" << min << ", " << max << "]" << endl;

     const double * bin_edges = binning.BinEdgeArr(idim);
     for(unsigned int i = 0; i <= nbin_dim; i++) {
       if(i == 0) cout << "\t\t Bin edges: [";
       cout << bin_edges[i];
       if(i < nbin_dim) cout << ", ";
       if(i == nbin_dim) cout << "]" << endl;
     }

     unsigned int idim_test = binning.Dimension(kv);
     assert(idim == idim_test);
     unsigned int nbin_dim_test = binning.NBins(kv);
     assert(nbin_dim == nbin_dim_test);
  }

  cout << "Testing bin finding: " << endl;
  const unsigned int np = 20;
  TRandom3 rand(0);
  for(unsigned int ip = 0; ip < np; ip++) {
     cout << "\t Random point = " << ip << endl;
     vector<unsigned int> binarr(ndim);
     for(unsigned int idim = 0; idim < ndim; idim++) {
        conv::KineVar_t kv = binning.KineVar(idim);
        double min = binning.Minimum(idim);
        double max = binning.Maximum(idim);
        double x = min + (max - min) * rand.Uniform();
        unsigned int bin = binning.FindBin(x,idim);
        binarr[idim] = bin;
        cout << "\t\t " << util::AsString(kv) << " = " << x << " --> FindBin returns: " << bin << endl;
     }
     unsigned int linbin = binning.LinearizedBin(binarr);
     cout << "\t\t LinearizedBin returns: " << linbin << endl;
     cout << "\t\t Inverting LinearizedBin returns: (";
     for(unsigned int idim = 0; idim < ndim; idim++) {
         cout << binning.Bin(linbin,idim);
         if(idim < ndim - 1) cout << ", ";
     }
     cout << ")" << endl;
          
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
