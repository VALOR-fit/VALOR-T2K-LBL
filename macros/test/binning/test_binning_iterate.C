/*
// ---------------------------------------------------------------------------------------------
// Test iterating over the Mondrian binning scheme
//
// F. Bench
// ---------------------------------------------------------------------------------------------
*/

void test_binning_iterate()
{
  using namespace valor;

  CompoundFitElements elems;
  elems.Read(util::CompleteConfigPath("analysis/t2k_analysis_configuration_oa2018_run1-9d_v1_postfit.xml").c_str(),"nomc nosplines");
  FitModelABC * modelabc = elems.FitModel();
  CompoundFitModel & model = *dynamic_cast<CompoundFitModel*>(modelabc);
  FitModelABC * submodelabc = model.GetSubModel(2);   // Get the first submodel with a 2D binning scheme (Ereco-thetareco and Etrue)
  GenericSingleSampleFitModel & submodel = *dynamic_cast<GenericSingleSampleFitModel*>(submodelabc);
  std::cout << "Using submodel " << submodel.Name() << std::endl;

  const SingleSampleKineBinning* sskb = submodel.GetBinning();
  KBCubHcmbMondrian  * binning_mond = dynamic_cast<KBCubHcmbMondrian*> (sskb->RecoOptBinning());

  // First test the Mondrian binning scheme, iterating over theta at fixed Ereco
  {
    std::cout << "Iterating over Mondrian thetareco at fixed Ereco = 0.075" << std::endl;
    std::cout << "Nbins (global) = " << binning_mond->NBins() << std::endl;
    conv::KineVar_t kv = conv::kvThetareco;
    std::map<conv::KineVar_t, double> other_kv_central_vals;
    other_kv_central_vals[conv::kvEreco] = 0.075;

    std::cout << "Nbins (local) = " << binning_mond->NBins(kv, other_kv_central_vals) << std::endl;

    int ibin_global = -1;
    int ix = 0;
    while(true) {
      ibin_global = binning_mond->GetNextBin(ibin_global, kv, other_kv_central_vals);
      if(ibin_global < 0) { break; }
      std::cout << "Found bin ix = " << ix << ", iglob = " << ibin_global << ": ["
                << binning_mond->BinMinimum(ibin_global, kv) << ", " << binning_mond->BinMaximum(ibin_global, kv)
                << "]" << std::endl;
      ix++;
    }
  }

  // Test the Mondrian binning scheme, iterating over Ereco at fixed theta
  {
    std::cout << "\n\nIterating over Mondrian Ereco at fixed thetareco = 45" << std::endl;
    conv::KineVar_t kv = conv::kvEreco;
    std::map<conv::KineVar_t, double> other_kv_central_vals;
    other_kv_central_vals[conv::kvThetareco] = 45;

    std::cout << "Nbins (local) = " << binning_mond->NBins(kv, other_kv_central_vals) << std::endl;

    int ibin_global = -1;
    int ix = 0;
    while(true) {
      ibin_global = binning_mond->GetNextBin(ibin_global, kv, other_kv_central_vals);
      if(ibin_global < 0) { break; }
      std::cout << "Found bin ix = " << ix << ", iglob = " << ibin_global << ": ["
                << binning_mond->BinMinimum(ibin_global, kv) << ", " << binning_mond->BinMaximum(ibin_global, kv)
                << "]" << std::endl;
      ix++;
    }
  }

  // Test the Mondrian binning scheme, iterating in 2D first over theta, then Ereco
  {
    std::cout << "\n\nIterating in 2D over Mondrian thetareco, Ereco" << std::endl;
    conv::KineVar_t kv_outer = conv::kvThetareco;
    conv::KineVar_t kv_inner = conv::kvEreco;
    std::map<conv::KineVar_t, double> other_kv_central_vals_outer;
    // TODO: So far this only works in 2D with rectangular binning
    // TODO: Modify this at the end of the Ereco loop?
    // TODO: Or need a loop around all this to iterate over starting points, but skipping any global bins that we've
    // TODO:   already visited? (otherwise weird shaped binnings may not work with this) (would need to keep tack of
    // TODO:   visited bins separately for each dimension, as there'll be some overlap)
    other_kv_central_vals_outer[kv_inner] = binning_mond->BinCentre(0, kv_inner);

    int ibin_global_outer = -1;
    int itheta = 0;
    while(true) {
      ibin_global_outer = binning_mond->GetNextBin(ibin_global_outer, kv_outer, other_kv_central_vals_outer);
      if(ibin_global_outer < 0) { break; }
      std::cout << "Found bin itheta = " << itheta << ", iglob = " << ibin_global_outer << ": ["
                << binning_mond->BinMinimum(ibin_global_outer, kv_outer) << ", " << binning_mond->BinMaximum(ibin_global_outer, kv_outer)
                << "]" << std::endl;

      // Loop over Ereco
      std::map<conv::KineVar_t, double> other_kv_central_vals_inner;
      other_kv_central_vals_inner[kv_outer] = binning_mond->BinCentre(ibin_global_outer, kv_outer);
      int ibin_global_inner = -1;
      int iereco = 0;
      while(true) {
        ibin_global_inner = binning_mond->GetNextBin(ibin_global_inner, kv_inner, other_kv_central_vals_inner);
        if(ibin_global_inner < 0) { break; }
        std::cout << "    Found bin iEreco = " << iereco << ", iglob = " << ibin_global_inner << ": ["
                  << binning_mond->BinMinimum(ibin_global_inner, kv_inner) << ", " << binning_mond->BinMaximum(ibin_global_inner, kv_inner)
                  << "]" << std::endl;
        iereco++;
      }
      itheta++;
    }
  }

}
