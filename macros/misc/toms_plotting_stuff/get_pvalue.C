{
  int ntoys = 0;
  int counter = 0;
  for(int i = 0; i < 50; i++) {
    TFile f(Form("goodness-of-fit-123good2Nsintheta-%d-run_3.root",i));
    if(f.IsZombie()) continue;
    goodness_params->Print();
    counter += goodness_params(2);
    ntoys   += goodness_params(3);
    cout << counter << "\t" << ntoys << endl
	 << "gives p-value of " << (double)counter/(double)ntoys << endl;
  }
}
