#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <stdlib.h>

void ApplyCutsToMinituples(TString in_file_name, TString out_file_name) {
    TFile *in_file = new TFile(in_file_name);
    if (in_file->IsZombie()) { std::cerr << "In file doesn't exist"; std::abort(); }
    TTree *mtuple = nullptr;
    if (in_file->Get("mtuple")) {
        mtuple = (TTree *) in_file->Get("mtuple");
    } else {
        std::cerr << "Error: TTree mtuple does not exist in input file." << std::endl;
        exit(0);
    }


    TFile *out_file = new TFile(out_file_name, "RECREATE");
    TTree *out_tree = mtuple->CloneTree(0);

    // Get the norm histogram from the input file and copy it to the output file
    TH1D *norm = (TH1D*) in_file->Get("norm");
    if (norm) {
        norm->Clone()->Write();
    } else {
        std::cerr << "Error: Norm TH1D does not exist in input file." << std::endl;
        exit(0);
    }

    Int_t iclass, ipp;
    Float_t fqevis, fqmreloss[200][6], evisible, fqmomm;
    Double_t erec;

    //Set Branch Addresses
    mtuple->SetBranchAddress("iclass", &iclass);
    mtuple->SetBranchAddress("erec", &erec);
    mtuple->SetBranchAddress("fqevis", &fqevis);
    mtuple->SetBranchAddress("evisible", &evisible);
    mtuple->SetBranchAddress("fqmomm", &fqmomm);
    mtuple->SetBranchAddress("fqmreloss", fqmreloss);
    mtuple->SetBranchAddress("ipp", &ipp);
    //copying+cutting code here
    double re_cuts=0, multire_cuts=0, rmu_cuts=0, multirmu_cuts_1de=0, multirmu_cuts_2de = 0, under_limits = 0, beyond_limits = 0;

    for(int ientry; ientry<mtuple->GetEntries(); ientry++) {
        mtuple->GetEntry(ientry);
        if ((fqevis < 100 || erec > 1250) && iclass == 11) {
            re_cuts++;
          //  std::cout << "Discarding event " << ientry << " due to 1Re cut." << std::endl;
        } else if ((fqevis < 100 || erec > 1250) && iclass == 14) {
            multire_cuts++;
           // std::cout << "Discarding event " << ientry << " due to multire cut." << std::endl;
        } else if (fqmomm < 200 && iclass == 13) {
            rmu_cuts++;
           // std::cout << "Discarding event " << ientry << " due to 1Rmu cut." << std::endl;
        } else if ((evisible < 30 || (fqmreloss[ipp][0] + fqmreloss[ipp][1]) < 300) && iclass == 31) {
            multirmu_cuts_1de++;
           // std::cout << "Discarding event " << ientry << " due to multirmu_cuts_1de cut." << std::endl;
        } else if (evisible < 30 && iclass == 32) {
            multirmu_cuts_2de++;
            // std::cout << "Discarding event " << ientry << " due to multirmu_cuts_2de cut." << std::endl;
        } else if (erec<0){
            under_limits++;
        } else if (erec>30000){
            beyond_limits++;
        }
         else {
            out_tree->Fill();
        }
    }

    std::cout << "For : " << out_file_name << " events cut in each sample are : " << std::endl;
    std::cout << "    1Re            : " << re_cuts << std::endl;
    std::cout << "    1Re1de         : " << multire_cuts << std::endl;
    std::cout << "    1Rmu           : " << rmu_cuts << std::endl;
    std::cout << "    numucc1pi, 1de : " << multirmu_cuts_1de << std::endl;
    std::cout << "    numucc1pi, 2de : " << multirmu_cuts_2de << std::endl;
    std::cout << "    <0 Gev : " << under_limits << std::endl;
    std::cout << "    >30 GeV : " << beyond_limits << std::endl;

    double event_loss = mtuple->GetEntries() - out_tree->GetEntries();
    if (event_loss != (re_cuts + multire_cuts + rmu_cuts + multirmu_cuts_1de + multirmu_cuts_2de + under_limits + beyond_limits)){
        std::cerr << "ERROR: Events cut vs event difference does not match!" << std::endl;
    }

    out_tree->SetName("mtuple");
    out_tree->Write("mtuple", TObject::kOverwrite);
    in_file->Close();
    out_file->Close();
}