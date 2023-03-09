#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

void validate_mtuples_OA2021_hard_cuts(TString in_file_name) {
	TFile *in_file = new TFile(in_file_name);
	TTree *mtuple = (TTree *) in_file->Get("mtuple");

	mtuple->SetBranchStatus("*", 1);

	double re_cuts = mtuple->GetEntries("(fqevis<100 || erec>1250) && iclass==11");
	double multire_cuts = mtuple->GetEntries("(fqevis<100 || erec>1250) && iclass==14");
	double rmu_cuts = mtuple->GetEntries("fqmomm<200 && iclass==13");
	double multirmu_cuts_1de = mtuple->GetEntries("(evisible<30 || (fqmreloss[ipp][0] + fqmreloss[ipp][1]) < 300) && iclass==31");
	double multirmu_cuts_2de = mtuple->GetEntries("evisible<30 && iclass==32");

	std::cout << "For : " << in_file_name << " events outside of sample range are : " << std::endl;
	std::cout << "    1Re            : " << re_cuts << std::endl;
	std::cout << "    1Re1de         : " << multire_cuts << std::endl;
	std::cout << "    1Rmu           : " << rmu_cuts << std::endl;
	std::cout << "    numucc1pi, 1de : " << multirmu_cuts_1de << std::endl;
	std::cout << "    numucc1pi, 2de : " << multirmu_cuts_2de << std::endl;

	if (re_cuts == 0 && re_cuts == 0 && rmu_cuts == 0 && multirmu_cuts_1de == 0 && multirmu_cuts_2de == 0){}
	else {std::cout << "WARNING! NON ZERO VALUE DETECTED!" << std::endl;}

	delete in_file;
}
