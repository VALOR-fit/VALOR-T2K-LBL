void change_prompt()
{
    TRint* rint = dynamic_cast<TRint*>(gApplication);
    if (!rint) return;

    rint->SetPrompt("RooT2K [%d] ");
}

void root2k()
{   
    TString script_dir = gSystem->Getenv("T2KOSCANA");
    script_dir += "./";

    TString curr_dir = gSystem->pwd();

    gSystem->cd(script_dir.Data());

    gROOT->ProcessLine(".x $T2KOSCANA/root/loadincs.C");    
    gROOT->ProcessLine(".x $T2KOSCANA/root/loadlibs.C");    

    gSystem->cd(curr_dir.Data());

    change_prompt();
}
