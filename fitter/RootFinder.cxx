//____________________________________________________________________________
/*
 
 Copyright (c) 2010-2014, The VALOR group. All rights reserved.
 
 ChangeLog:
 
 */
//____________________________________________________________________________
//#define __debug_output__
//#define __verbose_debug_output__
#include <fstream>

#include <TROOT.h>
#include <TGraph.h>
#include <TF1.h>
#include <TVector.h>

#include "RootFinder.h"

using std::cout;
using std::endl;
using std::ios;
using std::ifstream;

using namespace valor;

RootFinder::RootFinder()
{
    this->Init();
}
//___________________________________________________________________________
RootFinder::~RootFinder(void)
{
    this->CleanUp();
}
//___________________________________________________________________________
vector<double> RootFinder::GetRoots(void)
{
    
    //check if range is set
    if(GlobalXmax == -999 || GlobalXmin == -999) this->SetRange();
    
    double MinX = GlobalXmin;
    double MaxX = GlobalXmax;
    bool find_root = true;
    double root = -999;
#ifdef __debug_output__
    std::cout   << "RootFinder::GetRoots >> GlobalXmin: " << GlobalXmin
    << " GlobalXmax: " << GlobalXmax
    << std::endl;
#endif //debug output
    
    // Check if TGraph and TF1 have been set
    if(graph_eval == 0)
    {
        std::cerr << "TGraph not set!" << std::endl;
        find_root = false;
    }
    
    if(dchisq_val == 0)
    {
        std::cerr << "Dchisq not set!" << std::endl;
        find_root = false;
    }
    
    while(MinX < GlobalXmax){
        double x1 = MinX;
        double x2 = MinX + step_size;
        root_function = new TF1("f1",this,&RootFinder::find_root_function,x1,x2,0,"RootFinder","find_root_function");
        root = root_function->GetMinimumX();
#ifdef __verbose_debug_output__
        std::cout   << "RootFinder::GetRoots >> Xmin: " << MinX
        << " Xmax: " << MaxX << " Root: " << root
        << std::endl;
#endif //debug output
        
        if(TMath::Abs(root_function->Eval(root)) < 0.00001){
            roots.push_back(root);
        }

        MinX = MinX + step_size;
    }
    
    this->MakeRootTypes();
    
#ifdef __debug_output__

    std::cout << "RootFinder >> Roots Found: " << std::endl;
    
    for(int i=0;i<(int)roots.size();i++)
    {
        std::cout << roots[i] << "\t gradient: " << root_types[i] << std::endl;
    }
#endif //debug output


    return roots;
}
//___________________________________________________________________________
void RootFinder::MakeRootTypes(void)
{
    for(int i = 0; i< (int)roots.size();i++)
    {
        if((graph_eval->Eval(roots[i] + 0.000001,0,"S") - graph_eval->Eval(roots[i],0,"S")) > 0)
        {
            root_types.push_back(1);
        }
        else{
            root_types.push_back(-1);
        }
    }

}
//___________________________________________________________________________
double RootFinder::find_root_function(double *x, double*par)
{
    double a = graph_eval->Eval(*x,0,"S");
    double b = dchisq_val->EvalPar(x,par);
    
    return TMath::Abs(a - b);
}
//___________________________________________________________________________
void RootFinder::CleanUp(void)
{
    if ( graph_eval     )       delete graph_eval;
    if ( dchisq_val     )       delete dchisq_val;
    if ( root_function  )       delete root_function;
}
//___________________________________________________________________________
void RootFinder::Init(void)
{
    graph_eval      = 0;
    dchisq_val      = 0;
    root_function   = 0;
    GlobalXmin      = -999;
    GlobalXmax      = -999;
    step_size       = 1000;
}
//___________________________________________________________________________
double RootFinder::GetError()
{
    double start_min    = GlobalXmin;
    
    double error        = 0;
//    std::cout << "Nroots: " << (int)roots.size() << std::endl;
    for(int i=0; i<(int)roots.size(); i++)
    {
        if(root_types[i] < 0)
        {
        start_min = roots[i];
        }
        else if(root_types[i] > 0)
        {
            error += roots[i] - start_min;
        }
    }
    //Check gradient of last minium
    if(roots.size()>0){
    if(root_types[roots.size()-1] < 0) {
        error += GlobalXmax - roots[roots.size()-1];
#ifdef __show_rootfinder_debug__
        std::cout << "RootFinder >> Last root has neg gradient adding up to max " << GlobalXmax - roots[roots.size()-1] << std::endl;
#endif //rootfinder debug
    }

    }
    else error = GlobalXmax - GlobalXmin;

    
    return error;
}
//___________________________________________________________________________
void RootFinder::SetRange(void)
{
    int points = graph_eval->GetN();
    double min_x_val, max_x_val, dummy_y;
    graph_eval->GetPoint(0,min_x_val,dummy_y);
    graph_eval->GetPoint(points -1,max_x_val,dummy_y);
    
    GlobalXmin = min_x_val;
    GlobalXmax = max_x_val;
    
    step_size = (GlobalXmax - GlobalXmin) / 100;
    
}

void RootFinder::PrintRoots(void)
{
  cout << "RootFinder roots: " << endl;

  for (unsigned i=0; i < roots.size(); i++) {
    cout << "\t" << roots[i] << endl;
  }
}
