//____________________________________________________________________________
/*!
 
 \class     valor::utils
 
 \brief     Utilities for formatting and style
 
 \author    Raj Shah <raj.shah \at stfc.ac.uk>
 University of Oxford & STFC Rutherford Appleton Laboratory
 
 \created   Oct 10, 2015
 
 \copyright (c) 2010-2015, The VALOR group.
 
 */
//____________________________________________________________________________

#ifndef _valor_rootfinder_
#define _valor_rootfinder_

#include <iostream>
#include <string>
#include <vector>
#include <TF1.h>
#include <TH1.h>
#include <TGraph.h>
#include <TSystem.h>


using std::vector;

namespace valor {
    class  RootFinder {
        
    public:
        RootFinder();
        virtual ~RootFinder();
        double GetError                 ();
        void SetGraph                   (TGraph* g)             {graph_eval = new TGraph(*g); graph_eval->Sort();this->SetRange();}
        void SetGraph                   (TH1* h)                {graph_eval = new TGraph(h);  graph_eval->Sort();this->SetRange();}
        void SetRange                   (double x1, double x2)  {GlobalXmax = x2; GlobalXmin = x1; }
        void SetRange                   (void);
        void SetDchiSq                  (double chisq)          {dchisq_val = new TF1("f1",Form("%f",chisq),-20,20);}
        void SetDchiSqFunc              (TF1* f)                {dchisq_val = f;}
        vector<double> GetRoots         (void);
        vector<double> GetRootTypes     (void)                  {return root_types;}
        void PrintRoots                 (void);
    private:
        
        double find_root_function(double *x, double*par);
        TGraph* graph_eval;
        TF1* dchisq_val;
        TF1* root_function;
        void CleanUp (void);
        void Init (void);
        vector<double> roots;
        vector<double> root_types;
        void MakeRootTypes     (void);
        
        double GlobalXmin;
        double GlobalXmax;
        double step_size;
        ClassDef(RootFinder,1);
    }; // root finder
} // valor

#endif

