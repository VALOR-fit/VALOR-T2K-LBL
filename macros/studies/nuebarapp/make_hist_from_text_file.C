/*!
\macro   make_hist_from_text_file.C

\brief   Read in data from plaintext files into user-specified histograms.

\details The input data files should be arranged into whitespace-seperated columns (2 or 3 columns required),
         which will be stored as x y z into the histograms that the user specifies.
         WARNING: Make sure that there are no empty lines in the input file, as these will be stored in the histogram
         as (0,0,0) and will potentially cause problems with the first histogram bin.

\author  Francis Bench <fbench@hep.ph.liv.ac.uk>
         University of Liverpool

\created Nov 2017

\cpright Copyright (c) 2010-2018, The VALOR Fitting Group, All Rights Reserved
*/

#include "Riostream.h"
#include <stdlib.h>

void Read(TH1D *data, std::string data_file_name, std::string out_file_name, bool weighted, bool average_weighted, bool interpolate) ;
void Read(TH2D *data, std::string data_file_name, std::string out_file_name, bool weighted, bool average_weighted);

void make_hist_from_text_file(){
    TH1D *hist_Nobs_beta0 = new TH1D("Nobs_beta0", "Nobs_beta0", 66, -1, 65);
    hist_Nobs_beta0->SetXTitle("N_{obs}");
    hist_Nobs_beta0->SetYTitle("Entries");
    Read(hist_Nobs_beta0, "/hepstore/fbench/software/valor/t2xk/branches/2018_v0.0.1/P-theta_nuebarapp_2017-12-06/rate-only_beta=0.txt",
                          "/hepstore/fbench/software/valor/t2xk/branches/2018_v0.0.1/P-theta_nuebarapp_2017-12-06/rate-only_beta=0.root"
                          , true, true, true);
    TCanvas *canvas_Nobs_beta0 = new TCanvas("Nobs_beta0", "Nobs_beta0", 800, 600);
    hist_Nobs_beta0->Draw();
}


// Read in data and store in a 2D histogram, then write the histogram to a file. Assumes columns of whitespace-delimited data.
// weighted controls whether the points added to the histogram should increment the bin by 1 (weighted = false) or by a weight given the next column of data (weighted = true)
// If average_weighted = true, then the content of each bin is averaged by the number of points added to the bin
void Read(TH2D *data, std::string data_file_name, std::string out_file_name, bool weighted, bool average_weighted){
    ifstream in;
    in.open(data_file_name.c_str(), std::ifstream::in);

    double x,y,z;
    int nlines = 0;
    std::string line;
    TH2D *count = data->Clone("Counts");
    TFile *out_file = new TFile(out_file_name.c_str(),"RECREATE");

    cout << "For file " << data_file_name << ": \n";

    while (1) {
        getline(in, line);
        istringstream iss(line);
        std::string str;
        iss >> str;
        x = atof(str.c_str());
        iss >> std::str;
        y = atof(str.c_str());;
        if(weighted) {
            //in >> x >> y >> z;
            iss >> str;
            z = atof(str.c_str());
            if (!in.good()) break;
            if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n", x, y, z);
            data->Fill(x, y, z);
            if(average_weighted) count->Fill(x, y);
        }
        else{
            //in >> x >> y;
            if (!in.good()) break;
            if (nlines < 5) printf("x=%8f, y=%8f\n", x, y);
            data->Fill(x, y);
        }
        nlines++;
    }

    if(average_weighted && weighted){
        int n_bins = data->GetSize();
        for(int bin = 0; bin < n_bins; ++bin){
            if(data->IsBinUnderflow(bin) || data->IsBinOverflow(bin)) continue;
            double content = data->GetBinContent(bin);
            double counts = count->GetBinContent(bin);
            if(counts != 0) data->SetBinContent(bin, content / counts);
        }
    }

    printf(" found %d points\n \n ",nlines);
    in.close();
    data->Write(data->GetTitle());
    out_file->Close();
}


// Read in data and store in a 1D histogram, then write the histogram to a file. Assumes columns of whitespace-delimited data.
// weighted controls whether the points added to the histogram should increment the bin by 1 (weighted = false) or by a weight given the next column of data (weighted = true)
// If average_weighted = true, then the content of each bin is averaged by the number of points added to the bin
// If any bins have 0 content and the bins either side of it has content, assume that this is an error from the plot
// digitization and average the contents of the bins either side to interpolate a value for this bin (only for interpolate = true)
void Read(TH1D *data, std::string data_file_name, std::string out_file_name, bool weighted, bool average_weighted, bool interpolate){
    ifstream in;
    in.open(data_file_name.c_str(), std::ifstream::in);

    double x,y;
    int nlines = 0;
    std::string line;
    TH1D *count = data->Clone("Counts");
    TFile *out_file = new TFile(out_file_name.c_str(),"RECREATE");

    cout << "For file " << data_file_name << ": \n";

    while (1) {
        getline(in, line);
        istringstream iss(line);
        std::string str;
        iss >> str;
        x = atof(str.c_str());
        if(weighted) {
            iss >> str;
            y = atof(str.c_str());
            if (!in.good()) break;
            if (nlines < 5) printf("x=%8f, y=%8f\n", x, y);
            data->Fill(x, y);
            if(average_weighted) count->Fill(x);
            //printf("x=%8f, y=%8f\n", x, y);
        }
        else{
            if (!in.good()) break;
            if (nlines < 5) printf("x=%8f\n", x);
            data->Fill(x);
        }
        nlines++;
    }

    if(average_weighted && weighted){
        int n_bins = data->GetSize();
        for(int bin = 0; bin < n_bins; ++bin){
            if(data->IsBinUnderflow(bin) || data->IsBinOverflow(bin)) continue;
            double content = data->GetBinContent(bin);
            double counts = count->GetBinContent(bin);
            if(counts != 0) data->SetBinContent(bin, content / counts);
        }
    }

    // If any bins have 0 content (except for the first and last bins), assume that this is an error from the plot digitization and linearly interpolate that bin
    if(interpolate) {
        int n_bins = data->GetSize();
        for (int bin = 1; bin < n_bins - 1; ++bin) {
            if (data->IsBinUnderflow(bin) || data->IsBinOverflow(bin)) continue;
            double content = data->GetBinContent(bin);
            double counts_left  = count->GetBinContent(bin - 1);
            double counts_right = count->GetBinContent(bin + 1);
            if (content == 0 && counts_left != 0 && counts_right != 0) {
              double bin_centre = data->GetBinCenter(bin);
              //double new_value = data->Interpolate(bin_centre);
              double new_value = (data->GetBinContent(bin - 1) + data->GetBinContent(bin + 1)) / 2.0;
              data->SetBinContent(bin, new_value);
              cout << "Found zero content at x = " << bin_centre << "; Content set to a new value of " << new_value << endl;
            }
        }
    }

    printf(" found %d points\n \n ",nlines);
    in.close();
    data->Write(data->GetTitle());
    out_file->Close();
}