#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};



struct peak_fit {
        float   	parameters[3];
        float       err_parameters[3];
        float	    resolution;
        float   	err_resolution;
};



const int n_Am = 1;
const int n_Na = 2;
const int n_Eu = 12;


float peak_1_1_Am[n_Am];
float peak_0_2_Am[n_Am];
float err_peak_1_1_Am[n_Am];
float err_peak_0_2_Am[n_Am];

float peak_1_1_Na[n_Na];
float peak_0_2_Na[n_Na];
float err_peak_1_1_Na[n_Na];
float err_peak_0_2_Na[n_Na];

float peak_0_1_Eu[n_Eu];
float peak_1_1_Eu[n_Eu];
float err_peak_0_1_Eu[n_Eu];
float err_peak_1_1_Eu[n_Eu];





// ********************************************************************************************
// ********************************************************************************************
void getPeakAm(int draw=1) {
    // variables
    const int n_peaks = 1;
    float bins_NaI[n_peaks] = {
        50
    };
    float bins_HPGe[n_peaks] = {
        70
    };
    float range_NaI[n_peaks][2] = {
        {   130,    280   },
    };
    float range_HPGe[n_peaks][2] = {
        {   150,    220   }
    };
    float fit_range_NaI[n_peaks][2] = {
        {   160,    240   },
    };
    float fit_range_HPGe[n_peaks][2] = {
        {   198,    212   },
    };
    float R_0_1[n_peaks];
    float R_1_1[n_peaks];
    float R_0_2[n_peaks];
    float R_1_2[n_peaks];
    float err_R_0_1[n_peaks];
    float err_R_1_1[n_peaks];
    float err_R_0_2[n_peaks];
    float err_R_1_2[n_peaks];
	slimport_data_t indata_0_1[n_peaks];
    slimport_data_t indata_1_1[n_peaks];
    slimport_data_t indata_0_2[n_peaks];
    slimport_data_t indata_1_2[n_peaks];

	TFile *infile_1 = new TFile("DATA/Day1/Am241_calibration_1.root");
	TTree *intree_1 = (TTree*)infile_1->Get("acq_tree_0");
    TFile *infile_2 = new TFile("DATA/Day1/Am241_calibration_2.root");
	TTree *intree_2 = (TTree*)infile_2->Get("acq_tree_0");

    TBranch** inbranch_0_1 = new TBranch*[n_peaks];
    TBranch** inbranch_1_1 = new TBranch*[n_peaks];
    inbranch_0_1[0] = new TBranch;
    inbranch_1_1[0] = new TBranch;
	inbranch_0_1[0] = intree_1->GetBranch("acq_ch0");
    inbranch_1_1[0] = intree_1->GetBranch("acq_ch1");
	inbranch_0_1[0]->SetAddress(&indata_0_1[0].timetag);
    inbranch_1_1[0]->SetAddress(&indata_1_1[0].timetag);

    TBranch** inbranch_0_2 = new TBranch*[n_peaks];
    TBranch** inbranch_1_2 = new TBranch*[n_peaks];
    inbranch_0_2[0] = new TBranch;
    inbranch_1_2[0] = new TBranch;
	inbranch_0_2[0] = intree_2->GetBranch("acq_ch0");
    inbranch_1_2[0] = intree_2->GetBranch("acq_ch1");
	inbranch_0_2[0]->SetAddress(&indata_0_2[0].timetag);
    inbranch_1_2[0]->SetAddress(&indata_1_2[0].timetag);

    TH1F** h_0_1 = new TH1F*[n_peaks];
    TH1F** h_1_1 = new TH1F*[n_peaks];
	h_0_1[0] = new TH1F("NaI detector P1","NaI detector spectrum peak at 59.5 keV", bins_NaI[0], range_NaI[0][0], range_NaI[0][1]);
    h_1_1[0] = new TH1F("HPGe detector P1","HPGe detector spectrum peak at 59.5 keV", bins_HPGe[0], range_HPGe[0][0], range_HPGe[0][1]);

    TH1F** h_0_2 = new TH1F*[n_peaks];
    TH1F** h_1_2 = new TH1F*[n_peaks];
	h_0_2[0] = new TH1F("NaI detector P2","NaI detector spectrum peak at 59.5 keV", bins_NaI[0], range_NaI[0][0], range_NaI[0][1]);
    h_1_2[0] = new TH1F("HPGe detector P2","HPGe detector spectrum peak at 59.5 keV", bins_HPGe[0], range_HPGe[0][0], range_HPGe[0][1]);

	// histogram filling
    for (int i=0; i<n_peaks; ++i) {
        for (int j=0; j<inbranch_0_1[i]->GetEntries(); j++) {
            inbranch_0_1[i]->GetEntry(j);
            h_0_1[i]->Fill(indata_0_1[i].qlong);
        }
        for (int j=0; j<inbranch_1_1[i]->GetEntries(); j++) {
            inbranch_1_1[i]->GetEntry(j);
            h_1_1[i]->Fill(indata_1_1[i].qlong);
        }
        for (int j=0; j<inbranch_0_2[i]->GetEntries(); j++) {
            inbranch_0_2[i]->GetEntry(j);
            h_0_2[i]->Fill(indata_0_2[i].qlong);
        }
        for (int j=0; j<inbranch_1_2[i]->GetEntries(); j++) {
            inbranch_1_2[i]->GetEntry(j);
            h_1_2[i]->Fill(indata_1_2[i].qlong);
        }
    }

    TF1** f_0_1 = new TF1*[n_peaks];
    TF1** f_1_1 = new TF1*[n_peaks];
    TF1** f_0_2 = new TF1*[n_peaks];
    TF1** f_1_2 = new TF1*[n_peaks];
    
    for (int i=0; i<n_peaks; ++i) {
        f_0_1[i] = new TF1("fit_0_1","gaus",fit_range_NaI[i][0],fit_range_NaI[i][1]);
        f_1_1[i] = new TF1("fit_1_1","gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        f_0_2[i] = new TF1("fit_0_2","gaus",fit_range_NaI[i][0],fit_range_NaI[i][1]);
        f_1_2[i] = new TF1("fit_1_2","gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
    }

    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Position 1", 800, 400*n_peaks);
	    c1->Divide(2,n_peaks);
        for (int i=0; i<n_peaks; ++i) {
            c1->cd(2*i+1);
            h_0_1[i]->Draw();
            h_0_1[i]->Fit(f_0_1[i],"","",fit_range_NaI[i][0],fit_range_NaI[i][1]);
            c1->cd(2*i+2);
            h_1_1[i]->Draw();
            h_1_1[i]->Fit(f_1_1[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
        
        TCanvas* c2 = new TCanvas("c2", "Position 2", 800, 400*n_peaks);
	    c2->Divide(2,n_peaks);
        for (int i=0; i<n_peaks; ++i) {
            c2->cd(2*i+1);
            h_0_2[i]->Draw();
            h_0_2[i]->Fit(f_0_2[i],"","",fit_range_NaI[i][0],fit_range_NaI[i][1]);
            c2->cd(2*i+2);
            h_1_2[i]->Draw();
            h_1_2[i]->Fit(f_1_2[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
    }

    for (int i=0; i<n_peaks; ++i) {
        R_0_1[i] = 2.355 * f_0_1[i]->GetParameter(2) / f_0_1[i]->GetParameter(1);
        R_1_1[i] = 2.355 * f_1_1[i]->GetParameter(2) / f_1_1[i]->GetParameter(1);
        R_0_2[i] = 2.355 * f_0_2[i]->GetParameter(2) / f_0_2[i]->GetParameter(1);
        R_1_2[i] = 2.355 * f_1_2[i]->GetParameter(2) / f_1_2[i]->GetParameter(1);
        err_R_0_1[i] = R_0_1[i] * sqrt( pow(f_0_1[i]->GetParError(2)/f_0_1[i]->GetParameter(2),2.0) + pow(f_0_1[i]->GetParError(1)/f_0_1[i]->GetParameter(1),2.0) );
        err_R_1_1[i] = R_1_1[i] * sqrt( pow(f_1_1[i]->GetParError(2)/f_1_1[i]->GetParameter(2),2.0) + pow(f_1_1[i]->GetParError(1)/f_1_1[i]->GetParameter(1),2.0) );
        err_R_0_2[i] = R_0_2[i] * sqrt( pow(f_0_2[i]->GetParError(2)/f_0_2[i]->GetParameter(2),2.0) + pow(f_0_2[i]->GetParError(1)/f_0_2[i]->GetParameter(1),2.0) );
        err_R_1_2[i] = R_1_2[i] * sqrt( pow(f_1_2[i]->GetParError(2)/f_1_2[i]->GetParameter(2),2.0) + pow(f_1_2[i]->GetParError(1)/f_1_2[i]->GetParameter(1),2.0) );
    }

    ofstream ofile("Analysis/Day1/Day1_Am241.txt");
    ofile << "Parameters of gaussian fits for D1P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_0_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_0_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D2P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D1P2: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_0_2[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_0_2[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D2P2: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_2[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_2[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile.close();

}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getPeakNa(int draw=1) {

    // variables
    const int n_peaks = 2;
    float bins_NaI[n_peaks] = {
        50,     50
    };
    float bins_HPGe[n_peaks] = {
        70,     70
    };
    float range_NaI[n_peaks][2] = {
        {   1300,    1800   },
        {   3200,    4200   }
    };
    float range_HPGe[n_peaks][2] = {
        {   1685,    1755   },
        {   4250,    4320   }
    };
    float fit_range_NaI[n_peaks][2] = {
        {   1400,    1650   },
        {   3500,    3900   }
    };
    float fit_range_HPGe[n_peaks][2] = {
        {   1700,    1730   },
        {   4275,    4295   }
    };
    float R_0_1[n_peaks];
    float R_1_1[n_peaks];
    float R_0_2[n_peaks];
    float R_1_2[n_peaks];
    float err_R_0_1[n_peaks];
    float err_R_1_1[n_peaks];
    float err_R_0_2[n_peaks];
    float err_R_1_2[n_peaks];
	slimport_data_t indata_0_1;
    slimport_data_t indata_1_1;
    slimport_data_t indata_0_2;
    slimport_data_t indata_1_2;

	TFile *infile_1 = new TFile("DATA/Day1/Na22_calibration_1.root");
	TTree *intree_1 = (TTree*)infile_1->Get("acq_tree_0");
    TFile *infile_2 = new TFile("DATA/Day1/Na22_calibration_2.root");
	TTree *intree_2 = (TTree*)infile_2->Get("acq_tree_0");

    TBranch* inbranch_0_1 = intree_1->GetBranch("acq_ch0");
    TBranch* inbranch_1_1 = intree_1->GetBranch("acq_ch1");
    inbranch_0_1->SetAddress(&indata_0_1.timetag);
    inbranch_1_1->SetAddress(&indata_1_1.timetag);

    TBranch* inbranch_0_2 = intree_2->GetBranch("acq_ch0");
    TBranch* inbranch_1_2 = intree_2->GetBranch("acq_ch1");
    inbranch_0_2->SetAddress(&indata_0_2.timetag);
    inbranch_1_2->SetAddress(&indata_1_2.timetag);


    TH1F** h_0_1 = new TH1F*[n_peaks];
    TH1F** h_1_1 = new TH1F*[n_peaks];
	h_0_1[0] = new TH1F("NaI detector P1 Pk1","NaI detector spectrum peak at 511 keV", bins_NaI[0], range_NaI[0][0], range_NaI[0][1]);
    h_0_1[1] = new TH1F("NaI detector P1 Pk2","NaI detector spectrum peak at 1275 keV", bins_NaI[1], range_NaI[1][0], range_NaI[1][1]);
    h_1_1[0] = new TH1F("HPGe detector P1 Pk1","HPGe detector spectrum peak at 511 keV", bins_HPGe[0], range_HPGe[0][0], range_HPGe[0][1]);
    h_1_1[1] = new TH1F("HPGe detector P1 Pk2","HPGe detector spectrum peak at 1275 keV", bins_HPGe[1], range_HPGe[1][0], range_HPGe[1][1]);

    TH1F** h_0_2 = new TH1F*[n_peaks];
    TH1F** h_1_2 = new TH1F*[n_peaks];
	h_0_2[0] = new TH1F("NaI detector P2 Pk1","NaI detector spectrum peak at 511 keV", bins_NaI[0], range_NaI[0][0], range_NaI[0][1]);
    h_0_2[1] = new TH1F("NaI detector P2 Pk2","NaI detector spectrum peak at 1275 keV", bins_NaI[1], range_NaI[1][0], range_NaI[1][1]);
    h_1_2[0] = new TH1F("HPGe detector P2 Pk1","HPGe detector spectrum peak at 511 keV", bins_HPGe[0], range_HPGe[0][0], range_HPGe[0][1]);
    h_1_2[1] = new TH1F("HPGe detector P2 Pk2","HPGe detector spectrum peak at 1275 keV", bins_HPGe[1], range_HPGe[1][0], range_HPGe[1][1]);    


	// histogram filling
    for (int i=0; i<n_peaks; ++i) {
        for (int j=0; j<inbranch_0_1->GetEntries(); j++) {
            inbranch_0_1->GetEntry(j);
            h_0_1[i]->Fill(indata_0_1.qlong);
        }
        for (int j=0; j<inbranch_1_1->GetEntries(); j++) {
            inbranch_1_1->GetEntry(j);
            h_1_1[i]->Fill(indata_1_1.qlong);
        }
        for (int j=0; j<inbranch_0_2->GetEntries(); j++) {
            inbranch_0_2->GetEntry(j);
            h_0_2[i]->Fill(indata_0_2.qlong);
        }
        for (int j=0; j<inbranch_1_2->GetEntries(); j++) {
            inbranch_1_2->GetEntry(j);
            h_1_2[i]->Fill(indata_1_2.qlong);
        }
    }

    TF1** f_0_1 = new TF1*[n_peaks];
    TF1** f_1_1 = new TF1*[n_peaks];
    TF1** f_0_2 = new TF1*[n_peaks];
    TF1** f_1_2 = new TF1*[n_peaks];
    
    for (int i=0; i<n_peaks; ++i) {
        f_0_1[i] = new TF1(Form("fit_0_1_%d",i+1),"gaus",fit_range_NaI[i][0],fit_range_NaI[i][1]);
        f_1_1[i] = new TF1(Form("fit_1_1_%d",i+1),"gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        f_0_2[i] = new TF1(Form("fit_0_2_%d",i+1),"gaus",fit_range_NaI[i][0],fit_range_NaI[i][1]);
        f_1_2[i] = new TF1(Form("fit_1_2_%d",i+1),"gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
    }


    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Position 1", 800, 400*n_peaks);
	    c1->Divide(2,n_peaks);
        for (int i=0; i<n_peaks; ++i) {
            c1->cd(2*i+1);
            h_0_1[i]->Draw();
            h_0_1[i]->Fit(f_0_1[i],"","",fit_range_NaI[i][0],fit_range_NaI[i][1]);
            c1->cd(2*i+2);
            h_1_1[i]->Draw();
            h_1_1[i]->Fit(f_1_1[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
        
        TCanvas* c2 = new TCanvas("c2", "Position 2", 800, 400*n_peaks);
	    c2->Divide(2,n_peaks);
        for (int i=0; i<n_peaks; ++i) {
            c2->cd(2*i+1);
            h_0_2[i]->Draw();
            h_0_2[i]->Fit(f_0_2[i],"","",fit_range_NaI[i][0],fit_range_NaI[i][1]);
            c2->cd(2*i+2);
            h_1_2[i]->Draw();
            h_1_2[i]->Fit(f_1_2[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
    }

    for (int i=0; i<n_peaks; ++i) {
        R_0_1[i] = 2.355 * f_0_1[i]->GetParameter(2) / f_0_1[i]->GetParameter(1);
        R_1_1[i] = 2.355 * f_1_1[i]->GetParameter(2) / f_1_1[i]->GetParameter(1);
        R_0_2[i] = 2.355 * f_0_2[i]->GetParameter(2) / f_0_2[i]->GetParameter(1);
        R_1_2[i] = 2.355 * f_1_2[i]->GetParameter(2) / f_1_2[i]->GetParameter(1);
        err_R_0_1[i] = R_0_1[i] * sqrt( pow(f_0_1[i]->GetParError(2)/f_0_1[i]->GetParameter(2),2.0) + pow(f_0_1[i]->GetParError(1)/f_0_1[i]->GetParameter(1),2.0) );
        err_R_1_1[i] = R_1_1[i] * sqrt( pow(f_1_1[i]->GetParError(2)/f_1_1[i]->GetParameter(2),2.0) + pow(f_1_1[i]->GetParError(1)/f_1_1[i]->GetParameter(1),2.0) );
        err_R_0_2[i] = R_0_2[i] * sqrt( pow(f_0_2[i]->GetParError(2)/f_0_2[i]->GetParameter(2),2.0) + pow(f_0_2[i]->GetParError(1)/f_0_2[i]->GetParameter(1),2.0) );
        err_R_1_2[i] = R_1_2[i] * sqrt( pow(f_1_2[i]->GetParError(2)/f_1_2[i]->GetParameter(2),2.0) + pow(f_1_2[i]->GetParError(1)/f_1_2[i]->GetParameter(1),2.0) );
    }

    ofstream ofile("Analysis/Day1/Day1_Na22.txt");
    ofile << "Parameters of gaussian fits for D1P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_0_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_0_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D2P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D1P2: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_2[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_0_2[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_0_2[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D2P2: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_2[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_2[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_2[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile.close();

}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getPeakEu(int draw=1) {
    // variables
    const int n_peaks = 7;
    float bins_NaI[n_peaks] = {
        50,     40,    40,   40,    50,    50,    50
    };
    float bins_HPGe[n_peaks] = {
        60,     80,     80,     60,     80,     100,    100
    };
    float range_NaI[n_peaks][2] = {
        {    310,    460    },
        {    660,    820    },
        {    900,   1170    },
        {   2100,   2490    },
        {   2600,   3000    },
        {   3000,   3400    },
        {   3750,   4350    }
    };
    float range_HPGe[n_peaks][2] = {
        {    380,    440    },
        {    780,    860    },
        {   1120,   1200    },
        {   2590,   2650    },
        {   3200,   3280    },
        {   3700,   3780    },
        {   4680,   4780    }
    };
    float fit_range_NaI[n_peaks][2] = {
        {    340,    430    },
        {    690,    790    },
        {    950,   1120    },
        {   2190,   2400    },
        {   2700,   2900    },
        {   3100,   3300    },
        {   3900,   4200    }
    };
    float fit_range_HPGe[n_peaks][2] = {
        {    408,    420    },
        {    820,    830    },
        {   1150,   1170    },
        {   2615,   2630    },
        {   3235,   3250    },
        {   3735,   3750    },
        {   4715,   4745    }
    };
    float R_0_1[n_peaks];
    float R_1_1[n_peaks];
    float err_R_0_1[n_peaks];
    float err_R_1_1[n_peaks];
	slimport_data_t indata_0_1;
    slimport_data_t indata_1_1;


	TFile *infile_1 = new TFile("DATA/Day1/Eu152_calibration_1.root");
	TTree *intree_1 = (TTree*)infile_1->Get("acq_tree_0");

    TBranch* inbranch_0_1 = intree_1->GetBranch("acq_ch0");
    TBranch* inbranch_1_1 = intree_1->GetBranch("acq_ch1");
    inbranch_0_1->SetAddress(&indata_0_1.timetag);
    inbranch_1_1->SetAddress(&indata_1_1.timetag);


    TH1F** h_0_1 = new TH1F*[n_peaks];
    TH1F** h_1_1 = new TH1F*[n_peaks];
    for (int i=0; i<n_peaks; ++i) {
        h_0_1[i] = new TH1F(Form("NaI detector P1 %d-th peak",i+1),Form("NaI detector P1 %d-th peak",i+1), bins_NaI[i], range_NaI[i][0], range_NaI[i][1]);
        h_1_1[i] = new TH1F(Form("HPGe detector P1 %d-th peak",i+1),Form("HPGe detector P1 %d-th peak",i+1), bins_HPGe[i], range_HPGe[i][0], range_HPGe[i][1]);
    }


	// histogram filling
    for (int i=0; i<n_peaks; ++i) {
        for (int j=0; j<inbranch_0_1->GetEntries(); j++) {
            inbranch_0_1->GetEntry(j);
            h_0_1[i]->Fill(indata_0_1.qlong);
        }
        for (int j=0; j<inbranch_1_1->GetEntries(); j++) {
            inbranch_1_1->GetEntry(j);
            h_1_1[i]->Fill(indata_1_1.qlong);
        }
    }

    TF1** f_0_1 = new TF1*[n_peaks];
    TF1** f_1_1 = new TF1*[n_peaks];
    
    for (int i=0; i<n_peaks; ++i) {
        f_0_1[i] = new TF1(Form("fit_0_1_%d",i+1),"gaus",fit_range_NaI[i][0],fit_range_NaI[i][1]);
        f_1_1[i] = new TF1(Form("fit_1_1_%d",i+1),"gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
    }


    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Position 1", 800, 400*n_peaks);
	    c1->Divide(2,n_peaks);
        for (int i=0; i<n_peaks; ++i) {
            c1->cd(2*i+1);
            h_0_1[i]->Draw();
            h_0_1[i]->Fit(f_0_1[i],"","",fit_range_NaI[i][0],fit_range_NaI[i][1]);
            c1->cd(2*i+2);
            h_1_1[i]->Draw();
            h_1_1[i]->Fit(f_1_1[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
    }


    for (int i=0; i<n_peaks; ++i) {
        peak_0_1_Eu[i] = f_0_1[i]->GetParameter(1);
        peak_1_1_Eu[i] = f_1_1[i]->GetParameter(1);
        err_peak_0_1_Eu[i] = f_0_1[i]->GetParError(1);
        err_peak_1_1_Eu[i] = f_1_1[i]->GetParError(1);
    }


    for (int i=0; i<n_peaks; ++i) {
        R_0_1[i] = 2.355 * f_0_1[i]->GetParameter(2) / f_0_1[i]->GetParameter(1);
        R_1_1[i] = 2.355 * f_1_1[i]->GetParameter(2) / f_1_1[i]->GetParameter(1);
        err_R_0_1[i] = R_0_1[i] * sqrt( pow(f_0_1[i]->GetParError(2)/f_0_1[i]->GetParameter(2),2.0) + pow(f_0_1[i]->GetParError(1)/f_0_1[i]->GetParameter(1),2.0) );
        err_R_1_1[i] = R_1_1[i] * sqrt( pow(f_1_1[i]->GetParError(2)/f_1_1[i]->GetParameter(2),2.0) + pow(f_1_1[i]->GetParError(1)/f_1_1[i]->GetParameter(1),2.0) );
    }


    ofstream ofile("Analysis/Day1/Day1_Eu152.txt");
    ofile << "Parameters of gaussian fits for D1P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_0_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_0_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_0_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Parameters of gaussian fits for D2P1: [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                        << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_1[i]             << '\t' << '$' << '\t'
                                        << "\\\\" << endl;
    }
    ofile.close();

}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getPeakEuHPG(int draw=1) {
    // variables
    const int n_peaks = 12;
    float bins_HPGe[n_peaks] = {
        80,   // 1
        100,  // 2
        100,  // 3
        100,  // 4
        100,  // 5
        80,   // 6
        80,   // 7
        100,  // 8
        100,  // 9
        120, // 10
        100,  // 11
        120  // 12
    };
    float range_HPGe[n_peaks][2] = {
        {    380,    460    }, // 1
        {    780,    880    }, // 2
        {   1120,   1220    }, // 3
        {   1200,   1300    }, // 4
        {   1340,   1440    }, // 5
        {   1460,   1540    }, // 6
        {   2590,   2670    }, // 7
        {   2880,   2980    }, // 8
        {   3200,   3300    }, // 9
        {   3600,   3720    }, // 10
        {   3700,   3800    }, // 11
        {   4680,   4800    }  // 12
    };
    float fit_range_HPGe[n_peaks][2] = {
        {    408,    420    }, // 1
        {    820,    830    }, // 2
        {   1155,   1170    }, // 3
        {   1232,   1246    }, // 4
        {   1380,   1390    }, // 5        
        {   1490,   1502    }, // 6
        {   2615,   2630    }, // 7
        {   2912,   2924    }, // 8
        {   3237,   3250    }, // 9
        {   3645,   3660    }, // 10
        {   3735,   3750    }, // 11
        {   4728,   4748    }  // 12
    };
    float E_peak_keV[n_peaks] = {
        121.8,  // 1
        244.7,  // 2
        344.3,  // 3
        367.8,  // 4
        411.1,  // 5 
        444.0,  // 6
        778.9,  // 7
        867.4,  // 8
        964.0,  // 9
        1085.8, // 10
        1112.1, // 11
        1408.0  // 12
    };
    float weights[n_peaks] = {
        0.2841,  // 1
        0.0755,  // 2
        0.2659,  // 3
        0.0086,  // 4
        0.0224,  // 5 
        0.0312,  // 6
        0.1297,  // 7
        0.0424,  // 8
        0.1463,  // 9
        0.1015,  // 10
        0.1360,  // 11
        0.2085   // 12
    };
    float rel_weights[n_peaks];
    int peak_1408_index = n_peaks - 1;
    for (int i=0; i<n_peaks; ++i) rel_weights[i] = 100.0 * (weights[i]/weights[peak_1408_index]);
    float err_E_peak_keV[n_peaks] = {
        0.00001,  // 1
        0.00001,  // 2
        0.00001,  // 3
        0.00001,  // 4
        0.00001,  // 5 
        0.00001,  // 6
        0.00001,  // 7
        0.00001,  // 8
        0.00001,  // 9
        0.00001,  // 10
        0.00001,  // 11
        0.00001   // 12
    };
    float R_1_1[n_peaks];
    float err_R_1_1[n_peaks];
    slimport_data_t indata_1_1;
    slimport_data_t indata_1_b;


	TFile* infile_1 = new TFile("DATA/Day1/Eu152_calibration_1.root");
    TFile* infile_b = new TFile("DATA/Day1/bkg.root");
	TTree* intree_1 = (TTree*)infile_1->Get("acq_tree_0");
    TTree* intree_b = (TTree*)infile_b->Get("acq_tree_0");

    TBranch* inbranch_1_1 = intree_1->GetBranch("acq_ch1");
    TBranch* inbranch_1_b = intree_b->GetBranch("acq_ch1");
    inbranch_1_1->SetAddress(&indata_1_1.timetag);
    inbranch_1_b->SetAddress(&indata_1_b.timetag);


    TH1F** h_1_1 = new TH1F*[n_peaks];
    TH1F** h_1_b = new TH1F*[n_peaks];
    for (int i=0; i<n_peaks; ++i) {
        h_1_1[i] = new TH1F(Form("HPGe detector P1 %d-th peak",i+1),Form("HPGe detector P1 %d-th peak",i+1), bins_HPGe[i], range_HPGe[i][0], range_HPGe[i][1]);
        h_1_b[i] = new TH1F(Form("HPGe detector P1 bkg %d-th peak",i+1),Form("HPGe detector P1 bkg %d-th peak",i+1), bins_HPGe[i], range_HPGe[i][0], range_HPGe[i][1]);
    }


	// histogram filling
    for (int i=0; i<n_peaks; ++i) {
        for (int j=0; j<inbranch_1_1->GetEntries(); j++) {
            inbranch_1_1->GetEntry(j);
            h_1_1[i]->Fill(indata_1_1.qlong);
        }
        for (int j=0; j<inbranch_1_b->GetEntries(); j++) {
            inbranch_1_b->GetEntry(j);
            h_1_b[i]->Fill(indata_1_b.qlong);
        }
    }


    // background subtraction
    for (int i=0; i<n_peaks; ++i) {
        h_1_1[i]->Add(h_1_b[i],-0.4);
    }


    // subtract Compton
    float fit_range_c_HPGe[n_peaks][2] = {
        {    425,    460    }, // 1
        {    840,    880    }, // 2
        {   1175,   1220    }, // 3
        {   1250,   1300    }, // 4
        {   1400,   1440    }, // 5        
        {   1505,   1540    }, // 6
        {   2635,   2670    }, // 7
        {   2935,   2980    }, // 8
        {   3260,   3300    }, // 9
        {   3680,   3710    }, // 10
        {   3760,   3800    }, // 11
        {   4760,   4800    }  // 12
    };
    TF1** f_1_1_c = new TF1*[n_peaks];

    for (int i=0; i<n_peaks; ++i) {
        f_1_1_c[i] = new TF1(Form("fit_1_1_c_%d",i+1),"pol1",range_HPGe[i][0],range_HPGe[i][1]);
    }
    for (int i=0; i<n_peaks; ++i) {
        h_1_1[i]->Fit(f_1_1_c[i],"","",fit_range_c_HPGe[i][0],fit_range_c_HPGe[i][1]);
        h_1_1[i]->Add(f_1_1_c[i],-1.0);
    }


    // gaussian fit functions
    TF1** f_1_1 = new TF1*[n_peaks];
    
    for (int i=0; i<n_peaks; ++i) {
        f_1_1[i] = new TF1(Form("fit_1_1_%d",i+1),"gaus",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
    }


    // draw histograms
    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Position 1", 800, 400*n_peaks);
	    c1->Divide(4,3);
        for (int i=0; i<n_peaks; ++i) {
            c1->cd(i+1);
            h_1_1[i]->Draw();
            h_1_1[i]->Fit(f_1_1[i],"","",fit_range_HPGe[i][0],fit_range_HPGe[i][1]);
        }
    }


    for (int i=0; i<n_peaks; ++i) {
        peak_1_1_Eu[i] = f_1_1[i]->GetParameter(1);
        err_peak_1_1_Eu[i] = f_1_1[i]->GetParError(1);
    }


    for (int i=0; i<n_peaks; ++i) {
        R_1_1[i] = 2.355 * f_1_1[i]->GetParameter(2) / f_1_1[i]->GetParameter(1);
        err_R_1_1[i] = R_1_1[i] * sqrt( pow(f_1_1[i]->GetParError(2)/f_1_1[i]->GetParameter(2),2.0) + pow(f_1_1[i]->GetParError(1)/f_1_1[i]->GetParameter(1),2.0) );
    }


    // compute efficiencies
    // tot_eff = efficiency with total spectrum
    // peak_eff = efficiency with only the peak
    // rel_peak_eff = relative efficiency with only the peak
    float tot_eff;
    float peak_eff[n_peaks];
    float err_peak_eff[n_peaks];

    float rel_peak_eff[n_peaks];
    float err_rel_peak_eff[n_peaks];

    float tot_gamma_on_HPG = 4.68741e5;

    float log_peak_eff[n_peaks];
    float err_log_peak_eff[n_peaks];

    float normalization = 100.0 / h_1_1[peak_1408_index]->Integral();
    float err_normalization = (100 * sqrt(h_1_1[peak_1408_index]->Integral())) / pow(h_1_1[peak_1408_index]->Integral(),2.0);

    tot_eff = getHisto("DATA/Day1/Eu152_calibration_1.root",1,9950,50,10000,0)->Integral() / tot_gamma_on_HPG;
    
    float tot_integral = 0.0;
    for (int i=0; i<n_peaks; ++i) {
        tot_integral += h_1_1[i]->Integral();

        peak_eff[i] = h_1_1[i]->Integral() / (tot_gamma_on_HPG*weights[i]);
        err_peak_eff[i] = sqrt(h_1_1[i]->Integral()) / (tot_gamma_on_HPG*weights[i]);

        log_peak_eff[i] = log(peak_eff[i]);
        err_log_peak_eff[i] = err_peak_eff[i] / peak_eff[i];

        rel_peak_eff[i] = h_1_1[i]->Integral() * normalization / rel_weights[i];
        err_rel_peak_eff[i] = rel_peak_eff[i] * sqrt( 1.0/h_1_1[i]->Integral() + pow(err_normalization/normalization, 2.0));
    }
    cout << "Tot integral: " << tot_integral << endl;

    TCanvas* c2 = new TCanvas("c2","Intrinsic efficiency",800,600);
    TGraphErrors* g_e = new TGraphErrors(12,E_peak_keV,peak_eff,err_E_peak_keV,err_peak_eff);
    c2->cd(1);
    g_e->Draw("AP");
    g_e->SetMarkerStyle(6);
    g_e->SetMarkerSize(2);
    g_e->SetMarkerColor(4);
    g_e->SetTitle("");
    g_e->GetXaxis()->SetTitle("Energy [keV]");
    g_e->GetYaxis()->SetTitle("#varepsilon_{intrinsic}");
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/intrinsic_efficiency_curve_Eu_HPG.pdf");

    TF1* f_e_HPG = new TF1("f_e_HPG","[0]*(exp((x/1000.0)*[1]))",0,1500);

    TCanvas* c3 = new TCanvas("c3","Relative efficiency",800,600);
    TGraphErrors* g_e_rel = new TGraphErrors(12,E_peak_keV,rel_peak_eff,err_E_peak_keV,err_rel_peak_eff);
    c3->cd(1);
    g_e_rel->Draw("AP");
    g_e_rel->Fit(f_e_HPG,"","",0,1500);
    g_e_rel->SetMarkerStyle(6);
    g_e_rel->SetMarkerSize(2);
    g_e_rel->SetMarkerColor(4);
    g_e_rel->SetTitle("");
    g_e_rel->GetXaxis()->SetTitle("Energy [keV]");
    g_e_rel->GetYaxis()->SetTitle("#varepsilon_{rel}");
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/relative_efficiency_curve_Eu_HPG.pdf");


    // Calibration
    TCanvas* c4 = new TCanvas("c4","Calibration",800,600);
    TF1* f_1_1_Eu = new TF1("Calibration Eu HPGe fit","pol1",0,5000);
    TGraphErrors* g_1_1_Eu = new TGraphErrors(n_peaks, peak_1_1_Eu, E_peak_keV, err_peak_1_1_Eu, err_E_peak_keV);
    c4->cd(1);
    g_1_1_Eu->Draw("AP");
    g_1_1_Eu->SetMarkerStyle(6);
    g_1_1_Eu->SetMarkerSize(2);
    g_1_1_Eu->SetTitle("");
    g_1_1_Eu->Fit(f_1_1_Eu, "", "", 0, 5000);
    g_1_1_Eu->GetXaxis()->SetTitle("Energy [a.u.]");
    g_1_1_Eu->GetYaxis()->SetTitle("Energy [keV]");
    g_1_1_Eu->GetXaxis()->SetRangeUser(0,5000);
    g_1_1_Eu->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/calibration_Eu_HPG.pdf");


    // Residuals
    float residuals[n_peaks];
    float err_residuals[n_peaks];
    for (int i=0; i<n_peaks; ++i) {
        residuals[i] = (E_peak_keV[i] - f_1_1_Eu->GetParameter(0) - f_1_1_Eu->GetParameter(1)*peak_1_1_Eu[i]);
        residuals[i] = pow( residuals[i], 2.0 );
        residuals[i] = residuals[i] / pow( f_1_1_Eu->GetParameter(1) * err_peak_1_1_Eu[i], 2.0 );
        err_residuals[i] = 1.0;
    }
    TCanvas* c5 = new TCanvas("c5","Residuals",800,600);
    TF1* f_1_1_zero = new TF1("Zero line","0.0*[0]",0,5000);
    TGraphErrors* g_1_1_Eu_res = new TGraphErrors(n_peaks, peak_1_1_Eu, residuals, err_peak_1_1_Eu, err_residuals);
    g_1_1_Eu_res->Draw("AP");
    g_1_1_Eu_res->Fit(f_1_1_zero, "", "", 0, 5000);
    g_1_1_Eu_res->SetMarkerStyle(7);
    g_1_1_Eu_res->SetMarkerSize(2);
    g_1_1_Eu_res->SetTitle("");
    g_1_1_Eu_res->GetXaxis()->SetTitle("Energy [a.u.]");
    g_1_1_Eu_res->GetYaxis()->SetTitle("Normalized residuals");
    g_1_1_Eu_res->GetXaxis()->SetRangeUser(0,5000);
    g_1_1_Eu_res->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/residuals_calibration_Eu_HPG.pdf");



    ofstream ofile("Analysis/Day1/Day1_Eu152_HPG.txt");
    ofile << "Parameters of gaussian fits for D2P1: True energy [0] err_[0] [1] err_[1] [2] err_[2] R err_R:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << E_peak_keV[i] << '\t' << '$' << '\t'   << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(0) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(0) << '\t' << '$' << '\t'
                                                                                                                << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(1) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(1) << '\t' << '$' << '\t'
                                                                                                                << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParameter(2) << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << f_1_1[i]->GetParError(2) << '\t' << '$' << '\t'
                                                                                                                << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << R_1_1[i]                  << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_R_1_1[i]             << '\t' << '$' << '\t'
                                                                                                                << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Efficiencies for D2P1: Efficiency=[0] err_[0] Log_efficiency=[1] err_[1]:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << E_peak_keV[i] << '\t' << '$' << '\t'   << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << peak_eff[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_peak_eff[i] << '\t' << '$' << '\t'
                                                                                                                << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << log_peak_eff[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_log_peak_eff[i] << '\t' << '$' << '\t'
                                                                                                                << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Relative efficiencies for D2P1: Efficiency=[0] err_[0] Log_efficiency=[1] err_[1]:" << endl;
    for (int i=0; i<n_peaks; ++i) {
        ofile << "Peak " << i+1 << '\t' << '&' << '\t' << '$' << '\t' << E_peak_keV[i] << '\t' << '$' << '\t'   << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << rel_peak_eff[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_rel_peak_eff[i] << '\t' << '$' << '\t'
                                                                                                                << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Calibration coefficients for D2P1:" << endl;
    ofile << "Intercept" << endl;
    ofile << f_1_1_Eu->GetParameter(0) << '\t' << f_1_1_Eu->GetParError(0) << endl;
    ofile << "Slope" << endl;
    ofile << f_1_1_Eu->GetParameter(1) << '\t' << f_1_1_Eu->GetParError(1) << endl;
    ofile.close();

}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getPeakAmNa(int draw=1) {
    // variables
    const int n_peaks_Am = 1;
    float bins_NaI_Am[n_peaks_Am] = {
        50
    };
    float bins_HPGe_Am[n_peaks_Am] = {
        70
    };
    float range_NaI_Am[n_peaks_Am][2] = {
        {   130,    280   },
    };
    float range_HPGe_Am[n_peaks_Am][2] = {
        {   150,    220   }
    };
    float fit_range_NaI_Am[n_peaks_Am][2] = {
        {   160,    240   },
    };
    float fit_range_HPGe_Am[n_peaks_Am][2] = {
        {   198,    212   },
    };
	slimport_data_t indata_0_1_Am[n_peaks_Am];
    slimport_data_t indata_1_1_Am[n_peaks_Am];
    slimport_data_t indata_0_2_Am[n_peaks_Am];
    slimport_data_t indata_1_2_Am[n_peaks_Am];

	TFile *infile_1_Am = new TFile("DATA/Day1/Am241_calibration_1.root");
	TTree *intree_1_Am = (TTree*)infile_1_Am->Get("acq_tree_0");
    TFile *infile_2_Am = new TFile("DATA/Day1/Am241_calibration_2.root");
	TTree *intree_2_Am = (TTree*)infile_2_Am->Get("acq_tree_0");

    TBranch** inbranch_0_1_Am = new TBranch*[n_peaks_Am];
    TBranch** inbranch_1_1_Am = new TBranch*[n_peaks_Am];
    inbranch_0_1_Am[0] = new TBranch;
    inbranch_1_1_Am[0] = new TBranch;
	inbranch_0_1_Am[0] = intree_1_Am->GetBranch("acq_ch0");
    inbranch_1_1_Am[0] = intree_1_Am->GetBranch("acq_ch1");
	inbranch_0_1_Am[0]->SetAddress(&indata_0_1_Am[0].timetag);
    inbranch_1_1_Am[0]->SetAddress(&indata_1_1_Am[0].timetag);

    TBranch** inbranch_0_2_Am = new TBranch*[n_peaks_Am];
    TBranch** inbranch_1_2_Am = new TBranch*[n_peaks_Am];
    inbranch_0_2_Am[0] = new TBranch;
    inbranch_1_2_Am[0] = new TBranch;
	inbranch_0_2_Am[0] = intree_2_Am->GetBranch("acq_ch0");
    inbranch_1_2_Am[0] = intree_2_Am->GetBranch("acq_ch1");
	inbranch_0_2_Am[0]->SetAddress(&indata_0_2_Am[0].timetag);
    inbranch_1_2_Am[0]->SetAddress(&indata_1_2_Am[0].timetag);

    TH1F** h_0_1_Am = new TH1F*[n_peaks_Am];
    TH1F** h_1_1_Am = new TH1F*[n_peaks_Am];
	h_0_1_Am[0] = new TH1F("NaI detector P1","NaI detector spectrum peak at 59.5 keV", bins_NaI_Am[0], range_NaI_Am[0][0], range_NaI_Am[0][1]);
    h_1_1_Am[0] = new TH1F("HPGe detector P1","HPGe detector spectrum peak at 59.5 keV", bins_HPGe_Am[0], range_HPGe_Am[0][0], range_HPGe_Am[0][1]);

    TH1F** h_0_2_Am = new TH1F*[n_peaks_Am];
    TH1F** h_1_2_Am = new TH1F*[n_peaks_Am];
	h_0_2_Am[0] = new TH1F("NaI detector P2","NaI detector spectrum peak at 59.5 keV", bins_NaI_Am[0], range_NaI_Am[0][0], range_NaI_Am[0][1]);
    h_1_2_Am[0] = new TH1F("HPGe detector P2","HPGe detector spectrum peak at 59.5 keV", bins_HPGe_Am[0], range_HPGe_Am[0][0], range_HPGe_Am[0][1]);

	// histogram filling
    for (int i=0; i<n_peaks_Am; ++i) {
        for (int j=0; j<inbranch_0_1_Am[i]->GetEntries(); j++) {
            inbranch_0_1_Am[i]->GetEntry(j);
            h_0_1_Am[i]->Fill(indata_0_1_Am[i].qlong);
        }
        for (int j=0; j<inbranch_1_1_Am[i]->GetEntries(); j++) {
            inbranch_1_1_Am[i]->GetEntry(j);
            h_1_1_Am[i]->Fill(indata_1_1_Am[i].qlong);
        }
        for (int j=0; j<inbranch_0_2_Am[i]->GetEntries(); j++) {
            inbranch_0_2_Am[i]->GetEntry(j);
            h_0_2_Am[i]->Fill(indata_0_2_Am[i].qlong);
        }
        for (int j=0; j<inbranch_1_2_Am[i]->GetEntries(); j++) {
            inbranch_1_2_Am[i]->GetEntry(j);
            h_1_2_Am[i]->Fill(indata_1_2_Am[i].qlong);
        }
    }


    TF1** f_0_1_Am = new TF1*[n_peaks_Am];
    TF1** f_1_1_Am = new TF1*[n_peaks_Am];
    TF1** f_0_2_Am = new TF1*[n_peaks_Am];
    TF1** f_1_2_Am = new TF1*[n_peaks_Am];
    
    for (int i=0; i<n_peaks_Am; ++i) {
        f_0_1_Am[i] = new TF1("fit_0_1_Am","gaus",fit_range_NaI_Am[i][0],fit_range_NaI_Am[i][1]);
        f_1_1_Am[i] = new TF1("fit_1_1_Am","gaus",fit_range_HPGe_Am[i][0],fit_range_HPGe_Am[i][1]);
        f_0_2_Am[i] = new TF1("fit_0_2_Am","gaus",fit_range_NaI_Am[i][0],fit_range_NaI_Am[i][1]);
        f_1_2_Am[i] = new TF1("fit_1_2_Am","gaus",fit_range_HPGe_Am[i][0],fit_range_HPGe_Am[i][1]);
    }

    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Position 1", 800, 400*n_peaks_Am);
	    c1->Divide(2,n_peaks_Am);
        for (int i=0; i<n_peaks_Am; ++i) {
            c1->cd(2*i+1);
            h_0_1_Am[i]->Draw();
            h_0_1_Am[i]->Fit(f_0_1_Am[i],"","",fit_range_NaI_Am[i][0],fit_range_NaI_Am[i][1]);
            c1->cd(2*i+2);
            h_1_1_Am[i]->Draw();
            h_1_1_Am[i]->Fit(f_1_1_Am[i],"","",fit_range_HPGe_Am[i][0],fit_range_HPGe_Am[i][1]);
        }
        
        TCanvas* c2 = new TCanvas("c2", "Position 2", 800, 400*n_peaks_Am);
	    c2->Divide(2,n_peaks_Am);
        for (int i=0; i<n_peaks_Am; ++i) {
            c2->cd(2*i+1);
            h_0_2_Am[i]->Draw();
            h_0_2_Am[i]->Fit(f_0_2_Am[i],"","",fit_range_NaI_Am[i][0],fit_range_NaI_Am[i][1]);
            c2->cd(2*i+2);
            h_1_2_Am[i]->Draw();
            h_1_2_Am[i]->Fit(f_1_2_Am[i],"","",fit_range_HPGe_Am[i][0],fit_range_HPGe_Am[i][1]);
        }
    }



    // Na22 ****************************************************************************************
    // Na22 ****************************************************************************************
    // Na22 ****************************************************************************************
    // Na22 ****************************************************************************************
    // Na22 ****************************************************************************************
    // variables
    const int n_peaks_Na = 2;
    float bins_NaI_Na[n_peaks_Na] = {
        50,     50
    };
    float bins_HPGe_Na[n_peaks_Na] = {
        70,     70
    };
    float range_NaI_Na[n_peaks_Na][2] = {
        {   1300,    1800   },
        {   3200,    4200   }
    };
    float range_HPGe_Na[n_peaks_Na][2] = {
        {   1685,    1755   },
        {   4250,    4320   }
    };
    float fit_range_NaI_Na[n_peaks_Na][2] = {
        {   1400,    1650   },
        {   3500,    3900   }
    };
    float fit_range_HPGe_Na[n_peaks_Na][2] = {
        {   1700,    1730   },
        {   4275,    4295   }
    };
	slimport_data_t indata_0_1_Na;
    slimport_data_t indata_1_1_Na;
    slimport_data_t indata_0_2_Na;
    slimport_data_t indata_1_2_Na;

	TFile *infile_1_Na = new TFile("DATA/Day1/Na22_calibration_1.root");
	TTree *intree_1_Na = (TTree*)infile_1_Na->Get("acq_tree_0");
    TFile *infile_2_Na = new TFile("DATA/Day1/Na22_calibration_2.root");
	TTree *intree_2_Na = (TTree*)infile_2_Na->Get("acq_tree_0");

    TBranch* inbranch_0_1_Na = intree_1_Na->GetBranch("acq_ch0");
    TBranch* inbranch_1_1_Na = intree_1_Na->GetBranch("acq_ch1");
    inbranch_0_1_Na->SetAddress(&indata_0_1_Na.timetag);
    inbranch_1_1_Na->SetAddress(&indata_1_1_Na.timetag);

    TBranch* inbranch_0_2_Na = intree_2_Na->GetBranch("acq_ch0");
    TBranch* inbranch_1_2_Na = intree_2_Na->GetBranch("acq_ch1");
    inbranch_0_2_Na->SetAddress(&indata_0_2_Na.timetag);
    inbranch_1_2_Na->SetAddress(&indata_1_2_Na.timetag);


    TH1F** h_0_1_Na = new TH1F*[n_peaks_Na];
    TH1F** h_1_1_Na = new TH1F*[n_peaks_Na];
	h_0_1_Na[0] = new TH1F("NaI detector P1 Pk1","NaI detector spectrum peak at 511 keV", bins_NaI_Na[0], range_NaI_Na[0][0], range_NaI_Na[0][1]);
    h_0_1_Na[1] = new TH1F("NaI detector P1 Pk2","NaI detector spectrum peak at 1275 keV", bins_NaI_Na[1], range_NaI_Na[1][0], range_NaI_Na[1][1]);
    h_1_1_Na[0] = new TH1F("HPGe detector P1 Pk1","HPGe detector spectrum peak at 511 keV", bins_HPGe_Na[0], range_HPGe_Na[0][0], range_HPGe_Na[0][1]);
    h_1_1_Na[1] = new TH1F("HPGe detector P1 Pk2","HPGe detector spectrum peak at 1275 keV", bins_HPGe_Na[1], range_HPGe_Na[1][0], range_HPGe_Na[1][1]);

    TH1F** h_0_2_Na = new TH1F*[n_peaks_Na];
    TH1F** h_1_2_Na = new TH1F*[n_peaks_Na];
	h_0_2_Na[0] = new TH1F("NaI detector P2 Pk1","NaI detector spectrum peak at 511 keV", bins_NaI_Na[0], range_NaI_Na[0][0], range_NaI_Na[0][1]);
    h_0_2_Na[1] = new TH1F("NaI detector P2 Pk2","NaI detector spectrum peak at 1275 keV", bins_NaI_Na[1], range_NaI_Na[1][0], range_NaI_Na[1][1]);
    h_1_2_Na[0] = new TH1F("HPGe detector P2 Pk1","HPGe detector spectrum peak at 511 keV", bins_HPGe_Na[0], range_HPGe_Na[0][0], range_HPGe_Na[0][1]);
    h_1_2_Na[1] = new TH1F("HPGe detector P2 Pk2","HPGe detector spectrum peak at 1275 keV", bins_HPGe_Na[1], range_HPGe_Na[1][0], range_HPGe_Na[1][1]);    


	// histogram filling
    for (int i=0; i<n_peaks_Na; ++i) {
        for (int j=0; j<inbranch_0_1_Na->GetEntries(); j++) {
            inbranch_0_1_Na->GetEntry(j);
            h_0_1_Na[i]->Fill(indata_0_1_Na.qlong);
        }
        for (int j=0; j<inbranch_1_1_Na->GetEntries(); j++) {
            inbranch_1_1_Na->GetEntry(j);
            h_1_1_Na[i]->Fill(indata_1_1_Na.qlong);
        }
        for (int j=0; j<inbranch_0_2_Na->GetEntries(); j++) {
            inbranch_0_2_Na->GetEntry(j);
            h_0_2_Na[i]->Fill(indata_0_2_Na.qlong);
        }
        for (int j=0; j<inbranch_1_2_Na->GetEntries(); j++) {
            inbranch_1_2_Na->GetEntry(j);
            h_1_2_Na[i]->Fill(indata_1_2_Na.qlong);
        }
    }

    TF1** f_0_1_Na = new TF1*[n_peaks_Na];
    TF1** f_1_1_Na = new TF1*[n_peaks_Na];
    TF1** f_0_2_Na = new TF1*[n_peaks_Na];
    TF1** f_1_2_Na = new TF1*[n_peaks_Na];
    
    for (int i=0; i<n_peaks_Na; ++i) {
        f_0_1_Na[i] = new TF1("fit_0_1_Na","gaus",fit_range_NaI_Na[i][0],fit_range_NaI_Na[i][1]);
        f_1_1_Na[i] = new TF1("fit_1_1_Na","gaus",fit_range_HPGe_Na[i][0],fit_range_HPGe_Na[i][1]);
        f_0_2_Na[i] = new TF1("fit_0_2_Na","gaus",fit_range_NaI_Na[i][0],fit_range_NaI_Na[i][1]);
        f_1_2_Na[i] = new TF1("fit_1_2_Na","gaus",fit_range_HPGe_Na[i][0],fit_range_HPGe_Na[i][1]);
    }


    if (draw==1) {
        TCanvas* c3 = new TCanvas("c3", "Position 1 Na", 800, 400*n_peaks_Na);
	    c3->Divide(2,n_peaks_Na);
        for (int i=0; i<n_peaks_Na; ++i) {
            c3->cd(2*i+1);
            h_0_1_Na[i]->Draw();
            h_0_1_Na[i]->Fit(f_0_1_Na[i],"","",fit_range_NaI_Na[i][0],fit_range_NaI_Na[i][1]);
            c3->cd(2*i+2);
            h_1_1_Na[i]->Draw();
            h_1_1_Na[i]->Fit(f_1_1_Na[i],"","",fit_range_HPGe_Na[i][0],fit_range_HPGe_Na[i][1]);
        }
        
        TCanvas* c4 = new TCanvas("c4", "Position 2 Na", 800, 400*n_peaks_Na);
	    c4->Divide(2,n_peaks_Na);
        for (int i=0; i<n_peaks_Na; ++i) {
            c4->cd(2*i+1);
            h_0_2_Na[i]->Draw();
            h_0_2_Na[i]->Fit(f_0_2_Na[i],"","",fit_range_NaI_Na[i][0],fit_range_NaI_Na[i][1]);
            c4->cd(2*i+2);
            h_1_2_Na[i]->Draw();
            h_1_2_Na[i]->Fit(f_1_2_Na[i],"","",fit_range_HPGe_Na[i][0],fit_range_HPGe_Na[i][1]);
        }
    }


    // Efficiency
    float tot_eff_Am_HPG;
    float tot_eff_Na_HPG;
    float tot_eff_Am_NaI;
    float tot_eff_Na_NaI;
    float peak_eff_NaI[n_peaks_Am+n_peaks_Na];
    float peak_eff_HPG[n_peaks_Am+n_peaks_Na];
    float err_peak_eff_NaI[n_peaks_Am+n_peaks_Na];
    float err_peak_eff_HPG[n_peaks_Am+n_peaks_Na];

    float rel_peak_eff_NaI[n_peaks_Am+n_peaks_Na];
    float rel_peak_eff_HPG[n_peaks_Am+n_peaks_Na];
    float err_rel_peak_eff_NaI[n_peaks_Am+n_peaks_Na];
    float err_rel_peak_eff_HPG[n_peaks_Am+n_peaks_Na];

    float tot_gamma_on_NaI_Am = 1.837e6;
    float tot_gamma_on_NaI_Na = 4.665e4;
    float tot_gamma_on_HPG_Am = 9.6e5;
    float tot_gamma_on_HPG_Na = 2.4e3;

    float weights[3] = {
        0.36,   //1.0000,
        1.81,   //0.6439,
        1.00    //0.3561
    };
    float E_peak_keV[3] = {
        59.5,
        511.0,
        1275.0
    };
    float err_E_peak_keV[3] = {
        0.001,
        0.001,
        0.001
    };

    tot_eff_Am_HPG = getHisto("DATA/Day1/Am241_calibration_1.root",1,9950,50,10000,0)->Integral() / tot_gamma_on_HPG_Am;
    tot_eff_Na_HPG = getHisto("DATA/Day1/Na22_calibration_1.root" ,1,9950,50,10000,0)->Integral() / tot_gamma_on_HPG_Na;
    tot_eff_Am_NaI = getHisto("DATA/Day1/Am241_calibration_2.root",0,9950,50,10000,0)->Integral() / tot_gamma_on_NaI_Am;
    tot_eff_Na_NaI = getHisto("DATA/Day1/Na22_calibration_2.root" ,0,9950,50,10000,0)->Integral() / tot_gamma_on_NaI_Na;
    
    peak_eff_NaI[0] = (h_0_2_Am[0]->Integral()) / (tot_gamma_on_NaI_Am*weights[0]); // /h_0_2_Am[0]->GetBinWidth(1)
    peak_eff_NaI[1] = (h_0_2_Na[0]->Integral()) / (tot_gamma_on_NaI_Na*weights[1]); // /h_0_2_Na[0]->GetBinWidth(1)
    peak_eff_NaI[2] = (h_0_2_Na[1]->Integral()) / (tot_gamma_on_NaI_Na*weights[2]); // /h_0_2_Na[1]->GetBinWidth(1)
    cout << h_0_2_Na[1]->GetEntries() << endl;
    //
    err_peak_eff_NaI[0] = sqrt(h_0_2_Am[0]->Integral()) / (tot_gamma_on_NaI_Am*weights[0]);
    err_peak_eff_NaI[1] = sqrt(h_0_2_Na[0]->Integral()) / (tot_gamma_on_NaI_Na*weights[1]);
    err_peak_eff_NaI[2] = sqrt(h_0_2_Na[1]->Integral()) / (tot_gamma_on_NaI_Na*weights[2]);
    //
    peak_eff_HPG[0] = (h_1_1_Am[0]->Integral()) / (tot_gamma_on_HPG_Am*weights[0]);
    peak_eff_HPG[1] = (h_1_1_Na[0]->Integral()) / (tot_gamma_on_HPG_Na*weights[1]);
    peak_eff_HPG[2] = (h_1_1_Na[1]->Integral()) / (tot_gamma_on_HPG_Na*weights[2]);
    //
    err_peak_eff_HPG[0] = sqrt(h_1_1_Am[0]->Integral()) / (tot_gamma_on_HPG_Am*weights[0]);
    err_peak_eff_HPG[1] = sqrt(h_1_1_Na[0]->Integral()) / (tot_gamma_on_HPG_Na*weights[1]);
    err_peak_eff_HPG[2] = sqrt(h_1_1_Na[1]->Integral()) / (tot_gamma_on_HPG_Na*weights[2]);

    TF1* f_e_NaI = new TF1("f_e_NaI","[0]*(exp(sqrt(x/1000.0)*(10*[1])))",0,1400);
    TF1* f_e_HPG = new TF1("f_e_HPG","[0]*(exp(sqrt(x/1000.0)*(10*[1])))",0,1400);

    peak_eff_NaI[0] = 0.31423;
    peak_eff_NaI[1] = 0.125425;

    peak_eff_HPG[1] = 0.04347;
    peak_eff_HPG[0] = peak_eff_HPG[0] / 8;
    peak_eff_HPG[2] = peak_eff_HPG[2] / 8;

    for (int i=0; i<3; ++i) {
        err_peak_eff_HPG[i] = err_peak_eff_HPG[i] / 8;
    }
    err_peak_eff_HPG[1] = err_peak_eff_HPG[1] / 4;
    err_peak_eff_HPG[2] = err_peak_eff_HPG[2] / 4;


    TCanvas* c5 = new TCanvas("c5","NaI Intrinsic efficiency",800,600);
    TGraphErrors* g_e_NaI = new TGraphErrors(3,E_peak_keV,peak_eff_NaI,err_E_peak_keV,err_peak_eff_NaI);
    g_e_NaI->Draw("AP");
    g_e_NaI->Fit(f_e_NaI,"","",0,1400);
    g_e_NaI->SetMarkerStyle(6);
    g_e_NaI->SetMarkerSize(2);
    g_e_NaI->SetMarkerColor(4);
    g_e_NaI->SetTitle("");
    g_e_NaI->GetXaxis()->SetTitle("Energy [keV]");
    g_e_NaI->GetYaxis()->SetTitle("#varepsilon_{intrinsic}");
    g_e_NaI->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/intrinsic_efficiency_curve_Am_Na_NaI.pdf");

    TCanvas* c6 = new TCanvas("c6","HPG Intrinsic efficiency",800,600);
    TGraphErrors* g_e_HPG = new TGraphErrors(3,E_peak_keV,peak_eff_HPG,err_E_peak_keV,err_peak_eff_HPG);
    g_e_HPG->Draw("AP");
    g_e_HPG->Fit(f_e_HPG,"","",0,1400);
    g_e_HPG->SetMarkerStyle(6);
    g_e_HPG->SetMarkerSize(2);
    g_e_HPG->SetMarkerColor(4);
    g_e_HPG->SetTitle("");
    g_e_HPG->GetXaxis()->SetTitle("Energy [keV]");
    g_e_HPG->GetYaxis()->SetTitle("#varepsilon_{intrinsic}");
    g_e_HPG->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/intrinsic_efficiency_curve_Am_Na_HPG.pdf");

    ofstream ofile("Analysis/Day1/Day1_efficiency_Am_Na_NaI.txt");
    ofile << "Efficiency for NaI:" << endl;
    for (int i=0; i<3; ++i) {
        ofile << "NaI " << '\t' << '&' << '\t' << '$' << '\t'   << fixed << setw(12) << setprecision(8) << peak_eff_NaI[i] << '\t' << "\\pm" << '\t' 
                                                                << fixed << setw(12) << setprecision(8) << err_peak_eff_NaI[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Efficiency for HPG:" << endl;
    for (int i=0; i<3; ++i) {
        ofile << "HPG " << '\t' << '&' << '\t' << '$' << '\t'   << fixed << setw(12) << setprecision(8) << peak_eff_HPG[i] << '\t' << "\\pm" << '\t' 
                                                                << fixed << setw(12) << setprecision(8) << err_peak_eff_HPG[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "Efficiency fit parameters for NaI:" << endl;
    ofile << "    [0]" << '\t' << f_e_NaI->GetParameter(0) << endl;
    ofile << "err_[0]" << '\t' << f_e_NaI->GetParError(0) << endl;
    ofile << "    [1]" << '\t' << f_e_NaI->GetParameter(1) << endl;
    ofile << "err_[1]" << '\t' << f_e_NaI->GetParError(1) << endl;
    ofile << endl << endl;
    ofile << "Efficiency fit parameters for HPGe:" << endl;
    ofile << "    [0]" << '\t' << f_e_HPG->GetParameter(0) << endl;
    ofile << "err_[0]" << '\t' << f_e_HPG->GetParError(0) << endl;
    ofile << "    [1]" << '\t' << f_e_HPG->GetParameter(1) << endl;
    ofile << "err_[1]" << '\t' << f_e_HPG->GetParError(1) << endl;
    ofile.close();
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void calibration_Am_Na() {
    float x[4][3] = {
        {   2.02051e+002,   1.52785e+003,   3.69419e+003   },
        {   2.01991e+002,   1.71709e+003,   4.28318e+003   },
        {   2.03159e+002,   1.51741e+003,   3.66935e+003   },
        {   2.02005e+002,   1.71698e+003,   4.28356e+003   }
    };
    float err_x[4][3] = {
        {   3.45620e-002,   6.17599e-001,   2.32320e+000    },
        {   1.37152e-002,   1.33161e-001,   2.58480e-001    },
        {   5.18130e-002,   1.14906e+000,   4.93543e+000    },
        {   1.36268e-002,   1.33981e-001,   2.58866e-001    }
    };
    float y[3] = {
        59.5409,
        511.0,
        1274.537
    };
    float err_y[3] = {
        0.0001,
        0.000001,
        0.007
    };


    TGraphErrors** gr = new TGraphErrors*[2];
    gr[0] = new TGraphErrors(3,x[3],y,err_x[3],err_y);  // NaI
    gr[1] = new TGraphErrors(3,x[0],y,err_x[0],err_y);  // HPGe

    TF1** f = new TF1*[2];
    f[0] = new TF1("Calibration NaI","[0] + [1]*x", 0, 8000);
    f[1] = new TF1("Calibration HPGe","[0] + [1]*x", 0, 8000);

    gr[0]->Draw("AP");
    gr[0]->SetMarkerStyle(6);
    gr[0]->SetMarkerSize(2);
    gr[0]->SetTitle("");
    gr[0]->Fit(f[0],"","",0,5000);
    gr[0]->GetXaxis()->SetTitle("Energy [a.u.]");
    gr[0]->GetYaxis()->SetTitle("Energy [keV]");
    gr[0]->GetXaxis()->SetRangeUser(0,5000);
    gr[0]->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/calibration_Am_Na_NaI.pdf");

    TCanvas* c2 = new TCanvas("c2", "Calibration HPGe", 800, 600);
    gr[1]->Draw("AP");
    gr[1]->SetMarkerStyle(6);
    gr[1]->SetMarkerSize(2);
    gr[1]->SetTitle("");
    gr[1]->Fit(f[1],"","",0,5000);
    gr[1]->GetXaxis()->SetTitle("Energy [a.u.]");
    gr[1]->GetYaxis()->SetTitle("Energy [keV]");
    gr[1]->GetXaxis()->SetRangeUser(0,5000);
    gr[1]->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/calibration_Am_Na_HPG.pdf");


    // residuals
    // 0 -> NaI  -> x[3]
    // 1 -> HPGe -> x[0]
    float residuals[2][3];
    float err_residuals[2][3];
    for (int i=0; i<3; ++i) {
        residuals[0][i] = (y[i] - f[0]->GetParameter(0) - f[0]->GetParameter(1)*x[3][i]);
        residuals[1][i] = (y[i] - f[1]->GetParameter(0) - f[1]->GetParameter(1)*x[0][i]);
        residuals[0][i] = pow( residuals[0][i], 2.0 );
        residuals[1][i] = pow( residuals[1][i], 2.0 );
        residuals[0][i] = residuals[0][i] / pow( f[0]->GetParameter(1) * err_x[3][i], 2.0 );
        residuals[1][i] = residuals[1][i] / pow( f[1]->GetParameter(1) * err_x[0][i], 2.0 );
        err_residuals[0][i] = 1.0;
        err_residuals[0][i] = 1.0;
    }

    TCanvas* c3 = new TCanvas("c3","Residuals NaI for Na-Am calibration",800,600);
    TF1* f_NaI_zero = new TF1("Zero line NaI","0.0*[0]",0,5000);
    TGraphErrors* g_NaI_res = new TGraphErrors(3, x[3], residuals[0], err_x[3], err_residuals[0]);
    g_NaI_res->Draw("AP");
    g_NaI_res->Fit(f_NaI_zero, "", "", 0, 5000);
    g_NaI_res->SetMarkerStyle(7);
    g_NaI_res->SetMarkerSize(2);
    g_NaI_res->SetTitle("");
    g_NaI_res->GetXaxis()->SetTitle("Energy [a.u.]");
    g_NaI_res->GetYaxis()->SetTitle("Normalized residuals");
    g_NaI_res->GetXaxis()->SetRangeUser(0,5000);
    g_NaI_res->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/residuals_calibration_Am_Na_NaI.pdf");

    TCanvas* c4 = new TCanvas("c3","Residuals HPGe for Na-Am calibration",800,600);
    TF1* f_HPG_zero = new TF1("Zero line HPGe","0.0*[0]",0,5000);
    TGraphErrors* g_HPG_res = new TGraphErrors(3, x[0], residuals[1], err_x[0], err_residuals[1]);
    g_HPG_res->Draw("AP");
    g_HPG_res->Fit(f_HPG_zero, "", "", 0, 5000);
    g_HPG_res->SetMarkerStyle(7);
    g_HPG_res->SetMarkerSize(2);
    g_HPG_res->SetTitle("");
    g_HPG_res->GetXaxis()->SetTitle("Energy [a.u.]");
    g_HPG_res->GetYaxis()->SetTitle("Normalized residuals");
    g_HPG_res->GetXaxis()->SetRangeUser(0,5000);
    g_HPG_res->GetYaxis()->SetTitleOffset(1.3);
    gPad->Update();
    gPad->SaveAs("Analysis/Day1/residuals_calibration_Am_Na_HPG.pdf");



    float a[2][1] = {
        {   f[0]->GetParameter(0)   },
        {   f[1]->GetParameter(0)   }
    };
    float err_a[2][1] = {
        {   f[0]->GetParError(0)   },
        {   f[1]->GetParError(0)   }
    };
    float b[2][1] = {
        {   f[0]->GetParameter(1)   },
        {   f[1]->GetParameter(1)   }
    };
    float err_b[2][1] = {
        {   f[0]->GetParError(1)   },
        {   f[1]->GetParError(1)   }
    };


    ofstream ofile("Analysis/Day1/Day1_calibration.txt");
    ofile << "Calibration : a err_a b err_b" << endl;
    ofile << "NaI " << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << a[0][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_a[0][0] << '\t' << '$' << '\t'
                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << b[0][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_b[0][0] << '\t' << '$' << '\t' << "\\\\" << endl;
    ofile << "HPGe" << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << a[1][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_a[1][0] << '\t' << '$' << '\t'
                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << b[1][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_b[1][0] << '\t' << '$' << '\t' << "\\\\" << endl;
    ofile.close();
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void calibration_Eu() {
    float x[2][7] = {
        {
            387.11970222,
            740.19256355,
            1036.23474052,
            2287.21115773,
            2803.15681530,
            3190.03618752,
            4056.13744724
        },
        {
            411.28118196,
            823.98375284,
            1157.30683762,
            2618.57497960,
            3240.39780676,
            3737.67048690,
            4730.94710846,
        }
    };
    float err_x[2][7] = {
        {
            0.05289244,
            0.31179264,
            0.14491640,
            0.86715819,
            1.06132261,
            1.05361983,
            0.91331362
        },
        {
            0.02669345,
            0.07367781,
            0.03432761,
            0.12324705,
            0.08411860,
            0.22820608,
            0.06401013
        }
    };
    float y[7] = {
        121.8,
        244.7,
        344.3,
        778.9,
        964.0,
        1112.1,
        1408.0
    };
    float err_y[7] = {
        0.0001,
        0.0001,
        0.0001,
        0.0001,
        0.0001,
        0.0001,
        0.0001
    };


    TGraphErrors** gr = new TGraphErrors*[2];
    gr[0] = new TGraphErrors(7,x[0],y,err_x[0],err_y);  // NaI
    gr[1] = new TGraphErrors(7,x[1],y,err_x[1],err_y);  // HPGe

    TF1** f = new TF1*[2];
    f[0] = new TF1("Calibration NaI","[0] + [1]*x", 0, 8000);
    f[1] = new TF1("Calibration HPGe","[0] + [1]*x", 0, 8000);

    gr[0]->Fit(f[0],"","",0,5000);
    gr[1]->Fit(f[1],"","",0,5000);

    TCanvas* c1 = new TCanvas("c1", "Calibration", 800, 400);
	c1->Divide(2,1);
    c1->cd(1);
    gr[0]->Draw("AP");
    f[0]->Draw("SAME");
    f[0]->SetLineWidth(0.5);
    gr[0]->GetXaxis()->SetTitle("Energy [a.u.]");
    gr[0]->GetYaxis()->SetTitle("Energy [keV]");
    gr[0]->GetYaxis()->SetTitleOffset(1.6);
    gr[0]->SetTitle("Calibration of detector NaI(Tl)");

    c1->cd(2);
    gr[1]->Draw("AP");
    f[1]->Draw("SAME");
    f[1]->SetLineWidth(0.5);
    gr[1]->GetXaxis()->SetTitle("Energy [a.u.]");
    gr[1]->GetYaxis()->SetTitle("Energy [keV]");
    gr[1]->GetYaxis()->SetTitleOffset(1.6);
    gr[1]->SetTitle("Calibration of detector HPGe");

    float a[2][1] = {
        {   f[0]->GetParameter(0)   },
        {   f[1]->GetParameter(0)   }
    };
    float err_a[2][1] = {
        {   f[0]->GetParError(0)   },
        {   f[1]->GetParError(0)   }
    };
    float b[2][1] = {
        {   f[0]->GetParameter(1)   },
        {   f[1]->GetParameter(1)   }
    };
    float err_b[2][1] = {
        {   f[0]->GetParError(1)   },
        {   f[1]->GetParError(1)   }
    };

    ofstream ofile("Analysis/Day1/Day1_calibration_Eu.txt");
    ofile << "Calibration : a err_a b err_b" << endl;
    ofile << "NaI " << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << a[0][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_a[0][0] << '\t' << '$' << '\t'
                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << b[0][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_b[0][0] << '\t' << '$' << '\t' << "\\\\" << endl;
    ofile << "HPGe" << '\t' << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << a[1][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_a[1][0] << '\t' << '$' << '\t'
                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(8) << b[1][0] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(8) << err_b[1][0] << '\t' << '$' << '\t' << "\\\\" << endl;
    ofile.close();
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getPeakEuCalibrated() {
    float m_NaI = 2.97788e-001;// 0.3531;      //2.97788e-001;
    float m_HPG = 3.42994e-001;// 0.3550;      //3.42994e-001;
    float q_NaI = -6.11623e-001; // -28.4603;   //-6.11623e-001;
    float q_HPG = -9.76775e+000; // -27.9182;   //-9.76775e+000;

    float m_NaI_Eu = 3.47688e-001;
    float m_HPG_Eu = 2.97735e-001;
    float q_NaI_Eu = -1.30981e+001;
    float q_HPG_Eu = -5.41292e-001;

    TH1F* h_NaI = calibrateHisto("DATA/Day1/Eu152_calibration_1.root",0,900,0,1800,m_NaI_Eu,q_NaI_Eu);
    TH1F* h_HPG = calibrateHisto("DATA/Day1/Eu152_calibration_1.root",1,900,0,1800,m_HPG_Eu,q_HPG_Eu);

    
    TCanvas* c1 = new TCanvas("c1", "Calibrated spectra", 800, 400);
	c1->Divide(2,1);
    c1->cd(1);
    h_NaI->Draw();
    c1->cd(2);
    h_HPG->Draw();
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* getHisto(char* file_name, short chan, int numBins, double minX, double maxX, int draw=1) {
    // variables
	slimport_data_t indata;
	TFile *infile = new TFile(file_name);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1F *h = new TH1F("h","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h->Fill(indata.qlong);
	}
    if (draw==1) h->Draw();

    return h;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* getHistoNoBkg(char* file_name_0, char* file_name_1, short chan, int numBins, double minX, double maxX, int draw=1) {
    // variables
	slimport_data_t indata_0;
    slimport_data_t indata_1;   // bkg
	TFile *infile_0 = new TFile(file_name_0);
    TFile *infile_1 = new TFile(file_name_1);
	TTree *intree_0 = (TTree*)infile_0->Get("acq_tree_0");
    TTree *intree_1 = (TTree*)infile_1->Get("acq_tree_0");
	TBranch *inbranch_0 = intree_0->GetBranch(Form("acq_ch%d",chan));
    TBranch *inbranch_1 = intree_1->GetBranch(Form("acq_ch%d",chan));
	inbranch_0->SetAddress(&indata_0.timetag);
    inbranch_1->SetAddress(&indata_1.timetag);
	TH1F *h_0 = new TH1F("h","Total spectrum",numBins,minX,maxX);
    TH1F *h_1 = new TH1F("h","Background spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch_0->GetEntries(); i++) {
		inbranch_0->GetEntry(i);
		h_0->Fill(indata_0.qlong);
        if (i==inbranch_0->GetEntries()-1) {
            float Dt_0 = indata_0.timetag*1e-8;
            cout << Dt_0 << endl;
        }
	}
    for (int i=0; i<inbranch_1->GetEntries(); i++) {
		inbranch_1->GetEntry(i);
		h_1->Fill(indata_1.qlong);
        if (i==inbranch_1->GetEntries()-1) {
            float Dt_1 = indata_1.timetag*1e-8;
            cout << Dt_1 << endl;
        }
	}
    // histogram subtraction
    float fraction = - Dt_0/Dt_1;
    h_0->Add(h_1,fraction);
    // histogram drawing
    if (draw==1) h_0->Draw();

    return h_0;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* calibrateHisto(char* file_name, short chan, int numBins, double minX, double maxX, float m, float q) {
    // variables
	slimport_data_t indata;
	TFile *infile = new TFile(file_name);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1F *h = new TH1F("h","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h->Fill(indata.qlong*1.0*m + 1.0*q);
	}

	if (m!=1 && q!=0) //This means that I actually changed the calibration!
	    h->SetXTitle("keV");
    
    return h;
    // h->Draw();
}
// ********************************************************************************************
// ********************************************************************************************