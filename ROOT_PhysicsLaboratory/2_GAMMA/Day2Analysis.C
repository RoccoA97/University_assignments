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



// ReCalibration
float mu_Na22_NaI[2] = {
    1.27939e+003,
    3.09807e+003
};
float err_mu_Na22_NaI[2] = {
    1.39581e-001,
    5.78429e-001
};
float sigma_Na22_NaI[2] = {
    5.27325e+001,
    9.37685e+001
};
float err_sigma_Na22_NaI[2] = {
    1.38175e-001,
    6.48486e-001
};

float mu_Na22_HPG[2] = {
    1.71793e+003,
    4.28367e+003
};
float err_mu_Na22_HPG[2] = {
    1.29064e-001,
    1.96157e-001
};
float sigma_Na22_HPG[2] = {
    5.00231e+000,
    3.80469e+000
};
float err_sigma_Na22_HPG[2] = {
    8.15608e-002,
    1.46578e-001
};

// ReCalibration coefficients
float m_Na22_NaI = (1274.8-511.0) / (mu_Na22_NaI[1] - mu_Na22_NaI[0]);
float m_Na22_HPG = (1274.8-511.0) / (mu_Na22_HPG[1] - mu_Na22_HPG[0]);
float q_Na22_NaI = 511.0 - m_Na22_NaI*mu_Na22_NaI[0];
float q_Na22_HPG = 511.0 - m_Na22_HPG*mu_Na22_HPG[0];

// Efficiency fit parameters
float a_fit_eff_NaI     =  0.10865053;
float b_fit_eff_NaI     = -3.26162090;
float err_a_fit_eff_NaI =  0.00179308;
float err_b_fit_eff_NaI =  0.05992272;

float a_fit_eff_HPG     =  0.41182833;
float b_fit_eff_HPG     = -2.95312349;
float err_a_fit_eff_HPG =  0.01137090;
float err_b_fit_eff_HPG =  0.11124495;





// ********************************************************************************************
// ********************************************************************************************
void activities() {
    float Dt = 900.0;
    float e_geom_NaI = (1.0/4.0) * pow((37.5/25.0),2.0);
    float e_geom_NaI = (1.0/4.0) * pow((37.5/90.0),2.0);

    float mu_H2O_NaI = 6.61931e+002;
    float err_mu_H2O_NaI = 2.01928e+000;
    float sigma_H2O_NaI = 2.41441e+001;
    float err_sigma_H2O_NaI = 1.42471e+000;
    float int_H2O_NaI = 1.83009902000427250e+002;

    float mu_H2O_HPG = 6.59805e+002;
    float err_mu_H2O_HPG = 4.23423e-001;
    float sigma_H2O_HPG = 2.06843e+000;
    float err_sigma_H2O_HPG = 5.50942e-001;
    float int_H2O_HPG = 4.00004463195800780e+001;

    float mu_KCl_NaI = 6.61931e+002;
    float err_mu_KCl_NaI = 2.01928e+000;
    float sigma_KCl_NaI = 2.41441e+001;
    float err_sigma_KCl_NaI = 1.42471e+000;
    float int_KCl_NaI = 1.83009902000427250e+002;


    float act_NaI[]
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
TH1F* getHistoNoBkg(char* file_name_0, char* file_name_1, short chan, int numBins, double minX, double maxX, int subtract=1, int draw=1) {
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
    // histogram calibration
    if (chan==0) {
        TH1F* h_0_cal = calibrateHistoNaI(file_name_0, 250, 20, 2020);
        TH1F* h_1_cal = calibrateHistoNaI(file_name_1, 250, 20, 2020);
    }
    if (chan==1) {
        TH1F* h_0_cal = calibrateHistoHPG(file_name_0, 1000, 20, 2020);
        TH1F* h_1_cal = calibrateHistoHPG(file_name_1, 1000, 20, 2020);
    }
    // histogram subtraction
    float fraction = - Dt_0/Dt_1;
    cout << fraction << endl;
    if (subtract==1) h_0_cal->Add(h_1_cal,fraction);
    // histogram drawing
    if (draw==1) {
        h_0_cal->Draw();
        if (subtract!=1) {
            h_1_cal->Draw("SAME");
            h_1_cal->SetLineColor(2);
        }
    }

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
	    h->SetXTitle("Energy [keV]");
    
    return h;
    // h->Draw();
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* calibrateHistoNaI(char* file_name, int numBins, double minX, double maxX, float m=m_Na22_NaI, float q=q_Na22_NaI) {
    return calibrateHisto(file_name, 0, numBins, minX, maxX, m, q);
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* calibrateHistoHPG(char* file_name, int numBins, double minX, double maxX, float m=m_Na22_HPG, float q=q_Na22_HPG) {
    return calibrateHisto(file_name, 1, numBins, minX, maxX, m, q);
}
// ********************************************************************************************
// ********************************************************************************************.q