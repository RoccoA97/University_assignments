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
float mu_Co60_NaI[2] = {
    2.77605e+003,
    3.14230e+003
    //5.91359e+003
};
float err_mu_Co60_NaI[2] = {
    2.72881e-001,
    3.15863e-001
    //1.47463e+000
};
float sigma_Co60_NaI[2] = {
    8.72154e+001,
    9.31375e+001
    //1.22824e+002
};
float err_sigma_Co60_NaI[2] = {
    3.39030e-001,
    3.64973e-001
    //1.79592e+000
};

float mu_Co60_HPG[2] = {
    3.94247e+003,
    4.47781e+003
};
float err_mu_Co60_HPG[2] = {
    9.52439e-002,
    1.54134e-001
};
float sigma_Co60_HPG[2] = {
    4.01395e+000,
    3.89286e+000
};
float err_sigma_Co60_HPG[2] = {
    7.31152e-002,
    9.38602e-002
};

// ReCalibration coefficients
float m_Co60_NaI = (1332.508 - 1173.24) / (mu_Co60_NaI[1] - mu_Co60_NaI[0]);
float m_Co60_HPG = (1332.508 - 1173.24) / (mu_Co60_HPG[1] - mu_Co60_HPG[0]);
float q_Co60_NaI = 1173.24 - m_Co60_NaI*mu_Co60_NaI[0];
float q_Co60_HPG = 1173.24 - m_Co60_HPG*mu_Co60_HPG[0];





// ********************************************************************************************
// ********************************************************************************************
void getPeakAutunite() {
    const int n_NaI = 7;
    float mu_NaI[n_NaI] = {
        1.73975e+002,
        2.32588e+002,
        2.87564e+002,
        3.48202e+002,
        6.17580e+002,
        1.14835e+003,
        1.79068e+003
    };
    float err_mu_NaI[n_NaI] = {
        3.68326e-001,
        3.82010e-001,
        2.03853e-001,
        1.39888e-001,
        2.10206e-001,
        1.03159e+000,
        1.04182e+000
    };
    float sigma_NaI[n_NaI] = {
        2.12677e+001,
        2.54614e+001,
        2.12968e+001,
        2.01183e+001,
        2.77109e+001,
        6.65664e+001,
        7.37034e+001
    };
    float err_sigma_NaI[n_NaI] = {
        6.06544e-001,
        7.81909e-001,
        3.45329e-001,
        1.72912e-001,
        2.53615e-001,
        1.69814e+000,
        1.24816e+000
    };
    float counts_NaI[n_NaI] = {
        2.47460000000000000e+004,
        2.71620000000000000e+004,
        3.60490000000000000e+004,
        4.24160000000000000e+004,
        3.34740000000000000e+004,
        1.37420000000000000e+004,
        9.34700000000000000e+003
    };


    const int n_HPG = 9;
    float mu_HPG[n_HPG] = {
        1.86666e+002,
        2.42715e+002,
        2.96029e+002,
        3.52623e+002,
        6.09816e+002,
        7.68652e+002,
        9.34477e+002,
        1.12070e+003,
        1.23830e+003
    };
    float err_mu_HPG[n_HPG] = {
        7.05072e-002,
        6.21664e-002,
        4.35645e-002,
        3.46664e-002,
        3.38083e-002,
        1.41576e-001,
        2.00938e-001,
        8.75326e-002,
        1.21098e-001
    };
    float sigma_HPG[n_HPG] = {
        2.43853e+000,
        1.76816e+000,
        1.42866e+000,
        1.58416e+000,
        1.41862e+000,
        2.30793e+000,
        1.85112e+000,
        1.58169e+000,
        1.73632e+000
    };
    float err_sigma_HPG[n_HPG] = {
        8.79204e-002,
        6.40115e-002,
        4.18882e-002,
        4.37958e-002,
        2.29376e-002,
        1.58530e-001,
        2.21418e-001,
        8.50609e-002,
        1.15905e-001
    };
    float counts_HPG[n_HPG] = {
        2.52600000000000000e+003,
        2.05500000000000000e+003,
        3.15800000000000000e+003,
        4.27100000000000000e+003,
        3.13100000000000000e+003,
        4.77000000000000000e+002,
        2.80000000000000000e+002,
        7.11000000000000000e+002,
        2.99000000000000000e+002
    };
    float base_counts_HPG[n_HPG] = {
        {40,    30},
        {30,    30},
        {20,    40},
        {15,    40},
        {5,     50},
        {5,     50},
        {3,     40},
        {2,     50},
        {2,     40}
    };
}





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
        TH1F* h_0_cal = calibrateHistoNaI(file_name_0, 1000, 20, 2020);
        TH1F* h_1_cal = calibrateHistoNaI(file_name_1, 1000, 20, 2020);
    }
    if (chan==1) {
        TH1F* h_0_cal = calibrateHistoHPG(file_name_0, 2000, 20, 2020);
        TH1F* h_1_cal = calibrateHistoHPG(file_name_1, 2000, 20, 2020);
    }
    // histogram subtraction
    float fraction = - Dt_0/Dt_1;
    cout << fraction << endl;
    if (subtract==1) h_0_cal->Add(h_1_cal,fraction);
    // histogram drawing
    if (draw==1) {
        h_0_cal->Draw();
        if (subtract==1) {
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
		h->Fill(indata.qlong*1.0*m + 1.0*q + gRandom->Rndm());
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
TH1F* calibrateHistoNaI(char* file_name, int numBins, double minX, double maxX, float m=m_Co60_NaI, float q=q_Co60_NaI) {
    return calibrateHisto(file_name, 0, numBins, minX, maxX, m, q);
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* calibrateHistoHPG(char* file_name, int numBins, double minX, double maxX, float m=m_Co60_HPG, float q=q_Co60_HPG) {
    return calibrateHisto(file_name, 1, numBins, minX, maxX, m, q);
}
// ********************************************************************************************
// ********************************************************************************************