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





// Calibration parameters - ch_0
float q_ch_0 = -6.86330e+000;
float m_ch_0 =  2.64253e-001;
float err_q_ch_0 = 9.08190e-002;
float err_m_ch_0 = 2.44084e-005;

// Calibration parameters - ch_1
float q_ch_1 = -5.89484e+000;
float m_ch_1 =  2.33405e-001;
float err_q_ch_1 = 2.34874e-002;
float err_m_ch_1 = 1.83300e-005;





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
    float t_f;
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h->Fill(indata.qlong);
        if (i==inbranch->GetEntries()-1) t_f = indata.timetag * 4.0e-9;
	}
    if (draw==1) h->Draw();

    cout << t_f << endl;

    return h;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* getHistoCalibrated(char* file_name, short chan, int numBins, double minX, double maxX, int draw=1) {
    // variables
	slimport_data_t indata;
	TFile *infile = new TFile(file_name);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1F *h = new TH1F("h","Total spectrum",numBins,minX,maxX);
	// histogram filling
    if (chan==0) {
        for (int i=0; i<inbranch->GetEntries(); i++) {
		    inbranch->GetEntry(i);
		    h->Fill(indata.qlong*m_ch_0 + q_ch_0);
	    }
    }
    else {
        for (int i=0; i<inbranch->GetEntries(); i++) {
		    inbranch->GetEntry(i);
		    h->Fill(indata.qlong*m_ch_1 + q_ch_1);
	    }
    }
    if (draw==1) h->Draw();

    return h;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH2D* getCoinc(char* file_name, int numBinsX, int numBinsY, double minX, double maxX, double minY, double maxY, int draw=1) {
    // variables
	slimport_data_t indata_0;
    slimport_data_t indata_1;

	TFile *infile = new TFile(file_name);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch_0 = intree->GetBranch("acq_ch0");
    TBranch *inbranch_1 = intree->GetBranch("acq_ch1");
	inbranch_0->SetAddress(&indata_0.timetag);
    inbranch_1->SetAddress(&indata_1.timetag);

	TH2D *h = new TH2D("h","coincidence",numBinsX,minX,maxX,numBinsY,minY,maxY);
	// histogram filling
    cout << "Number of entries: " << inbranch_0->GetEntries() << endl;
	for (int i=0; i<inbranch_0->GetEntries(); i++) {
        // cout << indata_0.timetag * 4 * 1e-9 << '\t' << indata_1.timetag * 4 * 1e-9 << endl;
        inbranch_0->GetEntry(i);
        inbranch_1->GetEntry(i);
        
        if ((indata_0.qlong<maxX) && (indata_0.qlong>minX)) {
            if ((indata_1.qlong<maxY) && (indata_1.qlong>minY)) {
                h->Fill(indata_0.qlong,indata_1.qlong);
            }
        }
	}
    if (draw==1) h->Draw();

    return h;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void getAllCoincHisto(int numBinsX, int numBinsY, double minX, double maxX, double minY, double maxY) {
    // variables
    const int n_files = 12;

    TCanvas* c1 = new TCanvas("c1", "Coincidences", 800*4, 600*3);
    c1->Divide(4,3);
    for (int i=0; i<n_files; ++i) {
        c1->cd(i+1);
        TH2D* h = getCoinc(Form("DATA/Day2/run_%d_logic_or.root",i+1), numBinsX, numBinsY, minX, maxX, minY, maxY, 0);
        h->Draw("COLZ");
        gStyle->SetPalette(1);
        gStyle->SetOptStat(0);
        h->SetContour(10000);
    }
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void coincMethod(int numBinsX, int numBinsY, double minX, double maxX, double minY, double maxY, int draw=1) {
    // variables
    const int n_files = 12;
    float E_keV_0;
    float E_keV_1;
	slimport_data_t indata_0[n_files];
    slimport_data_t indata_1[n_files];
    slimport_data_t indata_b_0;
    slimport_data_t indata_b_1;

    // data retrieving
    TFile** infile = new TFile*[n_files];
    TTree** intree = new TTree*[n_files];
    TBranch** inbranch_0 = new TBranch*[n_files];
    TBranch** inbranch_1 = new TBranch*[n_files];

    for (int i=0; i<n_files; ++i) {
        intree[i] = new TTree;
        inbranch_0[i] = new TBranch;
        inbranch_1[i] = new TBranch;

        infile[i] = new TFile(Form("DATA/Day2/run_%d_logic_or.root",i+1));
        intree[i] = (TTree*)infile[i]->Get("acq_tree_0");
        inbranch_0[i] = intree[i]->GetBranch("acq_ch0");
        inbranch_1[i] = intree[i]->GetBranch("acq_ch1");
        inbranch_0[i]->SetAddress(&indata_0[i].timetag);
        inbranch_1[i]->SetAddress(&indata_1[i].timetag);
    }

	TH2D* h = new TH2D("h","coincidence",numBinsX,minX,maxX,numBinsY,minY,maxY);
	// histogram filling
    for (int i=0; i<n_files; ++i) {
        if (i!=1) {
            for (int j=0; j<inbranch_0[i]->GetEntries(); j++) {
                inbranch_0[i]->GetEntry(j);
                inbranch_1[i]->GetEntry(j);
                
                E_keV_0 = indata_0[i].qlong*m_ch_0 + q_ch_0 + gRandom->Rndm();
                E_keV_1 = indata_1[i].qlong*m_ch_1 + q_ch_1 + gRandom->Rndm();

                if ((E_keV_0<maxX) && (E_keV_0>minX)) {
                    if ((E_keV_1<maxY) && (E_keV_1>minY)) {
                        h->Fill(E_keV_0,E_keV_1);
                    }
                }
            }
        }
        cout << "File: " << i+1 << endl;
    }

    // background retrieving
    TFile* infile_b = new TFile("DATA/Day1/file_notturno_fondo.root");
    TTree* intree_b = (TTree*)infile_b->Get("acq_tree_0");
    TBranch* inbranch_b_0 = intree_b->GetBranch("acq_ch0");
    TBranch* inbranch_b_1 = intree_b->GetBranch("acq_ch1");
    inbranch_b_0->SetAddress(&indata_b_0.timetag);
    inbranch_b_1->SetAddress(&indata_b_1.timetag);

    TH2D* h_b = new TH2D("h_b","bkg",numBinsX,minX,maxX,numBinsY,minY,maxY);
    for (int i=0; i<inbranch_b_0->GetEntries(); ++i) {
        inbranch_b_0->GetEntry(i);
        inbranch_b_1->GetEntry(i);

        E_keV_0 = indata_b_0.qlong*m_ch_0 + q_ch_0 + gRandom->Rndm();
        E_keV_1 = indata_b_1.qlong*m_ch_1 + q_ch_1 + gRandom->Rndm();

        if ((E_keV_0<maxX) && (E_keV_0>minX)) {
            if ((E_keV_1<maxY) && (E_keV_1>minY)) {
                h_b->Fill(E_keV_0,E_keV_1);
            }
        }
    }

    // background subtraction
    float fraction = 110.0 / (24.0 * 60.0);
    h->Add(h_b, -fraction);
    

    if (draw==1) h->Draw();

    float c_11 = 9.37196e+04; // 1173 keV - ch_0
    float m_11 = 1.19095e+03;
    float s_11 = 2.79452e+01;
    float c_21 = 7.90883e+04; // 1333 keV - ch_0
    float m_21 = 1.35148e+03;
    float s_21 = 3.01299e+01;

    float c_12 = 4.05241e+04; // 1173 keV - ch_1
    float m_12 = 1.18916e+03;
    float s_12 = 2.83518e+01;
    float c_22 = 3.38967e+04; // 1333 keV - ch_1
    float m_22 = 1.34771e+03;
    float s_22 = 3.02504e+01;

    int binX_min_0 = h->GetXaxis()->FindBin(m_11-3.0*s_11);//1130);
    int binX_max_0 = h->GetXaxis()->FindBin(m_11+3.0*s_11);//1250);
    int binY_min_0 = h->GetYaxis()->FindBin(m_22-3.0*s_22);//1280);
    int binY_max_0 = h->GetYaxis()->FindBin(m_22+3.0*s_22);//1400);

    int binX_min_1 = h->GetXaxis()->FindBin(m_21-3.0*s_21);//1290);
    int binX_max_1 = h->GetXaxis()->FindBin(m_21+3.0*s_21);//1400);
    int binY_min_1 = h->GetYaxis()->FindBin(m_12-3.0*s_12);//1130);
    int binY_max_1 = h->GetYaxis()->FindBin(m_12+3.0*s_12);//1240);

    int C_1 = h->Integral(binX_min_0, binX_max_0, binY_min_0, binY_max_0);
    int C_2 = h->Integral(binX_min_1, binX_max_1, binY_min_1, binY_max_1);

    int sum_peak_counts = 4000;
    int N_11 = 6547168;//h->Integral(binX_min_0, binX_max_0, 1, numBinsY-1) + sum_peak_counts*11;
    int N_21 = 5956975;//h->Integral(binX_min_1, binX_max_1, 1, numBinsY-1) + sum_peak_counts*11;
    int N_22 = 2563329;//h->Integral(1, numBinsY-1, binY_min_0, binY_max_0) + sum_peak_counts*11;
    int N_12 = 2872168;//h->Integral(1, numBinsY-1, binY_min_1, binY_max_1) + sum_peak_counts*11;

    cout << "Counts for gamma_1 in D1 and gamma_2 in D2: " << C_1 << endl;
    cout << "Counts for gamma_2 in D1 and gamma_1 in D2: " << C_2 << endl;

    float d_1 = 12.0;
    float d_2 = 22.0;
    float r   = 3.75;

    float S_1 = 0.25 * pow((r/d_1),2.0);
    float S_2 = 0.25 * pow((r/d_2),2.0);

    float Abs_E_11 = (C_1*1.0) / (N_22*1.0);
    float Abs_E_12 = (C_2*1.0) / (N_21*1.0);
    float Abs_E_21 = (C_2*1.0) / (N_12*1.0);
    float Abs_E_22 = (C_1*1.0) / (N_11*1.0);

    float err_Abs_E_11 = Abs_E_11 * sqrt( (1.0/(C_1*1.0)) + (1.0/(N_22*1.0)) );
    float err_Abs_E_12 = Abs_E_12 * sqrt( (1.0/(C_2*1.0)) + (1.0/(N_21*1.0)) );
    float err_Abs_E_21 = Abs_E_21 * sqrt( (1.0/(C_2*1.0)) + (1.0/(N_12*1.0)) );
    float err_Abs_E_22 = Abs_E_22 * sqrt( (1.0/(C_1*1.0)) + (1.0/(N_11*1.0)) );

    float Int_E_11 = Abs_E_11 * (1.0/S_1);
    float Int_E_12 = Abs_E_12 * (1.0/S_2);
    float Int_E_21 = Abs_E_21 * (1.0/S_1);
    float Int_E_22 = Abs_E_22 * (1.0/S_2);

    float err_Int_E_11 = err_Abs_E_11 * (1.0/S_1);
    float err_Int_E_12 = err_Abs_E_12 * (1.0/S_2);
    float err_Int_E_21 = err_Abs_E_21 * (1.0/S_1);
    float err_Int_E_22 = err_Abs_E_22 * (1.0/S_2);

    ofstream ofile("Analysis/Day2/Day2_coinc_method.txt");
    ofile << "Day2 results:" << endl << endl;
    ofile << "Legend: " << endl;
    ofile << "Quantity_{index1,index2}" << endl;
    ofile << "index1 = energy of gamma (1=1173, 2=1333)" << endl;
    ofile << "index2 = index of detector (1=D1, 2=D2)" << endl;
    ofile << "d_index = distance from detector #index" << endl;
    ofile << "r = radius of detector surface" << endl;
    ofile << "S_index = geometric factor of detector #index" << endl;
    ofile << endl;

    ofile << "d_1 [cm]: " << d_1 << endl;
    ofile << "d_2 [cm]: " << d_2 << endl;
    ofile << "r   [cm]: " << r   << endl;
    ofile << endl;

    ofile << "S_1: " << S_1 << endl;
    ofile << "S_2: " << S_2 << endl;
    ofile << endl;

    ofile << "Counts N_11:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_11 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_11) << '\t' << '$' << endl;
    ofile << "Counts N_12:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_12 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_12) << '\t' << '$' << endl;
    ofile << "Counts N_21:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_21 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_21) << '\t' << '$' << endl;
    ofile << "Counts N_22:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_22 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_22) << '\t' << '$' << endl;
    ofile << endl;

    ofile << "Coincidences C_1:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << C_1 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(C_1) << '\t' << '$' << endl;
    ofile << "Coincidences C_2:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << C_2 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(C_2) << '\t' << '$' << endl;
    ofile << endl;

    ofile << "COINCIDENCE METHOD" << endl;
    ofile << "*********************************************************************************" << endl;
    ofile << "Absolute Efficiency Abs_E_11:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_11 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_11 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_12:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_12 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_12 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_21:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_21 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_21 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_22:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_22 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_22 << '\t' << '$' << endl;
    ofile << endl;

    ofile << "Intrinsic Efficiency Int_E_11:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_11 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_11 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_12:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_12 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_12 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_21:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_21 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_21 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_22:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_22 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_22 << '\t' << '$' << endl;
    ofile << "*********************************************************************************" << endl;
    ofile << endl;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void sumMethod(int numBinsX, int numBinsY, double minX, double maxX, double minY, double maxY, int draw=1) {
    // variables
    const int n_files = 12;
    float E_keV_0;
    float E_keV_1;
	slimport_data_t indata_0[n_files];
    slimport_data_t indata_1[n_files];
    slimport_data_t indata_b_0;
    slimport_data_t indata_b_1;

    // data retrieving
    TFile** infile = new TFile*[n_files];
    TTree** intree = new TTree*[n_files];
    TBranch** inbranch_0 = new TBranch*[n_files];
    TBranch** inbranch_1 = new TBranch*[n_files];

    for (int i=0; i<n_files; ++i) {
        intree[i] = new TTree;
        inbranch_0[i] = new TBranch;
        inbranch_1[i] = new TBranch;

        infile[i] = new TFile(Form("DATA/Day2/run_%d_logic_or.root",i+1));
        intree[i] = (TTree*)infile[i]->Get("acq_tree_0");
        inbranch_0[i] = intree[i]->GetBranch("acq_ch0");
        inbranch_1[i] = intree[i]->GetBranch("acq_ch1");
        inbranch_0[i]->SetAddress(&indata_0[i].timetag);
        inbranch_1[i]->SetAddress(&indata_1[i].timetag);
    }

	TH1F* h_0 = new TH1F("h_0","ch_0",numBinsX,minX,maxX);
    TH1F* h_1 = new TH1F("h_1","ch_1",numBinsY,minY,maxY);
	// histogram filling
    for (int i=0; i<n_files; ++i) {
        if (i!=1) {
            for (int j=0; j<inbranch_0[i]->GetEntries(); j++) {
                inbranch_0[i]->GetEntry(j);
                inbranch_1[i]->GetEntry(j);
                
                E_keV_0 = indata_0[i].qlong*m_ch_0 + q_ch_0 + gRandom->Rndm();
                E_keV_1 = indata_1[i].qlong*m_ch_1 + q_ch_1 + gRandom->Rndm();

                h_0->Fill(E_keV_0);
                h_1->Fill(E_keV_1);
            }
        }
        cout << "File: " << i+1 << endl;
    }

    // background retrieving
    TFile* infile_b = new TFile("DATA/Day1/file_notturno_fondo.root");
    TTree* intree_b = (TTree*)infile_b->Get("acq_tree_0");
    TBranch* inbranch_b_0 = intree_b->GetBranch("acq_ch0");
    TBranch* inbranch_b_1 = intree_b->GetBranch("acq_ch1");
    inbranch_b_0->SetAddress(&indata_b_0.timetag);
    inbranch_b_1->SetAddress(&indata_b_1.timetag);

    TH1F* h_b_0 = new TH1F("h_b_0","bkg ch_0",numBinsX,minX,maxX);
    TH1F* h_b_1 = new TH1F("h_b_1","bkg ch_1",numBinsY,minY,maxY);
    for (int i=0; i<inbranch_b_0->GetEntries(); ++i) {
        inbranch_b_0->GetEntry(i);
        inbranch_b_1->GetEntry(i);

        E_keV_0 = indata_b_0.qlong*m_ch_0 + q_ch_0 + gRandom->Rndm();
        E_keV_1 = indata_b_1.qlong*m_ch_1 + q_ch_1 + gRandom->Rndm();

        h_b_0->Fill(E_keV_0);
        h_b_1->Fill(E_keV_1);
    }

    // background subtraction
    float fraction = 110.0 / (24.0 * 60.0);
    h_0->Add(h_b_0, -fraction);
    h_1->Add(h_b_1, -fraction);
    

    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "ch_0", 800, 600);
        c1->cd();
        h_0->Draw();
        TCanvas* c2 = new TCanvas("c2", "ch_1", 800, 600);
        c2->cd();
        h_1->Draw();
    }

    ///*
    TF1* g_11 = new TF1("g_11","gaus",1000,1400);
    TF1* g_21 = new TF1("g_21","gaus",1150,1550);
    TF1* b_1  = new TF1("b_1" ,"pol3",1075,1500);
    TF1* p_1  = new TF1("p_1" ,"gaus",2300,2900);
    TF1* b_p1 = new TF1("b_p1","pol3",2300,2900);

    g_11->SetParameter(0, 9.37196e+04); // 1173 keV - ch_0
    g_11->SetParameter(1, 1.19095e+03);
    g_11->SetParameter(2, 2.79452e+01);
    g_21->SetParameter(0, 7.90883e+04); // 1333 keV - ch_0
    g_21->SetParameter(1, 1.35148e+03);
    g_21->SetParameter(2, 3.01299e+01);
    b_1 ->SetParameter(0, 1.36646e+06); // compton 1173, 1333 keV - ch_0
    b_1 ->SetParameter(1,-2.96508e+03);
    b_1 ->SetParameter(2, 2.14790e+00);
    b_1 ->SetParameter(3,-5.18865e-04);
    p_1 ->SetParameter(0, 3.11276e+02); // 2506 keV - ch_0
    p_1 ->SetParameter(1, 2.60196e+03);
    p_1 ->SetParameter(2, 5.12564e+01);
    b_p1->SetParameter(0, 8.02600e+04); // compton 2506 keV - ch_0
    b_p1->SetParameter(1,-9.00392e+01);
    b_p1->SetParameter(2, 3.37028e-02);
    b_p1->SetParameter(3,-4.20761e-06);

    float c_11 = 9.37196e+04; // 1173 keV - ch_0
    float m_11 = 1.19095e+03;
    float s_11 = 2.79452e+01;
    float c_21 = 7.90883e+04; // 1333 keV - ch_0
    float m_21 = 1.35148e+03;
    float s_21 = 3.01299e+01;
    float c_p1 = 3.11276e+02; // 2506 keV - ch_0
    float m_p1 = 2.60196e+03;
    float s_p1 = 5.12564e+01;

    TF1* g_12 = new TF1("g_12","gaus",1000,1400);
    TF1* g_22 = new TF1("g_22","gaus",1150,1550);
    TF1* b_2  = new TF1("b_2" ,"pol3",1075,1500);
    TF1* p_2  = new TF1("p_2" ,"gaus",2300,2900);
    TF1* b_p2 = new TF1("b_p2","pol3",2300,2900);

    g_12->SetParameter(0, 4.05241e+04); // 1173 keV - ch_1
    g_12->SetParameter(1, 1.18916e+03);
    g_12->SetParameter(2, 2.83518e+01);
    g_22->SetParameter(0, 3.38967e+04); // 1333 keV - ch_1
    g_22->SetParameter(1, 1.34771e+03);
    g_22->SetParameter(2, 3.02504e+01);
    b_2 ->SetParameter(0, 7.31748e+05); // compton 1173, 1333 keV - ch_1
    b_2 ->SetParameter(1,-1.62365e+03);
    b_2 ->SetParameter(2, 1.20250e+00);
    b_2 ->SetParameter(3,-2.96826e-04);
    p_2 ->SetParameter(0, 8.57618e+01); // 2506 keV - ch_1
    p_2 ->SetParameter(1, 2.56782e+03);
    p_2 ->SetParameter(2, 6.05132e+01);
    b_p2->SetParameter(0, 6.70838e+03); // compton 2506 keV - ch_1
    b_p2->SetParameter(1,-6.95624e+00);
    b_p2->SetParameter(2, 2.40342e-03);
    b_p2->SetParameter(3,-2.76674e-07);

    float c_12 = 4.05241e+04; // 1173 keV - ch_1
    float m_12 = 1.18916e+03;
    float s_12 = 2.83518e+01;
    float c_22 = 3.38967e+04; // 1333 keV - ch_1
    float m_22 = 1.34771e+03;
    float s_22 = 3.02504e+01;
    float c_p2 = 8.57618e+01; // 2506 keV - ch_1
    float m_p2 = 2.56782e+03;
    float s_p2 = 6.05132e+01;

    float N_11 = g_11->Integral(m_11 - 3.0*s_11, m_11 + 3.0*s_11);
    float N_21 = g_21->Integral(m_21 - 3.0*s_21, m_21 + 3.0*s_21);
    float N_22 = g_22->Integral(m_22 - 3.0*s_22, m_22 + 3.0*s_22);
    float N_12 = g_12->Integral(m_12 - 3.0*s_12, m_12 + 3.0*s_12);
    float P_1  = p_1 ->Integral(m_p1 - 3.0*s_p1, m_p1 + 3.0*s_p1);
    float P_2  = p_2 ->Integral(m_p2 - 3.0*s_p2, m_p2 + 3.0*s_p2);
    //*/
    /*
    int binX_min_2506_0 = h_0->GetXaxis()->FindBin(2500);
    int binX_max_2506_0 = h_0->GetXaxis()->FindBin(2800);
    int binX_min_2506_1 = h_1->GetXaxis()->FindBin(2430);
    int binX_max_2506_1 = h_1->GetXaxis()->FindBin(2800);

    int binX_min_1173_0 = h_0->GetXaxis()->FindBin(1120);
    int binX_max_1173_0 = h_0->GetXaxis()->FindBin(1270);
    int binX_min_1173_1 = h_1->GetXaxis()->FindBin(1120);
    int binX_max_1173_1 = h_1->GetXaxis()->FindBin(1270);

    int binX_min_1333_0 = h_0->GetXaxis()->FindBin(1270);
    int binX_max_1333_0 = h_0->GetXaxis()->FindBin(1450);
    int binX_min_1333_1 = h_1->GetXaxis()->FindBin(1270);
    int binX_max_1333_1 = h_1->GetXaxis()->FindBin(1450);


    int P_1 = h_0->Integral(binX_min_2506_0, binX_max_2506_0);
    int P_2 = h_1->Integral(binX_min_2506_1, binX_max_2506_1);

    int N_11 = h_0->Integral(binX_min_1173_0, binX_max_1173_0) - 1226032;
    int N_21 = h_0->Integral(binX_min_1333_0, binX_max_1333_0) - 392272;
    int N_22 = h_1->Integral(binX_min_1333_1, binX_max_1333_1) - 188079;
    int N_12 = h_1->Integral(binX_min_1173_1, binX_max_1173_1) - 610043;
    */

    cout << "Counts for gamma_1+gamma_2 in D1: " << P_1 << endl;
    cout << "Counts for gamma_1+gamma_2 in D2: " << P_2 << endl;

    float d_1 = 12.0;
    float d_2 = 22.0;
    float r   = 3.75;

    float S_1 = 0.25 * pow((r/d_1),2.0);
    float S_2 = 0.25 * pow((r/d_2),2.0);

    float Abs_E_11 = (P_1*1.0) / (N_21*1.0);
    float Abs_E_12 = (P_2*1.0) / (N_22*1.0);
    float Abs_E_21 = (P_1*1.0) / (N_11*1.0);
    float Abs_E_22 = (P_2*1.0) / (N_12*1.0);

    float err_Abs_E_11 = Abs_E_11 * sqrt( (1.0/(P_1*1.0)) + (1.0/(N_21*1.0)) );
    float err_Abs_E_12 = Abs_E_12 * sqrt( (1.0/(P_2*1.0)) + (1.0/(N_22*1.0)) );
    float err_Abs_E_21 = Abs_E_21 * sqrt( (1.0/(P_1*1.0)) + (1.0/(N_11*1.0)) );
    float err_Abs_E_22 = Abs_E_22 * sqrt( (1.0/(P_2*1.0)) + (1.0/(N_12*1.0)) );

    float Int_E_11 = Abs_E_11 * (1.0/S_1);
    float Int_E_12 = Abs_E_12 * (1.0/S_2);
    float Int_E_21 = Abs_E_21 * (1.0/S_1);
    float Int_E_22 = Abs_E_22 * (1.0/S_2);

    float err_Int_E_11 = err_Abs_E_11 * (1.0/S_1);
    float err_Int_E_12 = err_Abs_E_12 * (1.0/S_2);
    float err_Int_E_21 = err_Abs_E_21 * (1.0/S_1);
    float err_Int_E_22 = err_Abs_E_22 * (1.0/S_2);

    ofstream ofile("Analysis/Day2/Day2_sum_method.txt");
    ofile << "Day2 results:" << endl << endl;
    ofile << "Legend: " << endl;
    ofile << "Quantity_{index1,index2}" << endl;
    ofile << "index1 = energy of gamma (1=1173, 2=1333)" << endl;
    ofile << "index2 = index of detector (1=D1, 2=D2)" << endl;
    ofile << "d_index = distance from detector #index" << endl;
    ofile << "r = radius of detector surface" << endl;
    ofile << "S_index = geometric factor of detector #index" << endl;
    ofile << endl;

    ofile << "d_1 [cm]: " << d_1 << endl;
    ofile << "d_2 [cm]: " << d_2 << endl;
    ofile << "r   [cm]: " << r   << endl;
    ofile << endl;

    ofile << "S_1: " << S_1 << endl;
    ofile << "S_2: " << S_2 << endl;
    ofile << endl;

    ofile << "Counts N_11:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_11 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_11) << '\t' << '$' << endl;
    ofile << "Counts N_12:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_12 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_12) << '\t' << '$' << endl;
    ofile << "Counts N_21:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_21 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_21) << '\t' << '$' << endl;
    ofile << "Counts N_22:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << N_22 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(N_22) << '\t' << '$' << endl;
    ofile << endl;

    ofile << "Sum peak counts P_1:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << P_1 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(P_1) << '\t' << '$' << endl;
    ofile << "Sum peak counts P_2:" << '\t' << '$' << '\t' << fixed << setw(10) << setprecision(0) << P_2 << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(0) << sqrt(P_2) << '\t' << '$' << endl;
    ofile << endl;

    ofile << "SUM PEAK METHOD" << endl;
    ofile << "*********************************************************************************" << endl;
    ofile << "Absolute Efficiency Abs_E_11:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_11 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_11 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_12:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_12 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_12 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_21:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_21 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_21 << '\t' << '$' << endl;
    ofile << "Absolute Efficiency Abs_E_22:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Abs_E_22 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Abs_E_22 << '\t' << '$' << endl;
    ofile << endl;

    ofile << "Intrinsic Efficiency Int_E_11:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_11 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_11 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_12:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_12 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_12 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_21:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_21 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_21 << '\t' << '$' << endl;
    ofile << "Intrinsic Efficiency Int_E_22:" << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(6) << Int_E_22 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(6) << err_Int_E_22 << '\t' << '$' << endl;
    ofile << "*********************************************************************************" << endl;
    ofile << endl;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void geom( int numBins, double minX, double maxX, int draw=1) {
	ofstream out("Analysis/Day2/NewGeom.txt");
	// variables
	slimport_data_t indata;
    const char* name_file[6] = {
        "DATA/Day2/0_degree_logic_or.root",
        "DATA/Day2/20_degree_logic_or.root",
        "DATA/Day2/40_degree_logic_or.root",
        "DATA/Day2/50_degree_logic_or.root",
        "DATA/Day2/70_degree_logic_or.root",
        "DATA/Day2/90_degree_logic_or.root"
    };
	int Ang[6] ={0,20,40,50,70,90};
	float ext[4] = {4000.0,4800.0,5500.0,6300.0};
	float m, q, zero, fondo;
	long start, finish;

	//out << "Gamma1-fondo \t ErrGamma1 \t Gamma2 \t ErrGamma2\tGamma1+2\tErr1+2" << std::endl;
	out << "1/Gamma1" << std::endl;
	
	TF1 *fa1 = new TF1("fa1","gaus(0) + gaus(3) + pol3(6)",4000,6600);
	fa1->SetParameter(0 , 4.53306e+03);
	fa1->SetParameter(1 , 5.12719e+03);
	fa1->SetParameter(2 , 1.15037e+02);
	fa1->SetParameter(3 , 3.91710e+03);
	fa1->SetParameter(4 , 5.81651e+03);
	fa1->SetParameter(5 , 1.24694e+02);
	fa1->SetParameter(6 , 1.56872e+04);
	fa1->SetParameter(7 ,-6.72006e+00);
	fa1->SetParameter(8 , 9.78878e-04);
	fa1->SetParameter(9 ,-4.84504e-08);
	
	TCanvas** c = new TCanvas*[6];
	for (int j =0; j<6; j++){
		if (draw==1) {
			c[j] = new TCanvas(Form("c%d",j+1),Form("c%d",j+1), 800,600);
			c[j]->cd();
		}
			

		TFile *infile = new TFile(name_file[j]);
		TTree *intree = (TTree*)infile->Get("acq_tree_0");
		TBranch *inbranch = intree->GetBranch("acq_ch1");
		inbranch->SetAddress(&indata.timetag);
		
		TH1F *h_spectrum = new TH1F("h_spectrum","Total spectrum",numBins,minX,maxX);

		inbranch->GetEntry(1); 
		start = indata.timetag; 
		finish = start + 40000000000;
		// histogram filling
		for (int i=0; i<inbranch->GetEntries(); i++) {			
			inbranch->GetEntry(i);
			if(indata.timetag > start && indata.timetag < finish)
				h_spectrum->Fill(indata.qlong);
		}

		//Fit Lineare sul Fondo
		float  width = h_spectrum->GetXaxis()->GetBinWidth(10) ; 
		h_spectrum->Fit(fa1, "", "", 4000, 6600);

		TF1* g_12 = new TF1("g_12","gaus",4000,6600);
		TF1* g_22 = new TF1("g_22","gaus",4000,6600);
		
		float c_12 = fa1->GetParameter(0); // 1173 keV - ch_1
		float m_12 = fa1->GetParameter(1);
		float s_12 = fa1->GetParameter(2);
		float c_22 = fa1->GetParameter(3); // 1333 keV - ch_1
		float m_22 = fa1->GetParameter(4);
		float s_22 = fa1->GetParameter(5);

		g_12->SetParameter(0, c_12);
		g_12->SetParameter(1, m_12);
		g_12->SetParameter(2, s_12);
		g_22->SetParameter(0, c_22);
		g_22->SetParameter(1, m_22);
		g_22->SetParameter(2, s_22);


		//Integral1
    	float Gamma1 = g_12->Integral(m_12 - 3.0*s_12, m_12 + 3.0*s_12) / width;
		float Gamma2 = g_22->Integral(m_22 - 3.0*s_22, m_22 + 3.0*s_22) / width;
		
		//out << Gamma1-fondo << " \t " << sqrt(Gamma1-fondo) <<" \t " << Gamma2 << " \t " << sqrt(Gamma2) << " \t " << Gamma1+Gamma2-fondo << " \t " << sqrt(Gamma1-fondo + Gamma2 ) <<  std::endl;

		out << static_cast<float>(Gamma2 ) << std::endl;

		}	
}
// ********************************************************************************************
// ********************************************************************************************