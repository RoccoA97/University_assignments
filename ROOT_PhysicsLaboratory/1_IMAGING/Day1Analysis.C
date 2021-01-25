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



void plotSpectrum(int draw=1) {

    int nbin[8] = {
        200,
        200,
        200,
        200,
        200,
        200,
        200,
        200
    };
    int range[8][2] = {
        { 0,    14000},
        { 0,    10000},
        { 0,    10000},
        { 0,    10000},
        { 0,    10000},
        { 0,    10000},
        { 0,    10000},
        { 0,    10000}
    };

    float data[8][8] = {
        {4530,	10,     860,    40,     12160,     40,     1670,    60},
        {2710,	 4,     598,     6,      7010,     30,     1280,    70},
        {3091,	 3,     551,     3,      7870,     20,      980,    20},
        {1619,	 1,     244,     2,      4143,      7,      496,     9},
        {2311,	 2,     379,     3,      5940,     10,      620,    20},
        {1951,	 2,     291,     2,      4986,      8,      440,    10},
        {1712,	 3,     488,     5,      4120,      8,     1100,   100},
        {2654,	 4,     564,     7,      6810,     10,      710,    20}
    };
    float parameters[8][4] = {
        {57,	    1,	    0.1001,     0.0002},
        {29,	    2,	    0.1777,     0.0003},
        {17,	    2,	    0.1599,     0.0003},
        {21,	    2,	    0.3027,     0.0006},
        {24,	    2,	    0.2105,     0.0004},
        {20,	    2,	    0.2517,     0.0005},
        {-32,       2,	    0.3173,     0.0006},
        {23,	    2,	    0.1838,     0.0003}
    };

    float mean_peak_1[8];
    float mean_peak_2[8];
    float sigma_peak_1[8];
    float sigma_peak_2[8];
    float err_mean_peak_1[8];
    float err_mean_peak_2[8];
    float err_sigma_peak_1[8];
    float err_sigma_peak_2[8];
    float fwhm_1[8];
    float fwhm_2[8];
    float err_fwhm_1[8];
    float err_fwhm_2[8];
    float resolution_1[8];
    float resolution_2[8];
    float err_resolution_1[8];
    float err_resolution_2[8];
    float a[8];
    float b[8];
    float err_a[8];
    float err_b[8];

    for (int i=0; i<8; ++i) {
        mean_peak_1[i]      = data[i][0];
        mean_peak_2[i]      = data[i][4];
        sigma_peak_1[i]     = data[i][2];
        sigma_peak_2[i]     = data[i][6];
        err_mean_peak_1[i]  = data[i][1];
        err_mean_peak_2[i]  = data[i][5];
        err_sigma_peak_1[i] = data[i][3];
        err_sigma_peak_2[i] = data[i][7];

        a[i] = parameters[i][0];
        b[i] = parameters[i][2];
        err_a[i] = parameters[i][1];
        err_b[i] = parameters[i][3];
    }

    float y_peaks[2] = {
        511.0,
        1274.537
    };
    float err_y_peaks[2] = {
        0.0,
        0.007
    };
    float x[8][2] = {
        {mean_peak_1[0],mean_peak_2[0]},
        {mean_peak_1[1],mean_peak_2[1]},
        {mean_peak_1[2],mean_peak_2[2]},
        {mean_peak_1[3],mean_peak_2[3]},
        {mean_peak_1[4],mean_peak_2[4]},
        {mean_peak_1[5],mean_peak_2[5]},
        {mean_peak_1[6],mean_peak_2[6]},
        {mean_peak_1[7],mean_peak_2[7]}
    };
    float err_x[8][2] = {
        {err_mean_peak_1[0],err_mean_peak_2[0]},
        {err_mean_peak_1[1],err_mean_peak_2[1]},
        {err_mean_peak_1[2],err_mean_peak_2[2]},
        {err_mean_peak_1[3],err_mean_peak_2[3]},
        {err_mean_peak_1[4],err_mean_peak_2[4]},
        {err_mean_peak_1[5],err_mean_peak_2[5]},
        {err_mean_peak_1[6],err_mean_peak_2[6]},
        {err_mean_peak_1[7],err_mean_peak_2[7]}
    };


    slimport_data_t indata_1;
    slimport_data_t indata_2;
    slimport_data_t indata_3;
    slimport_data_t indata_4;
    slimport_data_t indata_5;
    slimport_data_t indata_6;
    slimport_data_t indata_7;
    slimport_data_t indata_8;

	TFile* infile_1 = new TFile("DATA/Day1/D1_calibration.root");
    TFile* infile_2 = new TFile("DATA/Day1/D2_calibration.root");
    TFile* infile_3 = new TFile("DATA/Day1/D3_calibration.root");
    TFile* infile_4 = new TFile("DATA/Day1/D4_calibration.root");
    TFile* infile_5 = new TFile("DATA/Day1/D5_calibration.root");
    TFile* infile_6 = new TFile("DATA/Day1/D6_calibration.root");
    TFile* infile_7 = new TFile("DATA/Day1/D7_calibration.root");
    TFile* infile_8 = new TFile("DATA/Day1/D8_calibration.root");

	TTree* intree_1 = (TTree*)infile_1->Get("acq_tree_0");
    TTree* intree_2 = (TTree*)infile_2->Get("acq_tree_0");
    TTree* intree_3 = (TTree*)infile_3->Get("acq_tree_0");
    TTree* intree_4 = (TTree*)infile_4->Get("acq_tree_0");
    TTree* intree_5 = (TTree*)infile_5->Get("acq_tree_1");
    TTree* intree_6 = (TTree*)infile_6->Get("acq_tree_1");
    TTree* intree_7 = (TTree*)infile_7->Get("acq_tree_1");
    TTree* intree_8 = (TTree*)infile_8->Get("acq_tree_1");

    TBranch* inbranch_1 = intree_1->GetBranch("acq_ch0");
    TBranch* inbranch_2 = intree_2->GetBranch("acq_ch1");
    TBranch* inbranch_3 = intree_3->GetBranch("acq_ch2");
    TBranch* inbranch_4 = intree_4->GetBranch("acq_ch3");
    TBranch* inbranch_5 = intree_5->GetBranch("acq_ch0");
    TBranch* inbranch_6 = intree_6->GetBranch("acq_ch1");
    TBranch* inbranch_7 = intree_7->GetBranch("acq_ch2");
    TBranch* inbranch_8 = intree_8->GetBranch("acq_ch3");

    inbranch_1->SetAddress(&indata_1.timetag);
    inbranch_2->SetAddress(&indata_2.timetag);
    inbranch_3->SetAddress(&indata_3.timetag);
    inbranch_4->SetAddress(&indata_4.timetag);
    inbranch_5->SetAddress(&indata_5.timetag);
    inbranch_6->SetAddress(&indata_6.timetag);
    inbranch_7->SetAddress(&indata_7.timetag);
    inbranch_8->SetAddress(&indata_8.timetag);

    TH1F** h = new TH1F*[8];
    for (int i=0; i<8; ++i) {
        h[i] = new TH1F(Form("Detector D%d",i+1), Form("Spectrum acquired by detector D%d",i+1), nbin[i], range[i][0], range[i][1]);
    }

    for (int i=0; i<inbranch_1->GetEntries(); i++) {  inbranch_1->GetEntry(i); h[0]->Fill(indata_1.qlong);    }
    for (int i=0; i<inbranch_2->GetEntries(); i++) {  inbranch_2->GetEntry(i); h[1]->Fill(indata_2.qlong);    }
    for (int i=0; i<inbranch_3->GetEntries(); i++) {  inbranch_3->GetEntry(i); h[2]->Fill(indata_3.qlong);    }
    for (int i=0; i<inbranch_4->GetEntries(); i++) {  inbranch_4->GetEntry(i); h[3]->Fill(indata_4.qlong);    }
    for (int i=0; i<inbranch_5->GetEntries(); i++) {  inbranch_5->GetEntry(i); h[4]->Fill(indata_5.qlong);    }
    for (int i=0; i<inbranch_6->GetEntries(); i++) {  inbranch_6->GetEntry(i); h[5]->Fill(indata_6.qlong);    }
    for (int i=0; i<inbranch_7->GetEntries(); i++) {  inbranch_7->GetEntry(i); h[6]->Fill(indata_7.qlong);    }
    for (int i=0; i<inbranch_8->GetEntries(); i++) {  inbranch_8->GetEntry(i); h[7]->Fill(indata_8.qlong);    }

    TF1** f = new TF1*[8];
    TF1* f[0] = new TF1("Calibration D1","[0] + [1]*x", range[0][0], range[0][1]);
    TF1* f[1] = new TF1("Calibration D2","[0] + [1]*x", range[1][0], range[1][1]);
    TF1* f[2] = new TF1("Calibration D3","[0] + [1]*x", range[2][0], range[2][1]);
    TF1* f[3] = new TF1("Calibration D4","[0] + [1]*x", range[3][0], range[3][1]);
    TF1* f[4] = new TF1("Calibration D5","[0] + [1]*x", range[4][0], range[4][1]);
    TF1* f[5] = new TF1("Calibration D6","[0] + [1]*x", range[5][0], range[5][1]);
    TF1* f[6] = new TF1("Calibration D7","[0] + [1]*x", range[6][0], range[6][1]);
    TF1* f[7] = new TF1("Calibration D8","[0] + [1]*x", range[7][0], range[7][1]);

    f[0]->SetParameters(a[0],b[0]);
    f[1]->SetParameters(a[1],b[1]);
    f[2]->SetParameters(a[2],b[2]);
    f[3]->SetParameters(a[3],b[3]);
    f[4]->SetParameters(a[4],b[4]);
    f[5]->SetParameters(a[5],b[5]);
    f[6]->SetParameters(a[6],b[6]);
    f[7]->SetParameters(a[7],b[7]);

    TGraphErrors** gr = new TGraphErrors*[8];
    TGraphErrors* gr[0] = new TGraphErrors(2,x[0],y_peaks,err_x[0],err_y_peaks);
    TGraphErrors* gr[1] = new TGraphErrors(2,x[1],y_peaks,err_x[1],err_y_peaks);
    TGraphErrors* gr[2] = new TGraphErrors(2,x[2],y_peaks,err_x[2],err_y_peaks);
    TGraphErrors* gr[3] = new TGraphErrors(2,x[3],y_peaks,err_x[3],err_y_peaks);
    TGraphErrors* gr[4] = new TGraphErrors(2,x[4],y_peaks,err_x[4],err_y_peaks);
    TGraphErrors* gr[5] = new TGraphErrors(2,x[5],y_peaks,err_x[5],err_y_peaks);
    TGraphErrors* gr[6] = new TGraphErrors(2,x[6],y_peaks,err_x[6],err_y_peaks);
    TGraphErrors* gr[7] = new TGraphErrors(2,x[7],y_peaks,err_x[7],err_y_peaks);

    if (draw==1) {
        for (int i=0; i<8; ++i) {
            h[i]->Draw();
            h[i]->GetXaxis()->SetTitle("Energy [a.u.]");
            h[i]->GetYaxis()->SetTitle("Number of photons");
            gStyle->SetStatX(0.89);
            gStyle->SetStatY(0.88);
            gPad->Update();
            gPad->SaveAs(Form("Analysis/Day1/D%d_spectrum.pdf",i+1));
        }
    }
    gPad->Clear();
    if (draw==1) {
        for (int i=0; i<8; ++i) {
            gr[i]->Draw("AP");
            gr[i]->SetTitle(Form("Calibration of detector D%i",i+1));
            f[i]->Draw("SAME");
            f[i]->SetLineWidth(0.5);
            gr[i]->GetXaxis()->SetTitle("Energy [a.u.]");
            gr[i]->GetYaxis()->SetTitle("Energy [keV]");
            gStyle->SetStatX(0.89);
            gStyle->SetStatY(0.88);
            gPad->Update();
            gPad->SaveAs(Form("Analysis/Day1/D%d_calibration.pdf",i+1));
            gPad->Clear();
        }
    }

    for (int i=0; i<8; ++i) {
        fwhm_1[i] = 2.355*sigma_peak_1[i];
        fwhm_2[i] = 2.355*sigma_peak_2[i];
        err_fwhm_1[i] = 2.355*err_sigma_peak_1[i];
        err_fwhm_2[i] = 2.355*err_sigma_peak_2[i];
        resolution_1[i] = fwhm_1[i] / mean_peak_1[i];
        resolution_2[i] = fwhm_2[i] / mean_peak_2[i];
        err_resolution_1[i] = resolution_1[i] * sqrt( pow(err_fwhm_1[i]/fwhm_1[i],2.0) + pow(err_mean_peak_1[i]/mean_peak_1[i],2.0) );
        err_resolution_2[i] = resolution_2[i] * sqrt( pow(err_fwhm_2[i]/fwhm_2[i],2.0) + pow(err_mean_peak_2[i]/mean_peak_2[i],2.0) );
    }

    ofstream ofile("Analysis/Day1/Day1.txt");
    cout << "*************" << endl;
    ofile << "LaTeX: parameters of the gaussian fits \\ pm error:" << endl;
    for (int i=0; i<8; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << mean_peak_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_mean_peak_1[i] << '\t' << '$' << '\t'
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << fwhm_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_fwhm_1[i] << '\t' << '$' << '\t'
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << resolution_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_resolution_1[i] << '\t' << '$' << '\t'
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << mean_peak_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_mean_peak_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << fwhm_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_fwhm_2[i] << '\t' << '$' << '\t'
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(12) << setprecision(5) << resolution_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(12) << setprecision(5) << err_resolution_2[i] << '\t' << '$' << "\t\\\\" << endl;
    }
}