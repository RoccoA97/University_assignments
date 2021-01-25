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





// Gaussian + linear background fits - ch_0
float mu_peaks_ch_0[3] = {
    2.50990e+002,   // 59.5 keV
    4.46658e+003,   // 1173 keV
    5.06913e+003    // 1333 keV
};
float err_mu_peaks_ch_0[3] = {
    3.26502e-001,   
    3.91850e-001,
    3.63903e-001
};
float sigma_peaks_ch_0[3] = {
    3.72595e+001,   // 59.5 keV
    1.08362e+002,   // 1173 keV
    1.19165e+002    // 1333 keV
};
float err_sigma_peaks_ch_0[3] = {
    5.37590e-001,  
    5.22952e-001,
    4.20916e-001
};

// Gaussian + linear background fits - ch_1
float mu_peaks_ch_1[3] = {
    2.80161e+002,   // 59.5 keV
    5.05633e+003,   // 1173 keV
    5.73275e+003    // 1333 keV
};
float err_mu_peaks_ch_1[3] = {
    9.31970e-002,
    6.03308e-001,
    5.23680e-001
};
float sigma_peaks_ch_1[3] = {
    4.32534e+001,   // 59.5 keV
    1.20050e+002,   // 1173 keV
    1.32981e+002    // 1333 keV
};
float err_sigma_peaks_ch_1[3] = {
    1.28363e-001,
    7.51249e-001,
    6.62424e-001
};

// Energies
float y_peaks[3] = {
    59.5,
    1173.0,
    1333.0
};
float err_y_peaks[3] = {
    0.001,
    0.001,
    0.001
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
void calibrate(int draw=1) {
    TGraphErrors* ge_ch_0 = new TGraphErrors(3, mu_peaks_ch_0, y_peaks, err_mu_peaks_ch_0, err_y_peaks);
    TGraphErrors* ge_ch_1 = new TGraphErrors(3, mu_peaks_ch_1, y_peaks, err_mu_peaks_ch_1, err_y_peaks);
    TF1* f_ch_0 = new TF1("Calibration ch_0", "pol1", 0, 7000);
    TF1* f_ch_1 = new TF1("Calibration ch_1", "pol1", 0, 7000);

    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "Calibration ch_0", 800, 600);
        c1->cd();
        ge_ch_0->Draw("AP");
        ge_ch_0->Fit(f_ch_0,"","",0,7000);
        ge_ch_0->GetXaxis()->SetTitle("Energy [a.u.]");
        ge_ch_0->GetYaxis()->SetTitle("Energy [keV]");
        ge_ch_0->GetYaxis()->SetRangeUser(0,1500);
        ge_ch_0->SetMarkerStyle(21);
        ge_ch_0->SetTitle("");
        gPad->Update();
        gPad->SaveAs("Analysis/Day1/D1_calibration.pdf");

        TCanvas* c2 = new TCanvas("c2", "Calibration ch_1", 800, 600);
        c2->cd();
        ge_ch_1->Draw("AP");
        ge_ch_1->Fit(f_ch_1,"","",0,7000);
        ge_ch_1->GetXaxis()->SetTitle("Energy [a.u.]");
        ge_ch_1->GetYaxis()->SetTitle("Energy [keV]");
        ge_ch_1->GetYaxis()->SetRangeUser(0,1500);
        ge_ch_1->SetMarkerStyle(21);
        ge_ch_1->SetTitle("");
        gPad->Update();
        gPad->SaveAs("Analysis/Day1/D2_calibration.pdf");
    }
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
TH1F* getHistoAtTimeT(char* file_name, short chan, int numBins, double minX, double maxX, double t, int draw=1) {
    // variables
	slimport_data_t indata;
    double time;
	TFile *infile = new TFile(file_name);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1F *h = new TH1F("h","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
        time = indata.timetag * 4.0e-9;
        if (time < t) {
            h->Fill(indata.qlong);
        }
	}
    if (draw==1) h->Draw();

    return h;
}
// ********************************************************************************************
// ********************************************************************************************





// ********************************************************************************************
// ********************************************************************************************
void resolution() {
    const int n_times = 8;
    /*
    ch_0
    EXT PARAMETER                APPROXIMATE        STEP         FIRST
    NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
    1  p0           9.56262e+03   2.56174e+01   3.39864e-02  -1.93391e-07
    2  p1           4.46677e+03   2.64512e-01   1.06496e-03   1.64636e-05
    3  p2           1.09321e+02   2.53774e-01   2.96072e-04  -1.37859e-04
    4  p3           8.06113e+03   2.15159e+01   2.85858e-02  -1.16555e-06
    5  p4           5.06943e+03   2.79322e-01   1.20865e-03  -2.78707e-05
    6  p5           1.19638e+02   2.49420e-01   2.98918e-04  -2.87900e-04
    7  p6           6.78050e+04   2.07290e+01   1.61660e-02  -1.04069e-05
    8  p7          -3.64473e+01   5.73951e-03   8.68972e-06  -5.49367e-02
    9  p8           6.61700e-03   1.04799e-06   1.57762e-09  -2.88814e+02
    10  p9         -4.04794e-07   1.56162e-10   9.65104e-14  -1.51401e+06
    */

    /*
    EXT PARAMETER                APPROXIMATE        STEP         FIRST
    NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
    1  p0           3.76547e+03   1.51484e+01   2.25009e-02   3.47517e-06
    2  p1           5.05874e+03   4.39026e-01   1.20610e-03   9.30344e-05
    3  p2           1.21247e+02   4.10222e-01   2.73696e-07   6.11761e-04
    4  p3           3.18799e+03   1.29367e+01   1.91647e-02  -6.28645e-06
    5  p4           5.73361e+03   4.65333e-01   1.36700e-03  -1.60165e-05
    6  p5           1.31890e+02   4.10209e-01   5.64279e-04   9.91698e-05
    7  p6           3.08676e+04   1.05020e+01   7.35941e-03  -5.95073e-06
    8  p7          -1.51358e+01   2.53455e-03   3.60866e-06  -4.93116e-02
    9  p8           2.50429e-03   4.00558e-07   5.97070e-10  -3.64501e+02
    10  p9         -1.39315e-07   5.20554e-11   3.32152e-14  -2.53637e+06
    */
    TH1F** h_0 = new TH1F*[n_times];
    h_0[0] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100,  30);
    h_0[1] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100,  60);
    h_0[2] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100,  90);
    h_0[3] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100, 120);
    h_0[4] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100, 150);
    h_0[5] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100, 180);
    h_0[6] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100, 210);
    h_0[7] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 0, 500, 100, 10100, 240);

    TH1F** h_1 = new TH1F*[n_times];
    h_1[0] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100,  30);
    h_1[1] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100,  60);
    h_1[2] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100,  90);
    h_1[3] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100, 120);
    h_1[4] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100, 150);
    h_1[5] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100, 180);
    h_1[6] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100, 210);
    h_1[7] = getHistoAtTimeT("DATA/Day1/Co60_new_calibration_ch1_ch0.root", 1, 500, 100, 10100, 240);

    TF1** f_0 = new TF1*[n_times];
    f_0[0] = new TF1("f_030_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[1] = new TF1("f_060_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[2] = new TF1("f_090_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[3] = new TF1("f_120_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[4] = new TF1("f_150_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[5] = new TF1("f_180_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[6] = new TF1("f_210_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);
    f_0[7] = new TF1("f_240_0", "gaus(0)+gaus(3)+pol3(6)", 3600, 5600);

    TF1** f_1 = new TF1*[n_times];
    f_1[0] = new TF1("f_030_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[1] = new TF1("f_060_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[2] = new TF1("f_090_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[3] = new TF1("f_120_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[4] = new TF1("f_150_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[5] = new TF1("f_180_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[6] = new TF1("f_210_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);
    f_1[7] = new TF1("f_240_1", "gaus(0)+gaus(3)+pol3(6)", 4400, 6500);

    for (int i=0; i<n_times; ++i) {
        f_0[i]->SetParameter(0,  9.56262e+03);
        f_0[i]->SetParameter(1,  4.46677e+03);
        f_0[i]->SetParameter(2,  1.09321e+02);
        f_0[i]->SetParameter(3,  8.06113e+03);
        f_0[i]->SetParameter(4,  5.06943e+03);
        f_0[i]->SetParameter(5,  1.19638e+02);
        f_0[i]->SetParameter(6,  6.78050e+04);
        f_0[i]->SetParameter(7, -3.64473e+01);
        f_0[i]->SetParameter(8,  6.61700e-03);
        f_0[i]->SetParameter(9, -4.04794e-07);

        f_1[i]->SetParameter(0,  3.76547e+03);
        f_1[i]->SetParameter(1,  5.05874e+03);
        f_1[i]->SetParameter(2,  1.21247e+02);
        f_1[i]->SetParameter(3,  3.18799e+03);
        f_1[i]->SetParameter(4,  5.73361e+03);
        f_1[i]->SetParameter(5,  1.31890e+02);
        f_1[i]->SetParameter(6,  3.08676e+04);
        f_1[i]->SetParameter(7, -1.51358e+01);
        f_1[i]->SetParameter(8,  2.50429e-03);
        f_1[i]->SetParameter(9, -1.39315e-07);
    }

    TCanvas** cc_0 = new TCanvas*[n_times];
    TCanvas** cc_1 = new TCanvas*[n_times];
    for (int i=0; i<n_times; ++i) {
        cc_0[i] = new TCanvas(Form("c%d_0",i), Form("c%d_0",i), 800, 600);
        cc_0[i]->cd();
        h_0[i]->Fit(f_0[i], "", "", 3600, 5600);
        cc_1[i] = new TCanvas(Form("c%d_1",i), Form("c%d_1",i), 800, 600);
        cc_1[i]->cd();
        h_1[i]->Fit(f_1[i], "", "", 4400, 6500);
    }

    float c_0[n_times];
    float s_0[n_times];
    float c_1[n_times];
    float s_1[n_times];
    float R_0[n_times];
    float R_1[n_times];
    float err_c_0[n_times];
    float err_s_0[n_times];
    float err_c_1[n_times];
    float err_s_1[n_times];
    float err_R_0[n_times];
    float err_R_1[n_times];
    for (int i=0; i<n_times; ++i) {
        c_0[i] = f_0[i]->GetParameter(4);
        c_1[i] = f_1[i]->GetParameter(4);
        s_0[i] = f_0[i]->GetParameter(5);
        s_1[i] = f_1[i]->GetParameter(5);
        R_0[i] = (s_0[i]*2.355) / c_0[i];
        R_1[i] = (s_1[i]*2.355) / c_1[i];
        err_c_0[i] = f_0[i]->GetParError(4);
        err_c_1[i] = f_1[i]->GetParError(4);
        err_s_0[i] = f_0[i]->GetParError(5);
        err_s_1[i] = f_1[i]->GetParError(5);
        err_R_0[i] = R_0[i] * sqrt( pow(err_c_0[i]/c_0[i], 2.0) + pow(err_s_0[i]/s_0[i], 2.0) );
        err_R_1[i] = R_1[i] * sqrt( pow(err_c_1[i]/c_1[i], 2.0) + pow(err_s_1[i]/s_1[i], 2.0) );
    }

    ofstream ofile("Analysis/Day1/Day1_time.txt");
    ofile << "Centroids ch_0" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(5) << c_0[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(5) << err_c_0[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;

    ofile << "Sigmas ch_0" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(5) << s_0[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(5) << err_s_0[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;

    ofile << "Resolutions ch_0" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(6) << R_0[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(6) << err_R_0[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;



    ofile << "Centroids ch_1" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(5) << c_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(5) << err_c_1[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;

    ofile << "Sigmas ch_1" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(5) << s_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(5) << err_s_1[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;

    ofile << "Resolutions ch_1" << endl;
    for (int i=0; i<n_times; ++i) {
            ofile << '$' << '\t' << fixed << setw(10) << setprecision(6) << R_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(10) << setprecision(6) << err_R_1[i] << '\t' << '$' << '\t' << "\\\\" << endl;
    }
    ofile << endl;
}
// ********************************************************************************************
// ********************************************************************************************