#include "TH1F.h"
#include "TF1.h"
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



void CalibrateHisto(TH1F *h_uncal, float m, float q) { //Re-scaling of axis, as in the slides

	int max_bin = h_uncal->GetNbinsX(); // This method returns the number of bins in x of the histogram
	float max_kev = h_uncal->GetBinCenter(max_bin)*m + q;
    float min_kev = h_uncal->GetBinCenter(1)*m + q;
	h_uncal->GetXaxis()->SetLimits(min_kev,max_kev);
	if (m!=1 && q!=0) //This means that I actually changed the calibration!
	    h_uncal->SetXTitle("keV");

}



vector<int> multiHisto(char *name_file, int draw=0, int column=0) {

    /*
    int nbin[8] = {
		100,
		30,
		30,
		100,
		100,
		100,
		100,
		30
	};
    int xmin[8] = {
		1000,
		1000,
		1000,
		1000,
		1000,
		1000,
		1000,
		2000
	};
    int xmax[8] = {
		20000,
		20000,
		20000,
		12000,
		6000,
		13000,
		9000,
		14000
	};
    */
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
    int xmin[8] = {
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0
	};
    int xmax[8] = {
		35000,
		35000,
		35000,
		35000,
		35000,
		35000,
		35000,
		35000
	};
    float par[8][2] = {
        {   0.1001,     57  },
        {   0.1777,     29  },
        {   0.1599,     17  },
        {   0.3027,     21  },
        {   0.2105,     24  },
        {   0.2517,     20  },
        {   0.3173,     -32 },
        {   0.1838,     23  }
    };
    /*
    int int_range[8][2] = {
        {   4000,       6000    },
        {   8000,       15000   },
        {   9500,       14500   },
        {   5500,       8500    },
        {   1800,       3100    },
        {   6400,       9600    },
        {   3000,       6000    },
        {   8000,       14000   }
    };
    int int_range_bin[8][2] = {
        {   14,     29  },
        {   11,     22  },
        {   14,     22  },
        {   40,     69  },
        {   15,     42  },
        {   44,     73  },
        {   24,     63  },
        {   16,     29  }
    };
    int int_range_bin[8][2] = {
        {   10,     100  },
        {   10,     100  },
        {   10,     100  },
        {   10,     100  },
        {    5,     100  },
        {    7,     100  },
        {   10,     100  },
        {   25,     100  }
    };
    */
    int int_range[8][2] = {
        {   1000,       35000   },
        {   1000,       35000   },
        {   1000,       35000   },
        {   1000,       35000   },
        {    500,       35000   },
        {    750,       35000   },
        {   1000,       35000   },
        {   2500,       35000   }
    };
    int int_range_bin[8][6] = {
        {  10,     10,		10,		10,		10,		200  },
        {  10,     10,		60,		10,		10,		200  },
        {  10,     10,		60,		10,		10,		200  },
        {   5,      5,		 5,		 5,		 5,		200  },
        {   3,      3,		 3,		 3,		 3,		200  },
        {   5,      5,		 5,		 5,		 5,		200  },
        {  40,     40,		40,		40,		40,		200  },
        {  20,     20,		60,		20,		20,		200  }
    };
    vector<int> integral(7);


    TBranch** inbranch_0 = new TBranch*[4];
    TBranch** inbranch_1 = new TBranch*[4];

    slimport_data_t indata_0[4];
    slimport_data_t indata_1[4];
	TFile *infile = new TFile(name_file);
	TTree *intree_0 = (TTree*)infile->Get("acq_tree_0");
    TTree *intree_1 = (TTree*)infile->Get("acq_tree_1");

    for (int i=0; i<4; ++i) {
        inbranch_0[i] = new TBranch;
        inbranch_1[i] = new TBranch;
        inbranch_0[i] = intree_0->GetBranch(Form("acq_ch%d",i));
        inbranch_1[i] = intree_1->GetBranch(Form("acq_ch%d",i));
        inbranch_0[i]->SetAddress(&indata_0[i].timetag);
        inbranch_1[i]->SetAddress(&indata_1[i].timetag);
    }

    TH1F** h = new TH1F*[8];
    for (int i=0; i<8; ++i) {
        h[i] = new TH1F(Form("h%d",i+1), Form("h%d",i+1), nbin[i], xmin[i], xmax[i]);
    }

    // histograms filling
    for (int i=0; i<4; ++i) {
    	for (int j=0; j<inbranch_0[i]->GetEntries(); ++j) {
    		inbranch_0[i]->GetEntry(j);
            // if (indata_0[i].timetag * 4e-9 <= 1000.0)
			h[i]->Fill(indata_0[i].qlong);
    	}
        for (int j=0; j<inbranch_1[i]->GetEntries(); ++j) {
    		inbranch_1[i]->GetEntry(j);
            // if (indata_1[i].timetag * 4e-9 <= 1000.0)
			h[4+i]->Fill(indata_1[i].qlong);
    	}
    }

    /*
    // calibrate histograms
    for (int i=0; i<8; ++i) {
        CalibrateHisto(h[i], par[i][0], par[i][1]);
    }
    */

    // integral
    for (int i=1; i<8; ++i) {
        integral[i-1] = h[i]->Integral(int_range_bin[i][column],int_range_bin[i][5]);
    }

    if (draw==1) {
        TCanvas* c1 = new TCanvas("c1", "c1", 800, 400);
	    c1->Divide(4,2);

        for (int i=0; i<8; ++i) {
            c1->cd(i+1);
            gPad->SetLogy();
            h[i]->Draw();
        }
    }

    return integral;
}



void analysisDay2() {
	// 1 = alum
	// 2 = lead
	// 3 = poly
	// 4 = iron
	// 5 = grph
	float rho_1 = 2.699;
	float rho_2 = 11.35;
	float rho_3 = 0.930;
	float rho_4 = 7.874;
	float rho_5 = 1.700;

	float Dx_1 = 2.0;
	float Dx_2 = 0.5;
	float Dx_3 = 2.0;
	float Dx_4 = 1.0;
	float Dx_5 = 2.0;

	float mu_1 = 0.08445;
	float mu_2 = 0.1614;
	float mu_3 = 0.09947;
	float mu_4 = 0.08414;
	float mu_5 = 0.08715;

    float exp_mu_1 = 0.;
    float exp_mu_2 = 0.;
    float exp_mu_3 = 0.;
    float exp_mu_4 = 0.;
    float exp_mu_5 = 0.;

    float err_exp_mu_1 = 0.;
    float err_exp_mu_2 = 0.;
    float err_exp_mu_3 = 0.;
    float err_exp_mu_4 = 0.;
    float err_exp_mu_5 = 0.;

    float denominator_err_exp_mu_1 = 0.;
    float denominator_err_exp_mu_2 = 0.;
    float denominator_err_exp_mu_3 = 0.;
    float denominator_err_exp_mu_4 = 0.;
    float denominator_err_exp_mu_5 = 0.;

    float final_compatibility_1;
    float final_compatibility_2;
    float final_compatibility_3;
    float final_compatibility_4;
    float final_compatibility_5;

    float int_ratio_1[7];
    float int_ratio_2[7];
    float int_ratio_3[7];
    float int_ratio_4[7];
    float int_ratio_5[7];

    float err_int_ratio_1[7];
    float err_int_ratio_2[7];
    float err_int_ratio_3[7];
    float err_int_ratio_4[7];
    float err_int_ratio_5[7];

	float atn_coef_1[7];
    float atn_coef_2[7];
    float atn_coef_3[7];
    float atn_coef_4[7];
    float atn_coef_5[7];

	float err_atn_coef_1[7];
    float err_atn_coef_2[7];
    float err_atn_coef_3[7];
    float err_atn_coef_4[7];
    float err_atn_coef_5[7];

	float compatibility_1[7];
	float compatibility_2[7];
	float compatibility_3[7];
	float compatibility_4[7];
	float compatibility_5[7];

    vector<int> integral_0_1 = multiHisto("DATA/Day2/I_0.root",0,0);
	vector<int> integral_0_2 = multiHisto("DATA/Day2/I_0.root",0,1);
	vector<int> integral_0_3 = multiHisto("DATA/Day2/I_0.root",0,2);
	vector<int> integral_0_4 = multiHisto("DATA/Day2/I_0.root",0,3);
	vector<int> integral_0_5 = multiHisto("DATA/Day2/I_0.root",0,4);

    vector<int> integral_1 = multiHisto("DATA/Day2/I_alum.root",0,0);
    vector<int> integral_2 = multiHisto("DATA/Day2/I_lead.root",0,1);
    vector<int> integral_3 = multiHisto("DATA/Day2/I_poly.root",0,2);
    vector<int> integral_4 = multiHisto("DATA/Day2/I_iron.root",0,3);
    vector<int> integral_5 = multiHisto("DATA/Day2/I_grph.root",0,4);

    vector<float> err_integral_0_1(7);
	vector<float> err_integral_0_2(7);
	vector<float> err_integral_0_3(7);
	vector<float> err_integral_0_4(7);
	vector<float> err_integral_0_5(7);

    vector<float> err_integral_1(7);
    vector<float> err_integral_2(7);
    vector<float> err_integral_3(7);
    vector<float> err_integral_4(7);
    vector<float> err_integral_5(7);

    for (int i=0; i<7; ++i) {
        err_integral_0_1[i] = sqrt(1.0*integral_0_1[i]);
		err_integral_0_2[i] = sqrt(1.0*integral_0_1[i]);
		err_integral_0_3[i] = sqrt(1.0*integral_0_1[i]);
		err_integral_0_4[i] = sqrt(1.0*integral_0_1[i]);
		err_integral_0_5[i] = sqrt(1.0*integral_0_1[i]);

        err_integral_1[i] = sqrt(1.0*integral_1[i]);
        err_integral_2[i] = sqrt(1.0*integral_2[i]);
        err_integral_3[i] = sqrt(1.0*integral_3[i]);
        err_integral_4[i] = sqrt(1.0*integral_4[i]);
        err_integral_5[i] = sqrt(1.0*integral_5[i]);
    }

    for (int i=0; i<7; ++i) {
        cout << integral_0_1[i] << '\t' << integral_1[i] << '\t' << integral_2[i] << '\t' << integral_3[i] << '\t' << integral_4[i] << '\t' << integral_5[i] << endl;
    }

    for (int i=0; i<7; ++i) {
        int_ratio_1[i] = (integral_1[i]*1.0) / (integral_0_1[i]*1.0);
        int_ratio_2[i] = (integral_2[i]*1.0) / (integral_0_2[i]*1.0);
        int_ratio_3[i] = (integral_3[i]*1.0) / (integral_0_3[i]*1.0);
        int_ratio_4[i] = (integral_4[i]*1.0) / (integral_0_4[i]*1.0);
        int_ratio_5[i] = (integral_5[i]*1.0) / (integral_0_5[i]*1.0);
    }

    for (int i=0; i<7; ++i) {
        err_int_ratio_1[i] = int_ratio_1[i] * sqrt( pow(err_integral_0_1[i]/integral_0_1[i], 2.0) + pow(err_integral_1[i]/integral_1[i], 2.0) );
        err_int_ratio_2[i] = int_ratio_2[i] * sqrt( pow(err_integral_0_2[i]/integral_0_2[i], 2.0) + pow(err_integral_2[i]/integral_2[i], 2.0) );
        err_int_ratio_3[i] = int_ratio_3[i] * sqrt( pow(err_integral_0_3[i]/integral_0_3[i], 2.0) + pow(err_integral_3[i]/integral_3[i], 2.0) );
        err_int_ratio_4[i] = int_ratio_4[i] * sqrt( pow(err_integral_0_4[i]/integral_0_4[i], 2.0) + pow(err_integral_4[i]/integral_4[i], 2.0) );
        err_int_ratio_5[i] = int_ratio_5[i] * sqrt( pow(err_integral_0_5[i]/integral_0_5[i], 2.0) + pow(err_integral_5[i]/integral_5[i], 2.0) );
    }

    float int_ratio_matrix[7][5];
    for (int i=0; i<7; ++i) {
        int_ratio_matrix[i][0] = int_ratio_1[i];
        int_ratio_matrix[i][1] = int_ratio_2[i];
        int_ratio_matrix[i][2] = int_ratio_3[i];
        int_ratio_matrix[i][3] = int_ratio_4[i];
        int_ratio_matrix[i][4] = int_ratio_5[i];
    }

    float err_int_ratio_matrix[7][5];
    for (int i=0; i<7; ++i) {
        err_int_ratio_matrix[i][0] = err_int_ratio_1[i];
        err_int_ratio_matrix[i][1] = err_int_ratio_2[i];
        err_int_ratio_matrix[i][2] = err_int_ratio_3[i];
        err_int_ratio_matrix[i][3] = err_int_ratio_4[i];
        err_int_ratio_matrix[i][4] = err_int_ratio_5[i];
    }

	for (int i=0; i<7; ++i) {
        atn_coef_1[i] = - log(int_ratio_1[i]) / (Dx_1*rho_1);
        atn_coef_2[i] = - log(int_ratio_2[i]) / (Dx_2*rho_2);
        atn_coef_3[i] = - log(int_ratio_3[i]) / (Dx_3*rho_3);
        atn_coef_4[i] = - log(int_ratio_4[i]) / (Dx_4*rho_4);
        atn_coef_5[i] = - log(int_ratio_5[i]) / (Dx_5*rho_5);
    }
    
	for (int i=0; i<7; ++i) {
        err_atn_coef_1[i] = err_int_ratio_1[i] / (int_ratio_1[i]*Dx_1*rho_1);
        err_atn_coef_2[i] = err_int_ratio_2[i] / (int_ratio_2[i]*Dx_2*rho_2);
        err_atn_coef_3[i] = err_int_ratio_3[i] / (int_ratio_3[i]*Dx_3*rho_3);
        err_atn_coef_4[i] = err_int_ratio_4[i] / (int_ratio_4[i]*Dx_4*rho_4);
        err_atn_coef_5[i] = err_int_ratio_5[i] / (int_ratio_5[i]*Dx_5*rho_5);
    }

	for (int i=0; i<7; ++i) {
        compatibility_1[i] = fabs(atn_coef_1[i] - mu_1) / err_atn_coef_1[i];
        compatibility_2[i] = fabs(atn_coef_2[i] - mu_2) / err_atn_coef_2[i];
        compatibility_3[i] = fabs(atn_coef_3[i] - mu_3) / err_atn_coef_3[i];
        compatibility_4[i] = fabs(atn_coef_4[i] - mu_4) / err_atn_coef_4[i];
        compatibility_5[i] = fabs(atn_coef_5[i] - mu_5) / err_atn_coef_5[i];
    }

    for (int i=0; i<7; ++i) {
        exp_mu_1 += atn_coef_1[i] / pow(err_atn_coef_1[i],2.0);
        exp_mu_2 += atn_coef_2[i] / pow(err_atn_coef_2[i],2.0);
        if (i!=5 && i!=6) exp_mu_3 += atn_coef_3[i] / pow(err_atn_coef_3[i],2.0);
        exp_mu_4 += atn_coef_4[i] / pow(err_atn_coef_4[i],2.0);
        exp_mu_5 += atn_coef_5[i] / pow(err_atn_coef_5[i],2.0);

        denominator_err_exp_mu_1 += 1.0 / pow(err_atn_coef_1[i],2.0);
        denominator_err_exp_mu_2 += 1.0 / pow(err_atn_coef_2[i],2.0);
        if (i!=5 && i!=6) denominator_err_exp_mu_3 += 1.0 / pow(err_atn_coef_3[i],2.0);
        denominator_err_exp_mu_4 += 1.0 / pow(err_atn_coef_4[i],2.0);
        denominator_err_exp_mu_5 += 1.0 / pow(err_atn_coef_5[i],2.0);
    }
    exp_mu_1 = exp_mu_1 / denominator_err_exp_mu_1;
    exp_mu_2 = exp_mu_2 / denominator_err_exp_mu_2;
    exp_mu_3 = exp_mu_3 / denominator_err_exp_mu_3;
    exp_mu_4 = exp_mu_4 / denominator_err_exp_mu_4;
    exp_mu_5 = exp_mu_5 / denominator_err_exp_mu_5;
    err_exp_mu_1 = sqrt(1.0 / denominator_err_exp_mu_1);
    err_exp_mu_2 = sqrt(1.0 / denominator_err_exp_mu_2);
    err_exp_mu_3 = sqrt(1.0 / denominator_err_exp_mu_3);
    err_exp_mu_4 = sqrt(1.0 / denominator_err_exp_mu_4);
    err_exp_mu_5 = sqrt(1.0 / denominator_err_exp_mu_5);

    final_compatibility_1 = fabs(exp_mu_1 - mu_1) / err_exp_mu_1;
    final_compatibility_2 = fabs(exp_mu_2 - mu_2) / err_exp_mu_2;
    final_compatibility_3 = fabs(exp_mu_3 - mu_3) / err_exp_mu_3;
    final_compatibility_4 = fabs(exp_mu_4 - mu_4) / err_exp_mu_4;
    final_compatibility_5 = fabs(exp_mu_5 - mu_5) / err_exp_mu_5;


    ofstream ofile("Analysis/Day2/Day2.txt");
	ofile << "Notation employed:" << endl;
	ofile << "[0] = bkg" << endl;
	ofile << "[1] = alum" << endl;
	ofile << "[2] = lead" << endl;
	ofile << "[3] = poly" << endl;
	ofile << "[4] = iron" << endl;
	ofile << "[5] = grph" << endl;
	ofile << endl;
    ofile << "Integral for [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << integral_0_1[i] << '\t' 
                                            << integral_1[i] << '\t' 
                                            << integral_2[i] << '\t' 
                                            << integral_3[i] << '\t' 
                                            << integral_4[i] << '\t' 
                                            << integral_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Error of integral for [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << err_integral_0_1[i] << '\t' 
                                            << err_integral_1[i] << '\t' 
                                            << err_integral_2[i] << '\t' 
                                            << err_integral_3[i] << '\t' 
                                            << err_integral_4[i] << '\t' 
                                            << err_integral_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Attenuation ratio for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << int_ratio_1[i] << '\t' 
                                            << int_ratio_2[i] << '\t' 
                                            << int_ratio_3[i] << '\t' 
                                            << int_ratio_4[i] << '\t' 
                                            << int_ratio_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Error of attenuation ratio for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << err_int_ratio_1[i] << '\t' 
                                            << err_int_ratio_2[i] << '\t' 
                                            << err_int_ratio_3[i] << '\t' 
                                            << err_int_ratio_4[i] << '\t' 
                                            << err_int_ratio_5[i] << endl;
    }
    ofile << endl << endl;
	ofile << "Attenuation coefficient for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << atn_coef_1[i] << '\t' 
                                            << atn_coef_2[i] << '\t' 
                                            << atn_coef_3[i] << '\t' 
                                            << atn_coef_4[i] << '\t' 
                                            << atn_coef_5[i] << endl;
    }
    ofile << endl << endl;
	ofile << "Errors of attenuation coefficient for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << err_atn_coef_1[i] << '\t' 
                                            << err_atn_coef_2[i] << '\t' 
                                            << err_atn_coef_3[i] << '\t' 
                                            << err_atn_coef_4[i] << '\t' 
                                            << err_atn_coef_5[i] << endl;
    }
    ofile << endl << endl;
	ofile << "Compatibility for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << compatibility_1[i] << '\t' 
                                            << compatibility_2[i] << '\t' 
                                            << compatibility_3[i] << '\t' 
                                            << compatibility_4[i] << '\t' 
                                            << compatibility_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Weighted attenuation coefficients for materials [1] [2] [3] [4] [5]:" << endl;
    ofile   << exp_mu_1 << '\t' 
            << exp_mu_2 << '\t' 
            << exp_mu_3 << '\t' 
            << exp_mu_4 << '\t' 
            << exp_mu_5 << endl;
    ofile << endl << endl;
    ofile << "Errors on weighted attenuation coefficients for materials [1] [2] [3] [4] [5]:" << endl;
    ofile   << err_exp_mu_1 << '\t' 
            << err_exp_mu_2 << '\t' 
            << err_exp_mu_3 << '\t' 
            << err_exp_mu_4 << '\t' 
            << err_exp_mu_5 << endl;
    ofile << endl << endl;
    ofile << endl << endl;
    ofile << "Final compatibility (on weighted attenuation coefficients) for materials [1] [2] [3] [4] [5]:" << endl;
    ofile   << err_exp_mu_1 << '\t' 
            << err_exp_mu_2 << '\t' 
            << err_exp_mu_3 << '\t' 
            << err_exp_mu_4 << '\t' 
            << err_exp_mu_5 << endl;
    ofile << endl << endl;
    ofile << "LaTeX: integral \\pm error for [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_0_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_0_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_3[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_4[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_5[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_5[i] << '\t' << '$' << "\t\\\\" << endl;

    }
    ofile << endl << endl;
    ofile << "LaTeX: I/I_0 \\pm error for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_3[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_4[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_5[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_5[i] << '\t' << '$' << "\t\\\\" << endl;
    }
	ofile << endl << endl;
    ofile << "LaTeX: attenuation coefficient \\pm error for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << atn_coef_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_atn_coef_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << atn_coef_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_atn_coef_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << atn_coef_3[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_atn_coef_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << atn_coef_4[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_atn_coef_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << atn_coef_5[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_atn_coef_5[i] << '\t' << '$' << "\t\\\\" << endl;
    }
	ofile << endl << endl;
    ofile << "LaTeX: compatibility for materials [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+2) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << compatibility_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << compatibility_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << compatibility_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << compatibility_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << compatibility_5[i] << '\t' << '$' << "\t\\\\" << endl;
    }
    ofile << endl << endl;
    ofile << "LaTeX: weighted attenuation coefficients for materials [1] [2] [3] [4] [5]:" << endl;
    ofile                  << '$' << '\t' << fixed << setw(8) << setprecision(4) << exp_mu_1 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_exp_mu_1 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << exp_mu_2 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_exp_mu_2 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << exp_mu_3 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_exp_mu_3 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << exp_mu_4 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_exp_mu_4 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << exp_mu_5 << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(4) << err_exp_mu_5 << '\t' << '$' << "\t\\\\" << endl;
    ofile << endl << endl;
    ofile << "LaTeX: final compatibility (on weighted attenuation coefficients) for materials [1] [2] [3] [4] [5]:" << endl;
    ofile                  << '$' << '\t' << fixed << setw(8) << setprecision(4) << final_compatibility_1 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << final_compatibility_2 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << final_compatibility_3 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << final_compatibility_4 << '\t' << '$' << '\t' 
            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(4) << final_compatibility_5 << '\t' << '$' << "\t\\\\";
    
    ofile.close();
}