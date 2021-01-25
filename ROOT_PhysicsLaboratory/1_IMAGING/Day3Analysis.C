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

int draw=0;


vector<int> multiHisto(char *name_file, int draw=0) {

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
		100,
		100,
		100,
		100,
		100,
		100,
		100,
		100
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
    int int_range_bin[8][2] = {
        {  10,     100  },
        {  25,     100  },
        {   5,     100  },
        {   5,     100  },
        {   4,     100  },
        {   5,     100  },
        {   5,     100  },
        {  10,     100  }
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
            if (indata_0[i].timetag * 4e-9 <= 1000.0) h[i]->Fill(indata_0[i].qlong);
    	}
        for (int j=0; j<inbranch_1[i]->GetEntries(); ++j) {
    		inbranch_1[i]->GetEntry(j);
            if (indata_1[i].timetag * 4e-9 <= 1000.0) h[4+i]->Fill(indata_1[i].qlong);
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
        integral[i-1] = h[i]->Integral(int_range_bin[i][0],int_range_bin[i][1]);
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



void attenuationRatio() {
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

    vector<int> integral_0 = multiHisto("DATA/Day3/I_bkg.root");
    vector<int> integral_1 = multiHisto("DATA/Day3/I_1.root");
    vector<int> integral_2 = multiHisto("DATA/Day3/I_2.root");
    vector<int> integral_3 = multiHisto("DATA/Day3/I_3.root");
    vector<int> integral_4 = multiHisto("DATA/Day3/I_4.root");
    vector<int> integral_5 = multiHisto("DATA/Day3/I_5.root");

    vector<float> err_integral_0(7);
    vector<float> err_integral_1(7);
    vector<float> err_integral_2(7);
    vector<float> err_integral_3(7);
    vector<float> err_integral_4(7);
    vector<float> err_integral_5(7);

    for (int i=0; i<7; ++i) {
        err_integral_0[i] = sqrt(1.0*integral_0[i]);
        err_integral_1[i] = sqrt(1.0*integral_1[i]);
        err_integral_2[i] = sqrt(1.0*integral_2[i]);
        err_integral_3[i] = sqrt(1.0*integral_3[i]);
        err_integral_4[i] = sqrt(1.0*integral_4[i]);
        err_integral_5[i] = sqrt(1.0*integral_5[i]);
    }

    for (int i=0; i<7; ++i) {
        cout << integral_0[i] << '\t' << integral_1[i] << '\t' << integral_2[i] << '\t' << integral_3[i] << '\t' << integral_4[i] << '\t' << integral_5[i] << endl;
    }

    for (int i=0; i<7; ++i) {
        int_ratio_1[i] = (integral_1[i]*1.0) / (integral_0[i]*1.0);
        int_ratio_2[i] = (integral_2[i]*1.0) / (integral_0[i]*1.0);
        int_ratio_3[i] = (integral_3[i]*1.0) / (integral_0[i]*1.0);
        int_ratio_4[i] = (integral_4[i]*1.0) / (integral_0[i]*1.0);
        int_ratio_5[i] = (integral_5[i]*1.0) / (integral_0[i]*1.0);
    }

    for (int i=0; i<7; ++i) {
        err_int_ratio_1[i] = int_ratio_1[i] * sqrt( pow(err_integral_0[i]/integral_0[i], 2.0) + pow(err_integral_1[i]/integral_1[i], 2.0) );
        err_int_ratio_2[i] = int_ratio_2[i] * sqrt( pow(err_integral_0[i]/integral_0[i], 2.0) + pow(err_integral_2[i]/integral_2[i], 2.0) );
        err_int_ratio_3[i] = int_ratio_3[i] * sqrt( pow(err_integral_0[i]/integral_0[i], 2.0) + pow(err_integral_3[i]/integral_3[i], 2.0) );
        err_int_ratio_4[i] = int_ratio_4[i] * sqrt( pow(err_integral_0[i]/integral_0[i], 2.0) + pow(err_integral_4[i]/integral_4[i], 2.0) );
        err_int_ratio_5[i] = int_ratio_5[i] * sqrt( pow(err_integral_0[i]/integral_0[i], 2.0) + pow(err_integral_5[i]/integral_5[i], 2.0) );
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
        for (int j=0; j<5; ++j) {
            cout << int_ratio_matrix[i][j] << "\t\t";
        }
        cout << endl;
    }

    TH2D* hh1 = new TH2D("hh1", "hh1", 5, 0.5, 5.5, 7, -8.5, -1.5);
    float x,y;
    for (int i=0; i<5; ++i) {
        for (int j=0; j<7; ++j) {
            x = i+1;
            y = -j-2;
            hh1->Fill(x,y,int_ratio_matrix[j][i]);
        }
    }
    gStyle->SetPalette(54);
    hh1->SetContour(1000);
    hh1->GetZaxis()->SetRangeUser(0.3, 1.2);
    hh1->Draw("colz");
    hh1->GetXaxis()->SetTitle("Position");
    hh1->GetYaxis()->SetTitle("Detector");
    hh1->GetXaxis()->SetTickLength(0.01);
    hh1->GetYaxis()->SetTickLength(0.01);
    hh1->SetTitle("Image");
    hh1->GetXaxis()->SetNdivisions(1005);
    gStyle->SetOptStat(0);
    gPad->Update();
    gPad->SaveAs("Analysis/Day3/Image.pdf");


    ofstream ofile("Analysis/Day3/Day3.txt");
    ofile << "Integral in positions [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << integral_0[i] << '\t' 
                                            << integral_1[i] << '\t' 
                                            << integral_2[i] << '\t' 
                                            << integral_3[i] << '\t' 
                                            << integral_4[i] << '\t' 
                                            << integral_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Error of integral in positions [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << err_integral_0[i] << '\t' 
                                            << err_integral_1[i] << '\t' 
                                            << err_integral_2[i] << '\t' 
                                            << err_integral_3[i] << '\t' 
                                            << err_integral_4[i] << '\t' 
                                            << err_integral_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Attenuation ratio in positions [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << int_ratio_1[i] << '\t' 
                                            << int_ratio_2[i] << '\t' 
                                            << int_ratio_3[i] << '\t' 
                                            << int_ratio_4[i] << '\t' 
                                            << int_ratio_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "Error of attenuation ratio in positions [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << err_int_ratio_1[i] << '\t' 
                                            << err_int_ratio_2[i] << '\t' 
                                            << err_int_ratio_3[i] << '\t' 
                                            << err_int_ratio_4[i] << '\t' 
                                            << err_int_ratio_5[i] << endl;
    }
    ofile << endl << endl;
    ofile << "LaTeX: integral \pm error in positions [0] [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_0[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_0[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_3[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_4[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(0) << 1.0*integral_5[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(0) << err_integral_5[i] << '\t' << '$' << "\t\\\\" << endl;

    }
    ofile << endl << endl;
    ofile << "LaTeX: I/I_0 \pm error in positions [1] [2] [3] [4] [5]:" << endl;
    for (int i=0; i<7; ++i) {
        ofile << Form("D%d",i+1) << '\t'    << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_1[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_1[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_2[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_2[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_3[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_3[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_4[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_4[i] << '\t' << '$' << '\t' 
                                            << '&' << '\t' << '$' << '\t' << fixed << setw(8) << setprecision(3) << int_ratio_5[i] << '\t' << "\\pm" << '\t' << fixed << setw(8) << setprecision(3) << err_int_ratio_5[i] << '\t' << '$' << "\t\\\\" << endl;
    }
    ofile.close();
}