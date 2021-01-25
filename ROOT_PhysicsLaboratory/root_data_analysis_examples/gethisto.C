//Digitizer data from the LAB

struct slimport_data_t {
	ULong64_t	timetag; //time stamp
	UInt_t		baseline;
	UShort_t	qshort; //integration with shorter time
	UShort_t	qlong; //integration with longer time
	UShort_t	pur;
	UShort_t	samples[4096];
};

TH1F* getHistoFromTree(char *name_file, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1F *h_spectrum = new TH1F("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}

TH1F* getHistoForChannelFromTree(char *name_file, short dgtz, short chan, int numBins, double minX, double maxX) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get(Form("acq_tree_%d",dgtz));
	TBranch *inbranch = intree->GetBranch(Form("acq_ch%d",chan));
	inbranch->SetAddress(&indata.timetag);
	TH1F *h_spectrum = new TH1F("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		h_spectrum->Fill(indata.qlong);
	}
	// return
	return h_spectrum;
}

void multiCanvas(char *name_file) {
	TCanvas* c1 = new TCanvas("c1");
	c1->Divide(4,2);
	for (int i=0; i<2; ++i) {
		for (int j=0; j<4; ++j) {
			if (i==0) {
				c1->cd(1+j);
				getHistoForChannelFromTree(name_file,i,j,100,1000,20000)->Draw();
			}
			else{
				c1->cd(5+j);
				getHistoForChannelFromTree(name_file,i,j,100,1000,20000)->Draw();
			}
		}
	}
}

TH1F* getHistoWithFilter(char *name_file, int numBins, double minX, double maxX, double lowThr = 0, double highThr = 999999) {
	// variables
	slimport_data_t indata;
	TFile *infile = new TFile(name_file);
	TTree *intree = (TTree*)infile->Get("acq_tree_0");
	TBranch *inbranch = intree->GetBranch("acq_ch0");
	inbranch->SetAddress(&indata.timetag);
	TH1F *h_spectrum = new TH1F("h_spectrum","Total spectrum",numBins,minX,maxX);
	// histogram filling
	for (int i=0; i<inbranch->GetEntries(); i++) {
		inbranch->GetEntry(i);
		if (indata.qlong>lowThr && indata.qlong<highThr) {
			h_spectrum->Fill(indata.qlong);
		}
	}
	// return
	return h_spectrum;
}
