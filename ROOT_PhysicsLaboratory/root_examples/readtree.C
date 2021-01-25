struct evt_data_t {
	UInt_t evcounter;
	UInt_t timestamp;
	UShort_t samples[8];
};

void readtree(const char *fileName = "exampletree.root") {
	
	evt_data_t	inc_data;	// declare the struct

	TFile *input_file = new TFile(fileName,"READ"); // open the file in reading mode
	TTree *input_tree = (TTree*)input_file->Get("datatree"); //read the tree

        // TFile->Get() returns a generic object which must
	// be explicitly converted putting (TTree*) before

	input_tree->Print(); //check on the screen the structure of the tree

	// I get the branch and associate it to the struct
	TBranch *branch = input_tree->GetBranch("evbranch");
	branch->SetAddress(&inc_data.evcounter);

	// I read the entries (leaves) and print the data
	cout << "Ev. #\tTime\tSamples" << endl;
	for (int i=0; i<18; i++) {
		branch->GetEntry(i); // this fills the struct
		cout << "#" << inc_data.evcounter << "\t";
		cout << inc_data.timestamp << "\t";
		for (int k=0; k<8; k++)
		    cout << inc_data.samples[k] << " ";
		cout << endl;
	}

	//input_file->Close();
}
