struct evt_data_t {
	UInt_t evcounter;
	UInt_t timestamp;
	UShort_t samples[8];
};

void writetree(const char *fileName = "exampletree.root") {
	
	evt_data_t out_data; // declare the struct
	
	TFile *output_file = new TFile(fileName,"RECREATE"); // declare an output file
    	output_file->cd();

	TTree *output_tree = new TTree("datatree","Tree di esempio"); // declare the tree

	// declare the branch; the leaflist corresponds to the struct.
	TBranch *branch = output_tree->Branch("evbranch", &out_data.evcounter, "counter/i:time:samples[8]/s");

	// fill the tree with invented data and print them
	cout << "Ev. #\tTime\tSamples" << endl;
	for(int i=0; i<18; i++){
		out_data.evcounter = i;
		cout << i << "\t";
		out_data.timestamp = 5+i*2+i*i;
		cout << 5+i*2+i*i << "\t";
		for(int k=0; k<8; k++){
			out_data.samples[k] = i+k*2;
			cout << out_data.samples[k] << " ";
		}
		cout << endl;
		branch->Fill();
	}

	// write, print the tree and close the file
	output_tree->Write();
	output_tree->Print();
	cout << endl;
	output_file->Close();

}
