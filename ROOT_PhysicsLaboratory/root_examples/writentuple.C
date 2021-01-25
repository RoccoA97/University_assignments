void writentuple(const char *fileName = "example_ntuple.root") {

	float	evcounter;
	float	timestamp;
	float	sample;

	// declare an output file
	TFile *output_file = new TFile(fileName,"RECREATE"); output_file->cd();

	// declare a N-tuple
	TNtuple *output_ntuple = new TNtuple("datantuple","Example N-tupla","counter:time:sample");

	// fill the tree with 18 fake events and print them
	cout << "Ev. #\tTime\tSample" << endl;
	for (int i=0; i<18; i++) {
		evcounter = i;					cout << i << "\t";
		timestamp = 5+i*2+i*i;		cout << 5+i*2+i*i << "\t";
		sample = i*35+44.5;			cout << i*35+44.5 << endl;
		output_ntuple->Fill(evcounter,timestamp,sample);
	}

	// write, print the tree and close the file
	output_ntuple->Write();
	output_ntuple->Print(); cout << endl;
	output_file->Close();

}
