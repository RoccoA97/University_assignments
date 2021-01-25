void readntuple(char *infilename) {

	float	evcounter;
	float	timestamp;
	float	sample;

	// declare an input file
	TFile *input_file = new TFile(infilename);

	// read the n-tuple in the input file
	TNtuple *input_ntuple = (TNtuple*)input_file->Get("datantuple");

	int nEntries = (int)input_ntuple->GetEntries();

	cout << "Acquisition loaded; number of entries: " << nEntries << endl;

	input_ntuple->SetBranchAddress("counter",&evcounter);
	input_ntuple->SetBranchAddress("time",&timestamp);
	input_ntuple->SetBranchAddress("sample",&sample);

	for (int i=0; i< nEntries; i++) {
		input_ntuple->GetEntry(i);
		cout << evcounter << "\t" << timestamp << "\t" << sample << endl;
	}

}

