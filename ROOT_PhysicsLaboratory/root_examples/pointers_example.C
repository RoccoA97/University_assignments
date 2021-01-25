void readEvent(int evNum, int *time, float *energy, int *samples) {
	// function pretending to read an event from file
	// energy is 0

	(*time)	= evNum*4;
	(*energy) = 0;
	for (int i=0; i<256; i++) 
	   samples[i] = i*evNum/16448;
}

void analyzeEvent(int time, float *energy, int samples[256]) {
	// function calculating the energy as sum of the measured samples
	for (int i=0; i<256; i++)
	   (*energy) += samples[i];
}

void pointers_example() {

	// we foresee a maximum of 1024 events

	// declaration of separate variables
	int	time[1024];
	float	energy[1024];
	int	samples[1024][256];

	// analyzing with a for loop over all the events
	for (int ev=0; ev<1024; ev++) {
		// normal variables
	   readEvent(ev, &time[ev], &energy[ev], samples[ev]);
	   analyzeEvent(time[ev], &energy[ev], samples[ev]);
	   cout << "Evt #" << ev << "; time = " << time[ev] << " \u03BCs; energy = " << energy[ev] << " eV" << endl;
	}

}

int increase(int *value, int addend) {
	(*value) += addend;
   return (*value);
}

int main() {
	int a = 3;
	a = increase(&a,4);
	return a;		// a is 7
}
