// declaration of a struct
struct event_t {
	int	time;
	float	energy;
	int	samples[256];
};

event_t readeventStruct(int evNum) {
	// function pretending to read an event from file
	// energy is 0
	event_t evt;
	evt.time = evNum*4;
	evt.energy = 0;
	for (int i=0; i<256; i++)
	   evt.samples[i] = i*evNum/16448;

	return evt;
}

event_t analyzeeventStruct(event_t evt) {
	// declaring an object of the type that will be sent as output
	event_t evtemp;

    /* I initialize the new event as the input one but I can't
	   use directly the input one because I passed it as
	   value and then it is available only for reading data,
	   I can't write the new values of energy. */
    
    evtemp.time	= evt.time;
	evtemp.energy	= evt.energy;
	for (int i=0; i<256; i++)
	   evtemp.samples[i] = evt.samples[i];

	// Function that calculates the energy as the sum of the measured samples
	for (int i=0; i<256; i++)
	   evtemp.energy += evtemp.samples[i];

	return evtemp;
}

void analyzeeventStructPuntatore(event_t *evt) {
	// Function that calculates the energy as the sum of the measured samples
	for (int i=0; i<256; i++)
        evt->energy += evt->samples[i];
}

void struct_example() {

	// we foresee maximum 1024 events

	// declaration of struct variable (in fact, an array of 1024 structs)
	event_t	event[1024];

	// reading events
	for (int ev=0; ev<1024; ev++) {

		// struct
		event[ev] = readeventStruct(ev);
		event[ev] = analyzeeventStruct(event[ev]);
		cout << "Evt #" << ev << ": energy = " << event[ev].energy << endl;
//		analyzeeventStructPuntatore(&event[ev]);

	}

}

