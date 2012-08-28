

Monitor Examples:
	Monitor* m0 = new Monitor(rrt, CONSTANT);
	m0->push_back(0.5);

	Monitor* m1 = new Monitor(rrt, ANALOG_BINARY_ADD );		// m1 = (x_2 + 0.5)
	m1->push_back(2);
	m1->push_back(m0);		

	Monitor* m2 = new Monitor(rrt, ANALOG_SHIFT ) ;			// m2 = x[t-(10e-7)]
	m2->push_back(10e-7);

	Monitor* m3 = new Monitor(rrt, ANALOG_BINARY_MUL) ;		// m3 = (x_2 + 0.5) * x[t-(10e-7)]
	m3->push_back(m1); 
	m3->push_back(m2);