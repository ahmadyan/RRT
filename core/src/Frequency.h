#pragma once
#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <complex>

#include <iostream>
#include <time.h>
#include <float.h>

using namespace std;

class Frequency{
public:
	//transformations

	//utility methosd
	vector<double> generatoreSweepWaveform(double f0, double f1, double interval, int steps);

	
	typedef std::complex<double> complex;

	int buffer_size;
	// input signal
	complex* in; //in[N]

	// frequencies of input signal after ft
	// Size increased by one because the optimized sdft code writes data to freqs[N]
	complex* freqs; //freqs[N+1];

	// output signal after inverse ft of freqs
	complex* out1; //out1[N];
	complex* out2; //out2[N];

	// forward coeffs -2 PI e^iw -- normalized (divided by N)
	complex* coeffs; //coeffs[N];
	// inverse coeffs 2 PI e^iw
	complex* icoeffs; //icoeffs[N];

	// global index for input and output signals
	int idx;


	// these are just there to optimize (get rid of index lookups in sdft)
	complex oldest_data, newest_data;
	void init_coeffs();
	void init();
	void sdft();
	void isdft();
	void ft();
	double mag(complex& c);
	void powr_spectrum(double *powr);
};