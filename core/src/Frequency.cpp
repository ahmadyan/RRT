#include "Frequency.h"

// for testing/debugging only: this function generates a sin waveform
// that ramps from frequency f0 to f1 
vector<double> Frequency::generatoreSweepWaveform(double f0, double f1, double interval, int steps){
	vector<double> signal;
	for(int i=0;i<steps;i++){
        double delta = i / (float)steps;
        double t = interval * delta;
        double phase = 2 * M_PI * t * (f0 + (f1 - f0) * delta / 2);
        while (phase > 2 * M_PI) phase -= 2 * M_PI; // optional
		signal.push_back(3 * sin(phase));
        //printf("%f %f %f", t, phase * 180 / M_PI, 3 * sin(phase));
	}
	return signal;
}

//initilaize e-to-the-i-thetas for theta = 0..2PI in increments of 1/N
void Frequency::init_coeffs(){
    for (int i = 0; i < buffer_size; ++i) {
        double a = -2.0 * M_PI * i  / double(buffer_size);
        coeffs[i] = complex(cos(a)/* / N */, sin(a) /* / N */);
    }
    for (int i = 0; i < buffer_size; ++i) {
        double a = 2.0 * M_PI * i  / double(buffer_size);
        icoeffs[i] = complex(cos(a),sin(a));
    }
}


// initialize all data buffers
void Frequency::init(){
    // clear data
    for (int i = 0; i < buffer_size; ++i)
        in[i] = 0;
    // seed rand()
    srand(857);
    init_coeffs();
    oldest_data = newest_data = 0.0;
    idx = 0;
}

// simulating adding data to circular buffer
/*void Frequency::add_data(){
    oldest_data = in[idx];
    newest_data = in[idx] = complex(rand() / double(N));
}*/


// sliding dft
void Frequency::sdft(){
    complex delta = newest_data - oldest_data;
    int ci = 0;
    for (int i = 0; i < buffer_size; ++i) {
        freqs[i] += delta * coeffs[ci];
        if ((ci += idx) >= buffer_size)
            ci -= buffer_size;
    }
}

// sliding inverse dft
void Frequency::isdft(){
    complex delta = newest_data - oldest_data;
    int ci = 0;
    for (int i = 0; i < buffer_size; ++i) {
        freqs[i] += delta * icoeffs[ci];
        if ((ci += idx) >= buffer_size)
            ci -= buffer_size;
    }
}

// "textbook" slow dft, nested loops, O(N*N)
void Frequency::ft(){
    for (int i = 0; i < buffer_size; ++i) {
        freqs[i] = 0.0;
        for (int j = 0; j < buffer_size; ++j) {
            double a = -2.0 * M_PI * i * j / double(buffer_size);
            freqs[i] += in[j] * complex(cos(a),sin(a));
        }
    }
}

double Frequency::mag(complex& c){
    return sqrt(c.real() * c.real() + c.imag() * c.imag());
}

void Frequency::powr_spectrum(double *powr){
    for (int i = 0; i < buffer_size/2; ++i) {
        powr[i] = mag(freqs[i]);
    }

}
