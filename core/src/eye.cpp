#include "eye.h"


EyeDiagram::EyeDiagram(double fc, int _sampleRate, int w){
	freq = fc;
	sampleRate = _sampleRate;

	period = 1 / freq;
	window = w;

	palpebraSuperior = vector< vector<node*> >(window);

	maxSuperior = new double[window];
	minSuperior = new double[window];
	
	palpebraInferior = vector< vector<node*> >(window);
	maxInferior = new double[window];	
	minInferior = new double[window];


}

EyeDiagram::~EyeDiagram(){
	delete maxSuperior;
	delete maxInferior;
	delete minSuperior;
	delete minInferior;
}

void EyeDiagram::push(){

}