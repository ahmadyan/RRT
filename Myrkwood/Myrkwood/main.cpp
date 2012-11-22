#include <iostream>
#include "config.h"
#include "circuit.h"
#include "Algorithm.h"
#include "set.h"
#include "trace.h"
#include "plot.h"

int main(int argc, const char * argv[]){
	Configuration* config = new Configuration("reach.conf") ;
	Circuit* circuit = new Circuit(vanderpol);
	Algorithm* alg = new Algorithm();
	Set* reachableSet = alg->compute_reachable_set();
	Trace* trace = alg->generate_counter_example();
	Plotter* plotter = new Plotter(gnuPlot, config);
	plotter->close();
	system("pause");
    return 0;
}

