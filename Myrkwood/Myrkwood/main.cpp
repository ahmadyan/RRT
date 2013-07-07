
#include <iostream>
#include "config.h"
#include "circuit.h"
#include "Algorithm.h"
#include "set.h"
#include "trace.h"
#include "plot.h"
#include "Voronoi.h"
#include "utility.h"
#include "Venice.h"


Circuit* loadCircuit(Configuration* config){
    Circuit* circuit;
    int k = 1;                  config->getParameter("edu.uiuc.crhc.core.system", &k);
    switch(k){
        case 1:
            circuit = new Circuit(tunneldiode, config);
            break;
        case 2:
            circuit = new Circuit(vanderpol, config);
            break;
        case 3:
            circuit = new Circuit(ringmodulator, config);
            break;
        case 4:
            circuit = new Circuit(dae_amplifier, config);
    }
    return circuit ;
}

void simulation_trace(){
    srand((int)time(0));
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Circuit* circuit = new Circuit(ringmodulator, config);
    Plotter* plotter = new Plotter(gnuPlot, config);
    //Trace* trace = circuit->simulate(config, false);
    //Trace* trace = circuit->backward_simulate(config);
    Trace* trace = circuit->sim_test();
    plotter->drawTrace(trace, 0, 1);
    //plotter->drawTrace(trace, 1e-3, 1);
	plotter->close();
    
    delete plotter;
    delete trace;
    delete config;
    delete circuit;    
}

void venice_experiment(){
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Plotter* plotter = new Plotter(gnuPlot, config);
    venice_demo();
    plotter->execute(venice_draw());
    delete plotter;
    delete config;
}

void verification_experiment(){
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Plotter* plotter = new Plotter(gnuPlot, config);
    Circuit* circuit = loadCircuit(config);
    utility::tick();
    
    Algorithm* alg = new Algorithm(circuit, config);
    Trace* trace = circuit->simulate_stiff(true);
	string reachableSet = alg->compute_reachable_set(trace);
	//Trace* trace = alg->generate_counter_example();
    utility::tock();
    
    Voronoi* v = new Voronoi(circuit);
    //v->generateDelaunay(config, circuit, alg->getForest() );
    v->generateVoronoi(config, circuit, alg->getForest() , 0, 1);
    
    plotter->execute(v->drawDelaunay());
    plotter->execute(v->drawVoronoi(alg->getForest()));
    //plotter->execute(circuit->generateVectorField());
    //plotter->execute(reachableSet);
    
}


void falsification_experiment(){
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Plotter* plotter = new Plotter(gnuPlot, config);
    Circuit* circuit = loadCircuit(config);
    utility::tick();
    Algorithm* alg = new Algorithm(circuit, config);
    Trace* trace = circuit->simulate_dae();
    plotter->execute( alg->generate_counter_example(trace) );
    
    utility::tock();
}

void classic_rrt_experiment(){
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Plotter* plotter = new Plotter(gnuPlot, config);
    Circuit* circuit = loadCircuit(config);
    //Trace* trace = circuit->simulate_stiff(true);
    //Trace* trace = circuit->sim_test();
    //Trace* trace = circuit->backward_simulate(config);
    Trace* trace = circuit->simulate_dae();
    utility::tick();
    Algorithm* alg = new Algorithm(circuit, config);
    alg->forward_exploration(trace);
    plotter->execute( alg->getForest()->draw(0,1)) ;
    //plotter->execute( alg->getForest()->draw(0)) ;
    utility::tock();
}

void     radau_experiment(){
    Configuration* config = new Configuration("/Users/adel/Dropbox/git/rrt/Myrkwood/Myrkwood/reach.conf") ;
    Plotter* plotter = new Plotter(gnuPlot, config);
    Circuit* circuit = loadCircuit(config);
    Trace* t = circuit->simulate_dae();
    cout << "sim finished" << endl ;
//    plotter->drawTrace(t, 0, 1);
    plotter->drawTrace(t, 0.1, 7);
}

int main(int argc, const char * argv[]){
    srand((int)time(0));
    //simulation_trace();
    //backward_simulation_trace();
    //venice_experiment();
    //verification_experiment();
    falsification_experiment();
    //classic_rrt_experiment();
    //radau_experiment();
    return 0;
}
