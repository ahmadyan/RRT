#pragma once

#include "set.h"
#include "trace.h"
#include "config.h"
#include "circuit.h"
#include "Forest.h"
class Algorithm{
    Configuration* config ;
    Circuit* circuit ;
    Forest* rrf;
public:
	Algorithm(Circuit*, Configuration*);
	~Algorithm();
	string compute_reachable_set(Trace*);
	string generate_counter_example(Trace* trace);
    string forward_exploration(Trace*);
    Forest* getForest() ;
};