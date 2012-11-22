#pragma once

#include "set.h"
#include "trace.h"
class Algorithm{
public:
	Algorithm();
	~Algorithm();
	Set* compute_reachable_set();
	Trace* generate_counter_example();
};