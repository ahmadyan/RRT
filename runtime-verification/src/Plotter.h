#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "transient.h"
#include "Rrt.h"
using namespace std;
class Plotter{
    string gnuPlotPath ;
public:
    Plotter(string path);
    ~Plotter();
    void plot();
    void plotRRT(string name, string title, string output, RRT rrt, string, string, string);
    void plot2D(string name, string title, string output, vector<transient> simData);
    void plotTrace(RRT rrt, int v1, int v2, int, double, double);
    void wait_for_key();
    
};
