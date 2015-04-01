#pragma once
#include <string>
#include <iostream>
#include <vector>
#include "transient.h"
#include "Rrt.h"
using namespace std;
#ifdef _WIN32
	const string gnuPlot = "C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist" ;
#elif  __linux__
	const string gnuPlot = "/usr/bin/gnuplot" ;
#elif  __APPLE__
	//const string gnuPlot = "/Applications/Gnuplot.app/Contents/Resources/bin/gnuplot" ; //Old gnuplot path (v4.0) on my Mac, does not support drawing a polygin, on aqua.
	const string gnuPlot = "/usr/local/bin/gnuplot"; //current gnuPlot path, v4.6 on x11
#else
	const string gnuPlot = "/usr/local/bin/gnuplot";
#error Undefined operating system/Unknown gnuplot path
#endif



class Plotter{
    string gnuPlotPath ;
	FILE *gnuplotPipe;
	 bool closed;
public:
    Plotter(string path);
    ~Plotter();

	void emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel,string zlabel);
    void emptyPlot(string title, double xmin, double xmax, double ymin, double ymax );

    void plot();
    void plotRRT(string name, string title, string output, RRT rrt, string, string, string);
    void plot2D(string name, string title, string output, vector<transient> simData);
	void plotTrace(RRT rrt, int v1, int v2, int, double, double, string);
    void wait_for_key();
    void drawTrace(vector<double>, double t);
	void drawArray(vector< vector<double> > trace, int, int);
    //void drawTrace(Trace* trace, int, int);
    //void drawTrace(Trace* trace, double t_max, int index);
    void drawLine(double iFromX, double iFromY, double iToX, double iToY);
    void drawLine(double iFromX, double iFromY, double iFromZ, double iToX, double iToY, double iToZ);
    void saveToPdf(string path);
    void waitForKey();
    void close();
	void execute(string str);
	void drawLine(double iFromX, double iFromY, double iToX, double iToY, string color);
	void test();
};
