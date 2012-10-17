//============================================================================
// Name        : Reach.cpp
// Author      : Adel Ahmadyan
// Version     :
// Copyright   : (C) 2012 UIUC.
// Description :
//============================================================================

#include <iostream>
#include <search.h>
#include <malloc.h>
#include "config.h"
#include "plot.h"
#include "utility.h"
#include "VoronoiDiagramGenerator.h"

using namespace std;

int main(int argc, char** argv) {
	if(argc<2){
		cout << "Configuration file is not specified." << endl ;
		return -1;
	}

	Configuration* config = new Configuration(argv[1]);
	Plotter* plot = new Plotter(gnuPlot, config);

	cout << "generating voronoi diagram..." << endl ;
	float* xValues = new float[100];
	float* yValues = new float[100];
	long count = 100;
	for(int i=0;i<count;i++){
		xValues[i] = utility::unifRand(-10, 10);
		yValues[i] = utility::unifRand(-10, 10);
	}

	VoronoiDiagramGenerator vdg;
	vdg.generateVoronoi(xValues,yValues,count, -10,10,-10,10,0);

	cout << "here it comes..." << endl ;
	vdg.resetIterator();

	float x1,y1,x2,y2;

	printf("\n-------------------------------\n");
	while(vdg.getNext(x1,y1,x2,y2))
	{
		plot->drawLine(x1,y1,x2, y2);
		printf("GOT Line (%f,%f)->(%f,%f)\n",x1,y1,x2, y2);

	}

	plot->close();
	delete config;
	delete plot;
	return 0;
}
