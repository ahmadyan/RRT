//
//  Voronoi.h
//  Myrkwood
//
//  Created by Adel Ahmadyan on 12/1/12.
//  Copyright (c) 2012 University of Illinois at Urbana-Champaign. All rights reserved.
//

#ifndef __Myrkwood__Voronoi__
#define __Myrkwood__Voronoi__

#include <iostream>
#include "cgal_adapter.h"
//#include "VoronoiDiagramGenerator.h"
#include "config.h"
#include "circuit.h"
#include "Forest.h"

class Voronoi{
    //VoronoiDiagramGenerator vdg;
    DT dt2;
    VD vd;
    Cropped_voronoi_from_delaunay vor;
    bool forest_enabled ;
public:
    Voronoi(Circuit* circuit);
    ~Voronoi();
    void generateDelaunay(Configuration* _config, Circuit* _circuit, Forest* _forest);
    void generateDelaunay(Configuration* _config, Circuit* _circuit, vector<geometry::Point*> p);
    string drawVoronoi(Forest* forest);
    string drawDelaunay();
    void generateVoronoi(Configuration* _config, Circuit* _circuit, vector<geometry::Point*> p);
    void generateVoronoi(Configuration* _config, Circuit* _circuit, Forest* forest, int, int);
};

#endif /* defined(__Myrkwood__Voronoi__) */
