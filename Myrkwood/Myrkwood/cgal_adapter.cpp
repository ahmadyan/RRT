//
//  cgal_adapter.cpp
//  Myrkwood
//
//  Created by Adel Ahmadyan on 12/1/12.
//  Copyright (c) 2012 University of Illinois at Urbana-Champaign. All rights reserved.
//

#include "cgal_adapter.h"
void test(){
    
    std::vector<Point_2> points;
    points.push_back(Point_2(0,0));
    points.push_back(Point_2(1,1));
    points.push_back(Point_2(0,1));
    
    DT dt2;
    //insert points into the triangulation
    dt2.insert(points.begin(),points.end());
    //construct a rectangle
    Iso_rectangle_2 bbox(-1,-1,2,2);
    Cropped_voronoi_from_delaunay vor(bbox);
    //extract the cropped Voronoi diagram
    dt2.draw_dual(vor);
    //print the cropped Voronoi diagram as segments
    
    
    std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(),
              std::ostream_iterator<Segment_2>(std::cout,"\n"));
    
}


void print_endpoint(Halfedge_handle e, bool is_src) {
    std::cout << "\t";
    if ( is_src ) {
        if ( e->has_source() )  std::cout << e->source()->point() << std::endl;
        else  std::cout << "point at infinity" << std::endl;
    } else {
        if ( e->has_target() )  std::cout << e->target()->point() << std::endl;
        else  std::cout << "point at infinity" << std::endl;
    }
}

