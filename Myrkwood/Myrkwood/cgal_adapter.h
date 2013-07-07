//
//  cgal_adapter.h
//  Myrkwood
//
//  Created by Adel Ahmadyan on 12/1/12.
//  Copyright (c) 2012 University of Illinois at Urbana-Champaign. All rights reserved.
//

#ifndef __Myrkwood__cgal_adapter__
#define __Myrkwood__cgal_adapter__


#include <iostream>
#include <iterator>

#define CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel                  K;
//typedef K::Point_2                                                           Point_2;
typedef K::Iso_rectangle_2                                                   Iso_rectangle_2;
typedef K::Segment_2                                                         Segment_2;
typedef K::Ray_2                                                             Ray_2;
typedef K::Line_2                                                            Line_2;
typedef CGAL::Delaunay_triangulation_2<K>                                    DT;
typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT>                 AT;
typedef CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT> AP;
typedef CGAL::Voronoi_diagram_2<DT,AT,AP>                                    VD;

typedef AT::Site_2                                                           Site_2;
typedef AT::Point_2                                                          Point_2;

typedef VD::Locate_result                                                    Locate_result;
typedef VD::Vertex_handle                                                    Vertex_handle;
typedef VD::Face_handle                                                      Face_handle;
typedef VD::Halfedge_handle                                                  Halfedge_handle;
typedef VD::Ccb_halfedge_circulator                                          Ccb_halfedge_circulator;







//A class to recover Voronoi diagram from stream.
//Rays, lines and segments are cropped to a rectangle
//so that only segments are stored
struct Cropped_voronoi_from_delaunay{
    std::list<Segment_2> m_cropped_vd;
    Iso_rectangle_2 m_bbox;
    
    Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}
    Cropped_voronoi_from_delaunay(){}
    
    template <class RSL>
    void crop_and_extract_segment(const RSL& rsl){
        CGAL::Object obj = CGAL::intersection(rsl,m_bbox);
        const Segment_2* s=CGAL::object_cast<Segment_2>(&obj);
        if (s) m_cropped_vd.push_back(*s);
    }
    
    void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
    void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
    void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
};

void test();
void print_endpoint(Halfedge_handle e, bool is_src);

#endif /* defined(__Myrkwood__cgal_adapter__) */
