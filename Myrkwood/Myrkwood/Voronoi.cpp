//
//  Voronoi.cpp
//  Myrkwood
//
//  Created by Adel Ahmadyan on 12/1/12.
//  Copyright (c) 2012 University of Illinois at Urbana-Champaign. All rights reserved.
//

#include "Voronoi.h"
#include <sstream>

using namespace std;
Voronoi::Voronoi(Circuit* circuit){
    Iso_rectangle_2 bbox(circuit->getMin(0),circuit->getMin(1),circuit->getMax(0),circuit->getMax(1));
    vor = Cropped_voronoi_from_delaunay(bbox);
}

Voronoi::~Voronoi(){}

void Voronoi::generateDelaunay(Configuration* config, Circuit* circuit, Forest* forest){
    forest_enabled = true;
    std::vector<Point_2> points;
    for(int i=0;i<forest->getSize();i++){
        points.push_back(  Point_2(forest->getNode(i)->getData(0), forest->getNode(i)->getData(1) ) );
    }
    dt2.insert(points.begin(),points.end());
    dt2.draw_dual(vor); //replaces rays with bounded segments
}

void Voronoi::generateDelaunay(Configuration* config, Circuit* circuit, vector<geometry::Point*> p){
    forest_enabled = false ;
    std::vector<Point_2> points;
    for(int i=0;i<p.size();i++){
        points.push_back(Point_2(p[i]->getData(0), p[i]->getData(1)));
    }
    
    for(int i=0;i<points.size();i++){
        cout << points[i] << endl ;
    }
    dt2.insert(points.begin(),points.end());
    
//    cout << "dt2" << endl ;
//    int i=0;

//    for(DT::Vertex_iterator it=dt2.vertices_begin(); it!=dt2.vertices_end(); it++){
//        cout << "i=" << i << endl ;
//        DT::Vertex_circulator vc = it->incident_vertices(), done(vc); ;
//        do{
//            cout << vc->point() << endl ;
//        }while(++vc!=done);
 //       i++;
 //       cout << "-------" << endl ;
 //   }
    
    

    dt2.draw_dual(vor); //replaces rays with bounded segments
    //std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(), std::ostream_iterator<Segment_2>(std::cout,"\n"));
}



/*
 for(VD::Vertex_iterator it = vd.vertices_begin(); it != vd.vertices_end(); it++){
 std::cout << it->point()  << std::endl;
 }
 cout << "-------" << endl ;
 
 std::cout << "\n===========\nNumber of Voronoi Sites : " << vd.number_of_vertices()  << std::endl;
 for(VD::Site_iterator it = vd.sites_begin(); it!=vd.sites_end(); it++){
 cout << *it << endl ;
 
 }
 std::cout << "\n===========\nNumber of Voronoi half edges : " << vd.number_of_halfedges()  << std::endl;
 for(VD::Halfedge_iterator it = vd.halfedges_begin(); it!=vd.halfedges_end(); it++){
 //cout << *it << endl ;
 cout << it->is_ray() << endl ;
 cout << it->is_segment() << endl ;
 //  cout << it->
 
 }*/


void Voronoi::generateVoronoi(Configuration* _config, Circuit* _circuit, vector<geometry::Point*> p){
    forest_enabled=false;
    Site_2 t;
    for(int i=0;i<p.size();i++){
        vd.insert( Site_2(p[i]->getData(0), p[i]->getData(1)) );
    }
    assert( vd.is_valid() );
}

void Voronoi::generateVoronoi(Configuration* _config, Circuit* _circuit, Forest* forest, int index1, int index2){
    forest_enabled = true;
    for(int i=0;i<forest->getSize();i++){
        vd.insert( Site_2(forest->getNode(i)->getData(index1), forest->getNode(i)->getData(index2) ) );
    }
    assert( vd.is_valid() );
}

string Voronoi::drawDelaunay(){
    stringstream str ;
    for(std::list<Segment_2>::iterator it = vor.m_cropped_vd.begin(); it!=vor.m_cropped_vd.end(); it++){
        float x1 = it->source().x() ;
        float y1 = it->source().y() ;
        float x2 = it->target().x() ;
        float y2 = it->target().y() ;
        str << " set arrow from " << x1 << "," << y1 <<  "   to     " << x2 << "," << y2  << "      nohead  lc rgb \"gold\" lw 1 \n" ;
    }
            str  << " replot \n" ;
    return str.str();
    
}

string Voronoi::drawVoronoi(Forest* forest){
    stringstream str ;
    int id=1;
    std::cout << "\n===========\nNumber of Voronoi faces : " << vd.number_of_faces ()  << std::endl;
    for(VD::Face_iterator it = vd.faces_begin(); it!=vd.faces_end(); it++, id++){
        if(forest->getNode(id-1)->getTreeType() == forward ){
            if(it->is_unbounded()){
                /*
                 half-edge iterator:
                 VD::Face::Ccb_halfedge_circulator c_halfedge = it->outer_ccb();
                 VD::Face::Ccb_halfedge_circulator c_halfedge_first = c_halfedge ;
                 do{
                 c_halfedge++;
                 }while (c_halfedge != c_halfedge_first);
                 */
            }else{//Voronoi cell is bounded
                VD::Face::Ccb_halfedge_circulator c_halfedge = it->outer_ccb();
                VD::Face::Ccb_halfedge_circulator c_halfedge_first = c_halfedge ;
                str << "set object " << id << " polygon from " ;
                
                double x0 = c_halfedge->target()->point().x();
                double y0 = c_halfedge->target()->point().y();
                
                do{
                    double x1 = c_halfedge->target()->point().x();
                    double y1 = c_halfedge->target()->point().y();
                    double x2 = c_halfedge->source()->point().x();
                    double y2 = c_halfedge->source()->point().y();
                    //cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl ;
                    str << x1 << "," << y1 << "  to " ;
                   
                    c_halfedge++;
                    
                } while (c_halfedge != c_halfedge_first);
                str << x0 << "," << y0 << " "  << endl ;
                str << "set object " << id <<" fc rgb \"red\" fillstyle solid 1.0 border lt -1" << endl ;
            }
        }
    }
        str  << " replot \n" ;
    return str.str();
}
