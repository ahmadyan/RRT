//  Venice is a code for determining the area of a circle using a MC method vs randomized reachability method
//  Venice.cpp
//  Myrkwood
//
//  Created by Adel Ahmadyan on 3/25/13.
//  Copyright (c) 2013 University of Illinois at Urbana-Champaign. All rights reserved.
//

/*  Results:
    Exact Area:  0.785398
    
    MC with 95% confidence interval
    Samples         Estimated Area          lower_bound         upper_bound
    1000            0.779                   0.739076            0.818924
 */

#include "Venice.h"
#include <math.h>
#include <sstream>
#include <boost/math/distributions/students_t.hpp>


stringstream mc_stream;
stringstream verify_stream;

int sample_count = 1000;
double xmin = 0 ;
double xmax = 1 ;
double ymin = 0 ;
double ymax = 1 ;


//returns true of x,y is inside circle of radius r
bool check(double x, double y, double r){
    if( (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) < r*r ){
        return true;
    }else{
        return false;
    }
}

void venice_demo(){
    double r = 0.5 ;
    double area = M_PI * r * r;

    double** samples = new double*[sample_count];
    for(int i=0;i<sample_count;i++){
        samples[i] = new double[2];
        samples[i][0] = 1.0*rand()/INT_MAX;
        samples[i][1] = 1.0*rand()/INT_MAX;
    }
    cout << "The following demo computes the area of a unit circle using " << endl;
    cout << "monte carlo sampling methods with given confidence interval and the new randomized verification" << endl;
    cout << "method" << endl;
    cout << "The exact area is " << area << endl ;
    venice_mc(r, samples);
    venice_randomized_verify(r, samples);
}

void venice_mc(double r, double** samples){
    cout << "Monte Carlo simulation" << endl ;
    using namespace boost::math;
    
    double sample_inside=0;
    cout << "sizeof(sample_count)=" << sizeof(samples) << endl ;
    for(int i=0;i<sample_count;i++){
        double x = samples[i][0];
        double y = samples[i][1];
        
        if(check(x,y,r)){
            sample_inside++;
            mc_stream << "set obj "<<  10+i<< " circle at " << x << " ," << y << "  size .001 fc rgb \"blue\" \n" << endl ;
        }else{
            mc_stream << "set obj "<<  10+i<< " circle at " << x << " ," << y << "  size .001 fc rgb \"red\" \n" << endl ;
        }
        
    }
    cout << " samples=" << sample_inside << " " << sample_count << endl ;
    //computing the confidence interval
    double area = (double)sample_inside / (double)sample_count ;
    double deviation = sqrt( area * (1-area));
    double alpha = 0.05;
    //students_t dist(deviation-1);
    //double T = quantile( complement( dist, alpha/ 2));
    double T = 1.96; //corresponds with 95% confidence interval
    cout << "T=" << T << endl ;
    double w = T * sqrt(deviation / (double)sample_count) ;
    //double w = T*deviation / sqrt(double(deviation));
    cout << "w=" << w << endl ;
    cout << "confidence interval= [" << area - w << " , " << area + w << " ] " << endl ;
    
    cout << "estimated area = " << area << endl ;
}

void venice_randomized_verify(double r, double** samples){
    DT dt2;
    VD vd;
    Cropped_voronoi_from_delaunay vor;
    
    Iso_rectangle_2 bbox(0,0,1,1);
    vor = Cropped_voronoi_from_delaunay(bbox);
    for(int i=0;i< sample_count;i++){
        vd.insert( Site_2( samples[i][0], samples[i][1] ) );
    }
    assert( vd.is_valid() );
    
    double area =0;
    double area2=0;
    double area3=0;
    int id=100;
    std::cout << "\n===========\nNumber of Voronoi faces : " << vd.number_of_faces ()  << std::endl;
    for(VD::Face_iterator it = vd.faces_begin(); it!=vd.faces_end(); it++, id++){
        vector<pair<double, double> > points ;
        if(it->is_unbounded()){
            cout << "hellp!" << endl ;
            VD::Face::Ccb_halfedge_circulator c_halfedge = it->outer_ccb();
            VD::Face::Ccb_halfedge_circulator c_halfedge_first = c_halfedge ;
            do{
				Point_2 raybase, raypoint; //base the origin, point is where the ray point to
				if(c_halfedge->is_unbounded()){
					if(c_halfedge->has_source()){
						Point_2 temp( c_halfedge->source()->point().x(), c_halfedge->source()->point().y()) ;
						raybase = temp;
						std::cout<<"Source : "<< raybase << std::endl;
					}
					if(c_halfedge->has_target()){
						Point_2 temp( c_halfedge->target()->point().x(), c_halfedge->target()->point().y()) ;
						raybase = temp;
						std::cout<<"Target : "<< raybase << std::endl;
						//calculate the mid point of two nearby site
						Point_2
                        midpoint(
                                 (c_halfedge->right()->point().x() + c_halfedge->opposite()->left()->point().x())/2.0,
                                 (c_halfedge->right()->point().y() + c_halfedge->opposite()->left()->point().y())/2.0
                                 );
						std::cout<<"Midpoint: " <<c_halfedge->up()->point()<<std::endl;
					}
				}
				c_halfedge++;
            }while (c_halfedge != c_halfedge_first);
            
        }else{//Voronoi cell is bounded
            
            bool polyIsInside = true;
            bool polyIsOutside = true;
            VD::Face::Ccb_halfedge_circulator c_halfedge = it->outer_ccb();
            VD::Face::Ccb_halfedge_circulator c_halfedge_first = c_halfedge ;
            verify_stream << "set object " << id << " polygon from " ;
            
            double x0 = c_halfedge->target()->point().x();
            double y0 = c_halfedge->target()->point().y();
            points.push_back(make_pair(x0, y0));
            if(!check(x0, y0, r)) polyIsInside = false;
            if(check(x0, y0, r)) polyIsOutside = false ;
            do{
                
                double x1 = c_halfedge->target()->point().x();
                double y1 = c_halfedge->target()->point().y();
                double x2 = c_halfedge->source()->point().x();
                double y2 = c_halfedge->source()->point().y();
                if(!check(x1, y1, r)) polyIsInside = false;
                if(!check(x2, y2, r)) polyIsInside = false;
                if(check(x1, y1, r)) polyIsOutside = false ;
                if(check(x2, y2, r)) polyIsOutside = false ;
                points.push_back(make_pair(x1, y1));
                verify_stream << x1 << "," << y1 << "  to " ;
                c_halfedge++;
                
            } while (c_halfedge != c_halfedge_first);
            verify_stream << x0 << "," << y0 << " "  << endl ;
            if( polyIsInside ){
                verify_stream << "set object " << id <<" fc rgb \"blue\" fillstyle solid border lt 0.5" << endl ;
            }else if (polyIsOutside){
                verify_stream << "set object " << id <<" fc rgb \"red\" fillstyle solid border -1" << endl ;
            }else{
                verify_stream << "set object " << id <<" fc rgb \"white\" fillstyle solid border lt 0.5" << endl ;
            }
            
            //Compute the area of the Voronoi cell, this is only valid for two-dimensions
            //In n-dimension, computing the volume of polytope is P-Complete
            double sum0=0, sum1=0, sum2=0;
            for(int i=0;i<points.size()-1;i++){
                sum0 += points[i].first * points[i+1].second;
            }
            sum0 += points[points.size()-1].first * points[0].second;
            
            for(int i=0;i<points.size()-1;i++){
                sum1 += points[i].second * points[i+1].first;
            }
            sum1 += points[points.size()-1].second * points[0].first;
            sum2=sum0-sum1;
            sum2/=2;
            double volume=abs(sum2);
            
            if ( polyIsInside )
                area += volume ;
            else if (polyIsOutside)
                area2 += volume;
            else
                area3 += volume ;
        }
    }
    cout << "area (lower bound) = " << area << endl ;
    cout << "area (upper bound) = " << 1-area2 << endl ;
    cout << "area (uncertain) = " << area3 << endl ;
    verify_stream  << " replot \n" ;
}

string venice_draw(){
    bool drawMC=false;
    stringstream str ;
    //str << "set style fill transparent solid  0.5" << endl ;
    str << "unset border ; unset tics ; unset ztics " << endl ;
    str << "set format y \"\"" << endl ;
    str << "set format x \"\"" << endl ;
    str << "set size square \n " << endl ;
   
    if(drawMC){
        str << mc_stream.str() << endl ;
    }else{
        str << verify_stream.str() << endl ;
    }

    str << "set object 1 rect from 0,0 to 1,1 fc rgb \"red\" fillstyle solid border -1 \n" ;
    str << "set obj 9999999 circle at  .5,.5  size .50 fc rgb \"black\" fs transparent border 2\n" ;
    
    str  << "replot \n" ;
    return str.str();
}