#include "circuit.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <sstream>
#include <math.h>
#include <typeinfo>

#include "point.h"
#include "polytope.h"
#include "geometryUtil.h"
#include "rootFinding.h"

using namespace std;
using namespace geometry;

Circuit::Circuit(SystemType type){
 if(type==vanderpol){
            dim  =   2;
            min = new double[dim];
            max = new double[dim];
            min[0] = -10 ;  //xmin
            min[1] = -10 ;  //ymin
            max[0] =  10 ;  //xmax
            max[1] =  10 ;  //ymax
            
            double* init_min = new double[dim];
            double* init_max = new double[dim];
            
            init_min[0] = -3 ;  //xmin
            init_max[0] = -2.8 ;  //xmax
            
            init_min[1] =  3 ;  //ymin
            init_max[1] =  3.1 ;  //ymax
            //TODO: should be replaced with a general object type
            
			//vector<Point*> v ;
           // Point* p0 = new Point(2);   p0->setData(0, init_min[0]); p0->setData(1, init_min[1]);   v.push_back(p0);
           // Point* p1 = new Point(2);   p1->setData(0, init_max[0]); p1->setData(1, init_min[1]);   v.push_back(p1);
          //  Point* p2 = new Point(2);   p2->setData(0, init_max[0]); p2->setData(1, init_max[1]);   v.push_back(p2);
          //  Point* p3 = new Point(2);   p3->setData(0, init_min[0]); p3->setData(1, init_max[1]);   v.push_back(p3);
          //  poly = new Polytope(2, v);
            
            mu = 1;
        }
}
Circuit::~Circuit(){}
  double Circuit::getMin(int i){
        if(i<dim) 
            return min[i];
        else
            return 0;
    }
    
    double Circuit::getMax(int i){
        if(i<dim)
            return max[i];
        else
            return 0;
    }
    
    int Circuit::getDimension(){
        return dim;
    }
    
    double* Circuit::getMinn(){
        return min;
    }
    
    double* Circuit::getMaxx(){
        return max ;
    }
    
	//Todo: implement this
    double* Circuit::getInitialState(){
        double* state = new double[dim];
        
		return state;
    }
    
    double Circuit::random(double a, double b){
        return (b-a)*(rand() / double(RAND_MAX)) + a;
    }
    
    double* Circuit::getVector(double* y){
        double* f = new double[dim];
        for(int i=0;i<dim;i++) f[i]=0;
        func(0, y, f, &mu);
        return f;
    }
    
    string Circuit::generateVectorField(){
        stringstream str;
        double coef = 20; 
        double dmax = 1 ;
        double* f = new double[dim];
        for(int i=0;i<dim;i++) f[i]=0;
        
        double* z = new double[dim];
        for(int i=0;i<dim;i++) z[i]=0;
        
        str << "plot '-' using 1:2:3:4 ti \"phase portrait\" with vectors lt 2 head filled " << endl ;
        int pointScale=2;
        for(int i=1; i<pointScale*max[0]; i++){
            for(int j=1; j<pointScale*max[1]; j++){
                double x = (double)(i+min[0]); ///(double)pointScale ;
                double y = (double)(j+min[1]); ///(double)pointScale ;
                z[0] = x;
                z[1] = y;
                func(0, z, f, &mu);
                double dx = f[0]/coef ;
                double dy = f[1]/coef ;
                
                if(dx>dmax) dx=dmax;
                if(dy>dmax) dy=dmax;
                if(dy<-dmax) dy=-dmax;
                if(dx<-dmax) dx=-dmax;
                
                str << x << " " << y << " " << dx << " " << dy << endl ;
            }
        }
        str << "e" << endl ;
        delete f;
        delete z;
        return str.str();
    }
    
    //------------------------------------------------------------
    //
    //    Transient Simulator for the System
    //
    //------------------------------------------------------------
    //Func will define the dynamics of the system. This is Van der Pol Oscillator
    int func (double t, const double y[], double f[], void *params){
        double mu = *(double *)params;
        f[0] = y[1];
        f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
        //f[0] = 2*y[0] -y[1] - y[0]*y[0]*y[0] ;  //x'=2x-y-x^3
        //f[1] = y[0] ;                           //y'=x
        return GSL_SUCCESS;
    }
    
    int jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
        double mu = *(double *)params;
        gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
        gsl_matrix * m = &dfdy_mat.matrix; 
        gsl_matrix_set (m, 0, 0, 0.0);
        gsl_matrix_set (m, 0, 1, 1.0);
        gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
        gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
        dfdt[0] = 0.0;
        dfdt[1] = 0.0;
        return GSL_SUCCESS;
    }
    
    vector<pair<double, double> > Circuit::simulate(){
        vector<pair<double, double> > trace ;
        gsl_odeiv2_system sys = {func, jac, 2, &mu};
        gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
        double t = 0.0, t1 = 20.0;
        double* state = getInitialState();
        double y[2] = { state[0], state[1] };
        for (int i = 1; i <= 1000; i++){
            double ti = i * t1 / 1000.0;
            int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
            if (status != GSL_SUCCESS){
                printf ("error, return value=%d\n", status);
                break;
            }
            trace.push_back(make_pair(y[0], y[1]));
            //printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
        }
        delete state;
        gsl_odeiv2_driver_free (d);
        return trace;
    }    
    
	/*
    //------------------------------------------------------------
    //
    //    Part II : Reachability Decision for Polytopes
    //
    //------------------------------------------------------------
    
    ///
    ///     Algorithm for decising wether or not the sink node is reachable from the source node
    ///     based on the direction of system trajectoreis at the borders of source and sink.
    ///
    ///     1. Find the common edge between source and sink
    ///     2. Find two points P1 and P2 such as the line p1-p2 is the biggest common segment between source and sink
    ///     3. Figure-out the left and right. 
    ///         a. We decide the reachability from the left of line P1-P2 to the right of the line P1-P2 in term of looking from P1 to P2
    ///             or the otherway.
    ///         b. We need to determine whether the source is at the left or at the right of P1-P2 line
    ///         c. Check the location of the centeroid of the polytope versus the P1-P2 line.
    ///         d. check the sign of the determint of the following matrix:
    ///                 | x_2 - x_1     x - x_1 |
    ///                 | y_2 - y_1     y - y_1 |  
    ///             So position = sign( (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) ), It is 0 on the line, and +1 on one side, -1 on the other side.
    ///     4. Asuming that source is at the left, if not change the sign of reachability, Decide on the reachability of the left
    ///         of P1_P2 to the right of P1-P2. For every dimension, Rotate the p1-p2 to the axis and then determine the reachability.
    ///         if sink is reachable for every dimension, then the sink is reachable.
    ///    5. Compute the slope and theta of line p1-p2
    ///    6. x-dim: rotate the system for (360-theta), check on reachablity
    ///    7. y-dim: rotate the system for (90-theta), check for reachablity
    ///    8. return the result
    bool Circuit::isReachable(Polytope* source, Polytope* sink){
    	if(source->id==7 && sink->id==10) return false;
        //Step 1: Find the adjacent edge.
        vector<Point*> interval = source->getSharedPoints(source, sink);
        Point* sourceCenter = source->getCentroid();
        Point* sinkCenter   = sink->getCentroid();
        int position1 = geometry::position(interval[0], interval[1], sourceCenter) ;
        int position2 = geometry::position(interval[0], interval[1], sinkCenter) ;
        functionSign sign = reachabilityUsingSampling(interval[0], interval[1]);
        //This is a buggy code.
        if( (sign==posnegative) ||
            ((sign==positive)&&(position2==+1))||
            ((sign==negative)&&(position2==-1))
           ){
            return true ;
        }else{
            return false;
        }
        return true;
    }
    
    functionSign Circuit::reachabilityUsingSampling(Point* p1, Point* p2){
        functionSign result = undefined; 
        int samples=100;
        double* f = new double[2];
        
        Point* tmp = new Point(2);
        double p1x = p1->getData(0);
        double p1y = p1->getData(1);
        double p2x = p2->getData(0);
        double p2y = p2->getData(1);
        for(int i=1;i<samples+1;i++){
            //sample point
            double y[2] = { p1x + i*(p2x-p1x)/(samples+1), p1y + i*(p2y-p1y)/(samples+1) };
            //code 1: sampling using differential trajectory:
            //func(0, y, f, &mu);
            //y[0] += f[0];
            //y[1] += f[1];
            //code 2: sampling using simulation
            gsl_odeiv2_system sys = {func, jac, 2, &mu};
            gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
            double t = 0.0, ti = 1;
            gsl_odeiv2_driver_apply (d, &t, ti, y);
            gsl_odeiv2_driver_free (d);
            tmp->setData(0, y[0]);
            tmp->setData(1, y[1]);
            
            //positition of the trejectory w.r.t p1-p2 line
            int position1 = geometry::position(p1, p2, tmp) ;

            //compute the direction of vector flow at the sampling point
            switch (result) {
                case positive:
                    result = (position1==1)? positive: posnegative;
                    break;
                case negative:
                    result = (position1==1)? posnegative: negative;
                    break;
                case posnegative:
                    result = (position1==1)? posnegative: posnegative;
                    break;
                case undefined:
                    result = (position1==1)? positive: negative;
                    break;
            }
            
        }
        delete f;
        return result ;
    }
    
    
    
    //axis 0 is x-axis, 1 is y, 2 ...
    functionSign Circuit::reachability(Point* p1, Point* p2, int axis){
        double p1x = p1->getData(0);
        double p1y = p1->getData(1);
        double p2x = p2->getData(0);
        double p2y = p2->getData(1);
        
        //Computing transformation angle theta (in radians)
        double dx  = p2x - p1x ;
        double dy  = p2y - p1y ;
        double theta = atan2( dy , dx ) ;     
        theta = (axis==0? -theta: -theta+ (M_PI/2) ) ; //This line took 2 days to develop!
        double sint = sin(theta) ;
        double cost = cos(theta) ;
        double p3x = (p1x+p2x)/2 ;
        double p3y = (p1y+p2y)/2 ;
        pair<double, double> p1_rotated = geometry::rotate(p1x, p1y, theta);
        pair<double, double> p2_rotated = geometry::rotate(p2x, p2y, theta);
        pair<double, double> p3_rotated = geometry::rotate(p3x, p3y, theta);
        
        if(axis==0){
            //case 1  --- ^
            // y is fixed, x is varied, check for f_y
            // y' = f(x)sint + f(y)cost
            // y' = ysint + ycost - yx^2cost - xcost
            
            double y = p1_rotated.second ;
            
            double a_0 = y*sint + y*cost ;
            double a_1 = -cost ;
            double a_2 = -y*cost ;
            
            double min = p1_rotated.first ;
            double max = p2_rotated.first ;
            vector<double> coeff;
            coeff.push_back(a_0);
            coeff.push_back(a_1);
            coeff.push_back(a_2);
            //Solution solution = findRoot(coeff, min, max);
            Solution solution = findRoot_fdf(coeff, min, max);            
            if(solution.existance){
                return posnegative ;
            }else{
                double x = p3_rotated.first ;
                double y = polynomial ( p3_rotated.first, coeff);
                if(y>0) return positive;
                else return negative;
            }
        }else{  
            // case 2  
            //  x'= xsint + ycost (wrong)
            //  x'=f(x)sin(t) + f(y)cos(t)
            //  x'=   ycos(t) - ((1-x^2)y-x)sin(t)
            //  x'= ycost - ysint + yx^2sint + xsint 
            //  | -> x is fixed, y is varied, check for f_x
            double x = p1_rotated.first ;
            
            double a_0 = x*sint ;
            double a_1 = cost - sint + x*x*sint  ;
            
            double min = p1_rotated.second ;
            double max = p2_rotated.second ;
            
            vector<double> coeff;
            coeff.push_back(a_0);
            coeff.push_back(a_1);
            //Solution solution = findRoot(coeff, min, max);
            Solution solution = findRoot_fdf(coeff, min, max);            
            if(solution.existance){
                return posnegative ;
            }else{
                double y = p3_rotated.second ;
                double x = polynomial ( y, coeff);
                if(x>0) return positive;
                else return negative;
            }
        }
    }*/