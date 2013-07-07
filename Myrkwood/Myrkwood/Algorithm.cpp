#include <iostream>
#include "Algorithm.h"
#include "sampling.h"
#include "Forest.h"

using namespace std;

Algorithm::Algorithm(Circuit* _circuit, Configuration* _config){
    circuit = _circuit;
    config = _config;
}

Algorithm::~Algorithm(){}

//Todo: implement this method
string Algorithm::compute_reachable_set(Trace* trace){
    return generate_counter_example(trace);
}

string Algorithm::generate_counter_example(Trace* trace){
    int k = 1;                  config->getParameter("edu.uiuc.crhc.core.sampling.number", &k);
    double temperature=0.5;     config->getParameter("edu.uiuc.crhc.core.sampling.temperature", &temperature);

    RPG* rpg = new RPG(circuit);
    rrf = new Forest(config, circuit);
    double* unsafe_point_data = new double[circuit->getDimension()];
    for(int i=0;i<circuit->getDimension();i++) unsafe_point_data[i]=0;
    
    unsafe_point_data[0]=    0;
    unsafe_point_data[1]=    5;
    unsafe_point_data[2]=    3;
    unsafe_point_data[3]=    3;
    unsafe_point_data[4]=    6;
    unsafe_point_data[5]=    3;
    unsafe_point_data[6]=    3;
    unsafe_point_data[7]=    0;
    geometry::Point* unsafe_point = new geometry::Point(circuit->getDimension(), unsafe_point_data );
    rrf->addRoot(unsafe_point, backward);
    
    cout << trace->getSampleSize() << endl ;
    double tmax = 1e-3;         config->getParameter("edu.uiuc.crhc.core.simulation.totalTime", &tmax);
    double time_resoulution = tmax / trace->getSampleSize();
    for(int i=1;i<trace->getSampleSize();i++){
        geometry::Point* q_new = new geometry::Point( trace->getSample(i));
        Node* node = rrf->addNode(q_new, rrf->getNode(i-1), forward);
        node->setTime(i*time_resoulution);
        //cout << "timing info: " << i * time_resoulution <<  "       parent->" << rrf->getNode(i-1)->getTime() << " " << rrf->getNode(i-1)->getTime()+time_resoulution << endl ;
    }
    
    for(int i=0;i<k; i++){
        cout << i << endl ;
        geometry::Point* q_sample = rpg->sample();
        Node* q_near   = rrf->find_nearest_point(q_sample) ;
        pair<geometry::Point*, Type> q_new    = rrf->find_optimum_trajectory(q_sample, q_near);
        if(q_new.second!=invalid){
            if(random(0,1)<=0.9){
                rrf->addNode(q_new.first, q_near, q_new.second);
            }else{
                rrf->addRoot(q_sample, coupling);
            }
        }
        delete q_sample ;
    }
    return rrf->draw(1,7) ;
}

string Algorithm::forward_exploration(Trace* trace){
    int k = 1;                  config->getParameter("edu.uiuc.crhc.core.sampling.number", &k);
    RPG* rpg = new RPG(circuit);
    rrf = new Forest(config, circuit);
    
    
    //Initializing an RRT with a single trace
    cout << trace->getSampleSize() << endl ;
    double tmax = 1e-3;         config->getParameter("edu.uiuc.crhc.core.simulation.totalTime", &tmax);
    double time_resoulution = tmax / trace->getSampleSize();
    for(int i=1;i<trace->getSampleSize();i++){
        geometry::Point* q_new = new geometry::Point( trace->getSample(i));
        Node* node = rrf->addNode(q_new, rrf->getNode(i-1), forward);
        node->setTime(i*time_resoulution);
        //cout << "timing info: " << i * time_resoulution <<  "       parent->" << rrf->getNode(i-1)->getTime() << " " << rrf->getNode(i-1)->getTime()+time_resoulution << endl ;
    }
    
    for(int i=0;i<k; i++){
        cout << i << endl ;
        geometry::Point* q_sample = rpg->sample();
        Node* q_near   = rrf->find_nearest_point(q_sample) ;
        pair<geometry::Point*, Type> q_new    = rrf->find_optimum_trajectory(q_sample, q_near);
        rrf->addNode(q_new.first, q_near, q_new.second);
        delete q_sample ;
    }
    return rrf->draw() ;
    
}

Forest* Algorithm::getForest(){
    return rrf;
}
