#include "TimedRRT.h"

TimedRRT::TimedRRT(string fileName): RRT(fileName){
}

TimedRRT::TimedRRT(int _d, int _k, double _simTime, string nam): RRT(_d+1, _k, nam){ //Add time dimension to the RRT
    sim_time = _simTime;
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

TimedRRT::TimedRRT(int _d, int _k, string nam): RRT(_d+1, _k, nam){ //Add time dimension to the RRT
    sim_time = 1 ; //default
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

void TimedRRT::build(double* initialState, double variation){
	root = new node(d);
    root->set(initialState);
	root->setRoot();
    
    for(int i=0;i<k; i++){
    	cout <<"#### " << i << endl ;
        //create a new sample
        node* q_sample = new node(d);
        q_sample->timed_randomize(i, k, min, max);
        //q_sample.randomize(min, max);
        
        //find nearest node in the tree
        node* q_near = getNearestNode(q_sample);

        double* state_near = q_near->get() ;
        double* state = new double[d];
        for(int i=0;i<d;i++) state[i]=state_near[i];
        
        //variation or input to the system
        //todo: should be defined in the main or system, not here
		double* param = new double[1];
        param[0] = unifRand(-variation, variation);
        double t_init = state[d-1];
        double t = system->simulate(state, param, dt);
        state[d-1] = t_init + t ;
        cout << endl <<  "TIME " << t_init << " " << t << " " << state[d-1] << endl ;
        node* q_new = new node(d);
        q_new->set(state);
        //add the new node to the tree
        q_near->addChildren( q_new ) ;
		q_new->setParent(q_near);		//We only make the parent-child releation ship during the tree build
		//for(int i=0;i<monitors.size();i++)
		//	monitors[i]->verify( q_new ) ;
    }
}

void TimedRRT::addMonitor(Monitor* m){
	monitors.push_back(m);
}

node* TimedRRT::getNearestNode(node* q_sample){
	return root->getNearestNode(q_sample, min, max, true).first ;
}

double TimedRRT::getSimTime(){
	return sim_time;
}