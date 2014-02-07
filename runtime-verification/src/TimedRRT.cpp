#include "TimedRRT.h"

TimedRRT::TimedRRT(string fileName): RRT(fileName){
}

TimedRRT::TimedRRT(int _d, int _k, int _var, double _simTime, string nam): RRT(_d+1, _k, _var, nam){ //Add time dimension to the RRT
    sim_time = _simTime;
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

TimedRRT::TimedRRT(int _d, int _k, int _var, string nam): RRT(_d+1, _k, _var, nam){ //Add time dimension to the RRT
    sim_time = 1 ; //default
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

void TimedRRT::build(double* initialState){
	double timeEnvlope=dt;		//time_envlope is the latest sampled discovered so far

	//first, we construct the root from the given initial state
	root = new node(d);
    root->set(initialState);
	root->setRoot();
	nodes.push_back(root);
	root->setIndex(0);

	int i = 1;
	bool simulationFinished = false;
	while (!simulationFinished){
		cout << "#################################################  Iteration  " << i << endl;
		//create a new sample
		node* q_sample = new node(d);
		q_sample->randomize(min, max);
		//ensuring forward progress, initially we push the timeEnvlope to the simTime to gain depth, 
		//when we reached that time, we use it normally to ensure breath. 
		if (timeEnvlope >= sim_time){
			//double t = q_sample->unifRand(0, sim_time);
			q_sample->set(d - 1, -1);
		}
		else{
			q_sample->set(d - 1, timeEnvlope);
		}

		//find nearest node in the tree
		vector<node*> q_near_vec = getNearestNode(q_sample, -1, true);
		node* q_near = q_near_vec[0];
		double* state_near = q_near->get();
		double* ic = new double[d];
		for (int j = 0; j<d; j++){
			ic[j] = state_near[j];
		}

		vector<double> param;	//variation or input to the system 
		for (int j = 0; j < var; j++){
			if (j < 2){
				param.push_back(unifRand(variationMin[j], variationMax[j]));
			}
			else{
				double tnow = i*dt;
				double tperiod = 100e-12;	//100ps -> 10GHz
				double freq = 1 / tperiod;
				int cycles = ceilf(tnow / tperiod) - 1;
				double tp = tnow - cycles * tperiod;

				double vin = 0; int vinIndex = 2;
				double vdd = 1;
				double gnd = 0;
				if (tp < 0.1*tperiod){
					//pulse rising
					vin = gnd + vdd* ((tp / tperiod) / 0.1) + unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				else if (tp<0.5*tperiod){
					//pulse is 1
					vin = vdd - unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				else if (tp < 0.6*tperiod){
					//pulse falling
					double x = vdd* (((tp / tperiod) - 0.5) / 0.1);
					vin = vdd - x - unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				else{
					//pulse is 0
					vin = gnd + unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				param.push_back(vin);
			}

		}

		vector<string> settings;
		stringstream icInputFileName; icInputFileName << "ic_" << q_near->getIndex() << ".ic0";
		stringstream icOutputFileName; icOutputFileName << "ic_" << i << ".ic";
		settings.push_back(icOutputFileName.str());
		settings.push_back(icInputFileName.str());
		settings.push_back("transient"); // or settings.push_back("dc");


		double t_init = ic[d - 1];
		vector<double> result = system->simulate(ic, param, settings, 0, dt);


		result[0] += t_init;		//result[0] contains the time-stamp, in the simulation it is stamped as dt, however we have to add the time of the parrent node as well.

		if ((result[0]>timeEnvlope) && (timeEnvlope <= sim_time)){
			cout << "Time envlope pushed to " << timeEnvlope << endl;
			timeEnvlope = result[0];
		}

		cout << endl << "Next state is " << result[1] << " " << result[2] << "      / simEnvlope=" << timeEnvlope << endl;

		node* q_new = new node(d);
		q_new->set(result);
		cout << "New node is " << q_new->toString() << endl;
		cout << "Parrent Node is" << q_near->toString() << endl;
		q_near->addChildren(q_new);		//add the new node to the tree
		q_new->setParent(q_near);		//We only make the parent-child releation ship during the tree build
		q_new->setIndex(i);
		nodes.push_back(q_new);
		//casting
		//vector<node*> cast = getNearestNode(q_sample, 0.1, true);
		//cout << endl << "CAST size=" << cast.size() << endl ;
		//for(int j=0;j<cast.size();j++){
		//	cout << "Cast id = " << cast[j]->getID() << endl ;
		//	q_new->addCast(cast[j]);
		//	cast[j]->addCast(q_new);
		//}

		for (int j = 0; j<monitors.size(); j++){
			monitors[j]->check(q_new);
		}
		delete ic;

		i++;
		if (i >= k || timeEnvlope >= sim_time){
			simulationFinished = true;
		}

	}
}

void TimedRRT::simulate(double* initialState){
	root = new node(d);
	root->set(initialState);
	root->setRoot();
	nodes.push_back(root);
	root->setIndex(0);
	for(int i=1;i<k; i++){
    	cout <<"[s] #### " << i << endl ;
		node* q_near = nodes[ nodes.size()-1 ]; //last inserted node for the linear simulation
        double* state_near = q_near->get() ;
        double* state = new double[d];
        for(int j=0;j<d;j++){
			state[j]=state_near[j];
		}

		vector<double> param;	//variation or input to the system 
		for (int j = 0; j < var; j++){
			if (j < 2){
				param.push_back(unifRand(variationMin[j], variationMax[j]));
			}
			else{
				double tnow = i*dt;
				double tperiod = 100e-12;	//100ps -> 10GHz
				double freq = 1 / tperiod;
				int cycles = ceilf(tnow / tperiod) - 1;
				double tp = tnow - cycles * tperiod;
				
				double vin = 0; int vinIndex = 2;
				double vdd = 1;
				double gnd = 0;
				if (tp < 0.1*tperiod){
					//pulse rising
					vin = gnd+ vdd* ((tp/tperiod)/0.1) + unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				else if (tp<0.5*tperiod){
					//pulse is 1
					vin = vdd - unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				else if (tp < 0.6*tperiod){
					//pulse falling
					double x=vdd* ( ((tp / tperiod) - 0.5)/0.1); 
					vin = vdd -x- unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}else{
					//pulse is 0
					vin = gnd+ unifRand(variationMin[vinIndex], variationMax[vinIndex]);
				}
				param.push_back(vin);
			}
			
		}

		
		vector<string> settings;
		stringstream icInputFileName; icInputFileName << "ic_" << q_near->getIndex() << ".ic0";
		stringstream icOutputFileName; icOutputFileName << "ic_" << i << ".ic";
		settings.push_back(icOutputFileName.str());
		settings.push_back(icInputFileName.str());
		settings.push_back("transient");

		double t_init = state[d - 1];
		vector<double> result = system->simulate(state, param, settings, 0, dt);
		result[0] += t_init;		//result[0] contains the time-stamp, in the simulation it is stamped as dt, however we have to add the time of the parrent node as well.

		node* q_new = new node(d);

		for (int j = 0; j < d; j++){
			cout << j << "->" << result[j] << endl;
		}

        q_new->set(result);
        q_near->addChildren(q_new);		//add the new node to the tree
		q_new->setParent(q_near);		//We only make the parent-child releation ship during the tree build
		q_new->setIndex(i);
		nodes.push_back(q_new);

		for(int i=0;i<monitors.size();i++){
			monitors[i]->check(q_new);
		}	
    }
}

void TimedRRT::addMonitor(Monitor* m){
	monitors.push_back(m);
}

double TimedRRT::getSimTime(){
	return sim_time;
}