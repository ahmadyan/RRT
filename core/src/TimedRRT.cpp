#include "TimedRRT.h"
#include <stdio.h>
#include <chrono>
#include <ctime>
#include <random>

TimedRRT::TimedRRT(Configuration* c, string fileName): RRT(c, fileName){
}

TimedRRT::TimedRRT(Configuration* c, int _d, int _k, int _var, double _simTime, string nam) : RRT(c, _d + 1, _k, _var, nam){ //Add time dimension to the RRT
    sim_time = _simTime;
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

TimedRRT::TimedRRT(Configuration* c, int _d, int _k, int _var, string nam) : RRT(c, _d + 1, _k, _var, nam){ //Add time dimension to the RRT
    sim_time = 1 ; //default
    setBound(_d, 0, sim_time); //maximum time for sim, minimum is 0 and it is assigned in RRT constructor
}

// We will only use this function for MC simulation of Inverter for now. 
// 1st: we determine the actual number of bits that should be applied in this simulation
// 2nd: we generate the boot sequence for getting initial eye diagram
// 3rd: we generate a random bit sequence
// 4th: we generate jitter and rise/fall time data
void TimedRRT::generateMonteCarloInputSequence(){
	double freq = 0; config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
	double period = 1 / freq;
	int IterationsPerBit = period / dt;
	int size = k / IterationsPerBit;

	bits = new int[size];
	jitter = new double[size];
	transition = new double[size];

	double jitterMean = 0; config->getParameter("edu.uiuc.csl.system.param.jitter.mean", &jitterMean);
	double jitterStdDev = 0; config->getParameter("edu.uiuc.csl.system.param.jitter.stddev", &jitterStdDev);
	double jitterMax = 0; config->getParameter("edu.uiuc.csl.system.param.jitter.max", &jitterMax);
	double transitionMean = 0; config->getParameter("edu.uiuc.csl.system.param.transition.mean", &transitionMean);
	double transitionStdDev = 0; config->getParameter("edu.uiuc.csl.system.param.transition.stddev", &transitionStdDev);
	double transitionMax = 0; config->getParameter("edu.uiuc.csl.system.param.transition.max", &transitionMax);

	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed1);
	std::normal_distribution<double> jitter_generator(jitterMean, jitterStdDev);
	std::normal_distribution<double> transition_generator(transitionMean, transitionStdDev);

	bits[0] = 0;		jitter[0] = 0;		transition[0] = 0;
	bits[1] = 0;		jitter[1] = 0;		transition[1] = 0;
	bits[2] = 1;		jitter[2] = 0;		transition[2] = 0;
	bits[3] = 1;		jitter[3] = 0;		transition[3] = 0;
	bits[4] = 0;		jitter[4] = 0;		transition[4] = 0;

	for (int i = 5; i < size; i++){
		bits[i] = rand() % 2;

		double x = jitter_generator(generator);
		if (x<0) x = 0;
		if (x>jitterMax) x = jitterMax;
		jitter[i] = x*1e-12;

		double y = transition_generator(generator);
		if (y<0) y = 0;
		if (y>transitionMax) y = transitionMax;
		transition[i] = y*1e-12;
	}

	ofstream file;
	file.open(config->get("edu.uiuc.crhc.core.options.mc.simdata"), std::ofstream::out);

	for (int i = 0; i < size; i++){
		cout << i << " " << bits[i] << " " << jitter[i] << " " << transition[i] << endl;
		file << i << " " << bits[i] << " " << jitter[i] << " " << transition[i] << endl;
	}

	file.close();
}

vector<double> TimedRRT::generateSimulationParameters(node* q_near){
	vector<double> param;	//variation or input to the system 
	for (int j = 0; j < var; j++){
		if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "sin")){
			double tnow = q_near->getTime() + dt;
			double Fc = 60;	config->getParameter("edu.uiuc.csl.system.param.freq", j, &Fc);	//Freq in Hertz
			double max = 0; config->getParameter("edu.uiuc.csl.system.param.max", j, &max);
			double min = 0; config->getParameter("edu.uiuc.csl.system.param.min", j, &min);
			double scale = max - min;
			double v = scale*sin(2 * 3.14159265358979323846 *tnow*Fc);	//Base signal

			double noise = 0;
			double dv = 0; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);	//noise
			if (config->checkParameter("edu.uiuc.csl.system.param.dist.type", j, "normal")){
				double mean = 0;  config->getParameter("edu.uiuc.csl.system.param.dist.mean", j, &mean);
				double var = 0;  config->getParameter("edu.uiuc.csl.system.param.dist.var", j, &var);
				noise = 0;
			}else{//uniform
				noise = generateUniformSample(-dv, +dv);
			}

			cout << "------------" << endl;
			cout << "TNow=" << tnow << endl;
			cout << "TNow=" << Fc << endl;
			cout << "Scale=" << scale << endl;
			cout << "V=" << v << endl;
			cout << "DV=" << dv << endl;
			cout << "noise=" << noise << endl;
			cout << "------------" << endl;
			cout << endl << endl;
			v += noise;
			param.push_back(v);
		}else if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "digital")){
			double dv = 0; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);
			double noise = generateUniformSample(-dv, dv);

			if (GammaSimMode == 1){
				double max = 0; config->getParameter("edu.uiuc.csl.system.param.max", j, &max);
				double min = 0; config->getParameter("edu.uiuc.csl.system.param.min", j, &min);
				double freq = 0; config->getParameter("edu.uiuc.csl.system.param.freq", j, &freq);
				double period = 1 / freq;
				double t = q_near->getTime();
				int cycles = floorf(t / period);
				double t1 = t - period*cycles;

				int bit = -1;	//the next bit
				double v = 0;
				//modeling jitter
				if ((t1 >= GammaJitter) || (cycles == 0)){
					bit = bits[cycles];
				}
				else{
					bit = bits[cycles - 1];
				}

				//modeling transition
				//double tTran = jitter[cycles] 
				if ((cycles != 0) && (t1 >= GammaJitter) && (t1 < GammaJitter + GammaTransition)){
					double tTran = t1 - GammaJitter;
					double transitionScale = tTran / GammaTransition;
					v = min + max*(transitionScale*bits[cycles] + (1 - transitionScale)*bits[cycles - 1]);
				}
				else{
					v = min + bit*max;
				}
				double dv = 0; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);
				double noise = generateUniformSample(-dv, dv);

				param.push_back(v + noise);


				cout << t1 << " " << bit << " " << GammaJitter << " " << GammaTransition << " " <<  v << "     ";
			}else{
				double currentInputValue = q_near->getInput(j);
				int currentDigit;
				int nextDigit = 0;
				if (currentInputValue > 0.5)
					currentDigit = 1;
				else
					currentDigit = 0;

				int counter = q_near->getCounter();
				int minPeriod = 50;
				int maxPeriod = 100;
				double switchProbability = 0.1;
				if (counter < minPeriod){
					//do not switch, it is too quick to switch now
					nextDigit = currentDigit;
				}
				else if (counter>maxPeriod){
					//definitely switch now
					nextDigit = 1 - currentDigit;
				}
				else{
					//switch with some probability
					double p = generateUniformSample(0, 1);
					if (p < switchProbability){
						nextDigit = 1 - currentDigit;
					}
					else{
						nextDigit = currentDigit;
					}
				}

				param.push_back(nextDigit + noise);
			}
		}
		else if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "pulse")){
			//cout << "Generating pulse!" << endl;
			double dv = 0; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);
			double noise = generateUniformSample(-dv, dv);
			double tnow = q_near->getTime() + dt;
			double tperiod = 100e-12;	//100ps -> 10GHz
			double freq = 1 / tperiod;
			int cycles = ceilf(tnow / tperiod) - 1;
			double tp = tnow - cycles * tperiod;

			double vin = 0; int vinIndex = 2;
			double vdd = 1;
			double gnd = 0;
			if (tp < 0.1*tperiod){
				cout << 1 << endl;
				//pulse rising
				vin = gnd + vdd* ((tp / tperiod) / 0.1) + noise;
			}
			else if (tp<0.5*tperiod){
				cout << 2 << endl;
				//pulse is 1
				vin = vdd - noise;
			}
			else if (tp < 0.6*tperiod){
				cout << 3 << endl;
				//pulse falling
				double x = vdd* (((tp / tperiod) - 0.5) / 0.1);
				vin = vdd - x - noise;
			}
			else{
				cout << 4 << endl;
				//pulse is 0
				vin = gnd + noise;
			}

			//cout << "tnow=" << tnow << endl;
			//cout << "v->" << vin << endl;
			param.push_back(vin);
		}else if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "mc")){		
			double max = 0; config->getParameter("edu.uiuc.csl.system.param.max", j, &max);
			double min = 0; config->getParameter("edu.uiuc.csl.system.param.min", j, &min);
			double freq = 0; config->getParameter("edu.uiuc.csl.system.param.freq", j, &freq);
			double period = 1 / freq;
			double t = q_near->getTime();
			int cycles = floorf(t / period);
			double t1 = t - period*cycles;

			int bit = -1;	//the next bit
			double v=0;
			//modeling jitter
			if ( (t1 >= jitter[cycles]) || (cycles==0)){
				bit = bits[cycles];
			}else{
				bit = bits[cycles - 1];
			}

			//modeling transition
			//double tTran = jitter[cycles] 
			if ((cycles!=0)&&(t1 > jitter[cycles]) && (t1 < jitter[cycles] + transition[cycles])){
				double tTran = t1 - jitter[cycles];
				double transitionScale = tTran / transition[cycles];
				v = min + max*(transitionScale*bits[cycles] + (1 - transitionScale)*bits[cycles - 1]);
			}
			else{
				v = min + bit*max;
			}
			double dv = 0; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);
			double noise = generateUniformSample(-dv, dv);

			param.push_back(bit + noise);

		}else if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "dc")){
			if (config->checkParameter("edu.uiuc.csl.system.param.dist.type", j, "gaussian") || config->checkParameter("edu.uiuc.csl.system.param.dist.type", j, "normal")){
				double mean = 0; config->getParameter("edu.uiuc.csl.system.param.dist.mean", j, &mean);
				double var = 1; config->getParameter("edu.uiuc.csl.system.param.dist.variance", j, &var);
				double std = sqrt(var);
				double v = generateTruncatedNormalSample(mean, std, variationMin[j], variationMax[j]);
				param.push_back(v); 
			}
			else{
				param.push_back(generateUniformSample(variationMin[j], variationMax[j]));
			}
		}else if (config->checkParameter("edu.uiuc.csl.system.param.type", j, "boot")){
			//We are going to generate the 00110 signal
			int boot[5] = { 0, 0, 1, 1, 0 };
			double tnow = q_near->getTime() + dt;
			double freq; config->getParameter("edu.uiuc.csl.system.param.freq", j, &freq);
			double period = 1 / freq;
			int cycles = ceilf(tnow / period) - 1 ;
			cycles = cycles % 5;
			double vin = boot[cycles] * 0.9;

			double dv; config->getParameter("edu.uiuc.csl.system.param.dv", j, &dv);
			double noise=0;
			if (config->checkParameter("edu.uiuc.csl.system.param.dist.type", j, "gaussian") || config->checkParameter("edu.uiuc.csl.system.param.dist.type", j, "normal")){
				double mean = 0; config->getParameter("edu.uiuc.csl.system.param.dist.mean", j, &mean);
				double var = 1; config->getParameter("edu.uiuc.csl.system.param.dist.variance", j, &var);
				double std = sqrt(var);
				noise=generateTruncatedNormalSample(mean, std, variationMin[j], variationMax[j]);
			}else{
				noise = generateUniformSample(-dv, dv);

			}
			vin += noise;
			param.push_back(vin);
		}else{
			param.push_back(generateUniformSample(variationMin[j], variationMax[j]));
		}
	}

	return param;
}

void TimedRRT::build(double* initialState){
	//first, we construct the root from the given initial state
	root = new node(d);
	root->set(initialState);
	root->setRoot();
	nodes.push_back(root);
	root->setIndex(0);
	root->setCounter(0);
	vector<double> p;
	for (int i = 0; i < var; i++)
		p.push_back(0);
	root->setInputVector(p);
	build();
}


void TimedRRT::build(){
	double timeEnvlope = dt;		//time_envlope is the latest sampled discovered so far
	timeEnvlope = 501e-12;
	bool forceIteration = false;
	if (config->checkParameter("edu.uiuc.csl.core.sampling.iteration.force", "1")){
		forceIteration = true;
	}
	config->getParameter("edu.uiuc.csl.core.sampling.iteration", &k);
	int i = nodes.size();
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
		vector<node*> q_near_vec = getNearestNode(q_sample);
		node* q_near = q_near_vec[0];
		
		double* state_near = q_near->get();
		double* ic = new double[d];
		for (int j = 0; j<d; j++){
			ic[j] = state_near[j];
		}


		vector<string> settings;
		stringstream icInputFileName; icInputFileName << "ic_" << q_near->getIndex() << ".ic0";
		stringstream icOutputFileName; icOutputFileName << "ic_" << i << ".ic";
		settings.push_back(icOutputFileName.str());
		settings.push_back(icInputFileName.str());
		settings.push_back("transient"); // or settings.push_back("dc");

		vector<double> param = generateSimulationParameters(q_near);
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
		q_new->setInputVector(param);


		

		/*
		int q_near_digit4;
		int q_new_digit4;

		if (q_near->getInput(4) > 0.5)
			q_near_digit4 = 1;
		else
			q_near_digit4 = 0;

		if (q_new->getInput(4) > 0.5)
			q_new_digit4 = 1;
		else
			q_new_digit4 = 0;

		if (q_new_digit4 == q_near_digit4){//if the input value is the same as before
			int cc = q_near->getCounter();
			q_new->setCounter(cc + 1);
		}
		else{
			q_new->setCounter(0);
		}
		*/
		Transition transition = tboot;
		nodes.push_back(q_new);
		eye->push(q_new, transition);
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
		if (forceIteration){
			if (i >= k){
				simulationFinished = true;
			}
		}
		else{
			if (i >= k || timeEnvlope >= sim_time){
				simulationFinished = true;
			}
		}
	}
}


void TimedRRT::worstCaseEyeDiagram(){
	generateMonteCarloInputSequence();
	config->getParameter("edu.uiuc.csl.core.sampling.iteration", &k);
	int i = nodes.size();
	while (i < k){
		GammaSimMode = false;
		node* q_near = eye->getNode(0);
		deltaSimulation(q_near);
		i++;
	}

	//worstCaseJitter();
}



void TimedRRT::deltaSimulation(node* q_near){
	GammaSimMode = 0;
	double* state_near = q_near->get();
	double* ic = new double[d];
	for (int j = 0; j<d; j++){
		ic[j] = state_near[j];
	}

	vector<string> settings;
	stringstream icInputFileName; icInputFileName << "ic_" << q_near->getIndex() << ".ic0";
	stringstream icOutputFileName; icOutputFileName << "ic_" << nodes.size() << ".ic";
	settings.push_back(icOutputFileName.str());
	settings.push_back(icInputFileName.str());
	settings.push_back("transient"); // or settings.push_back("dc");

	vector<double> param = generateSimulationParameters(q_near);
	double t_init = ic[d - 1];
	vector<double> result = system->simulate(ic, param, settings, 0, dt);

	result[0] += t_init;		//result[0] contains the time-stamp, in the simulation it is stamped as dt, however we have to add the time of the parrent node as well.

	node* q_new = new node(d);
	q_new->set(result);
	q_near->addChildren(q_new);		//add the new node to the tree
	q_new->setParent(q_near);		//We only make the parent-child releation ship during the tree build
	q_new->setIndex(nodes.size());
	q_new->setInputVector(param);

	Transition transition = tboot;
	nodes.push_back(q_new);
	eye->push(q_new, transition);
	delete ic;
}

void TimedRRT::worstCaseJitter(){
	for (int i = 0; i < 3;i++){
		cout << "jitter simulation #" << i << endl;
		node* q_near = eye->getNode(0);
		//generate a jitter data and transition data
		GammaSimMode = 1;
		//GammaJitter = generateUniformSample(0, 1e-11);
		//GammaTransition = generateUniformSample(0, 5e-11);

		double freq = 0; config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
		double period = 1 / freq;
		double jitterMax = 0; config->getParameter("edu.uiuc.csl.system.param.jitter.max", &jitterMax);
		jitterMax = jitterMax*1e-12;
		GammaJitter = generateUniformSample(0, 1e-11);
		GammaTransition = generateUniformSample(0, 2e-11);

		double t0 = q_near->getTime();
		int cycles = floorf(t0 / period);
		double t1 = t0 - period*cycles;

		if (t1 + GammaJitter >= jitterMax){
			GammaJitter = jitterMax - t1;
		}
		
		//determine how long the transition will take
		double eta = GammaJitter + GammaTransition;
		eta = 2e-11;
		int gammaIterations = 2*ceil(eta / dt);
		simulate(gammaIterations, q_near);
	}
}

void TimedRRT::simulate(double* initialState){
	root = new node(d);
	root->set(initialState);
	root->setRoot();
	nodes.push_back(root);
	root->setIndex(0);
	vector<double> p;
	for (int i = 0; i < var; i++)
		p.push_back(0);
	root->setInputVector(p);
	simulate(k, root);
}

//simulate the circuit for iter numbers from q_start
void TimedRRT::simulate(int iter, node* q_start){
	int i = 0;
	node* q_near = q_start;
	while (i < iter){
		double* state_near = q_near->get();
		double* state = new double[d];
		for (int j = 0; j<d; j++){
			state[j] = state_near[j];
		}
		vector<double> param = generateSimulationParameters(q_near);

		vector<string> settings;
		stringstream icInputFileName; icInputFileName << "ic_" << q_near->getIndex() << ".ic0";
		cout << endl <<  icInputFileName.str() << endl;
		stringstream icOutputFileName; icOutputFileName << "ic_" << nodes.size() << ".ic";
		settings.push_back(icOutputFileName.str());
		settings.push_back(icInputFileName.str());
		settings.push_back("transient");

		double t_init = state[d - 1];
		vector<double> result = system->simulate(state, param, settings, 0, dt);
		result[0] += t_init;		//result[0] contains the time-stamp, in the simulation it is stamped as dt, however we have to add the time of the parrent node as well.

		node* q_new = new node(d);

		q_new->set(result);
		q_near->addChildren(q_new);		//add the new node to the tree
		q_new->setParent(q_near);		//We only make the parent-child releation ship during the tree build
		q_new->setIndex(nodes.size());	//used to be setindex(i), which was buggy when loading rrt, took 2 hours to fix!
		q_new->setInputVector(param);

		Transition transition = tboot;
		nodes.push_back(q_new);

		if (GammaSimMode == 1){


			cout << "Sim#" << i << " " << q_new->getTime() << " vout=" << q_new->get(2) << endl;
			double freq = 0; config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
			double period = 1 / freq;
			double t = q_near->getTime();
			int cycles = floorf(t / period);
			double t1 = t - period*cycles;

			if (t1 < GammaJitter){
				//this node is still a candidate to start jittering
				q_new->setJitter(1);
			}
			else{
				q_new->setJitter(0);
			}
		}
		else{
			q_new->setJitter(0);
		}
		eye->push(q_new, transition);

		q_near = q_new;
		i++;
	}
}

void TimedRRT::addMonitor(Monitor* m){
	monitors.push_back(m);
}

double TimedRRT::getSimTime(){
	return sim_time;
}