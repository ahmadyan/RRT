#include "RRT.h"
#include <stdio.h>
#include <ctime>
#include <random>
#include <iostream>
#include <limits>
#include <stack>
#include <chrono>
using namespace std;

RRT::RRT(Configuration* c, int _d, int _k, int _var, string nam){
	d = _d;
	k = _k;
	var = _var;
	name = nam;
	min = new double[d];
	max = new double[d];
	variationMin = new double[var];
	variationMax = new double[var];
	//default values for min-max
	for (int i = 0; i < d; i++){
		min[i] = 0;
		max[i] = 1;
	}
	for (int i = 0; i < var; i++){
		variationMin[i] = 0;
		variationMax[i] = 0;
	}

	eye = new EyeDiagram(c);
}

RRT::RRT(Configuration* c, string fileName){
	eye = new EyeDiagram(c);
	load(fileName);
}

RRT::~RRT(){
	//if(min!=NULL) delete min ;
	//if(max!=NULL) delete max ;
}


double RRT::generateNormalSample(double mean, double std){
	//A Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed1);
	std::normal_distribution<double> normal(mean, std);
	return normal(generator);

}

double RRT::generateTruncatedNormalSample(double mean, double std, double min, double max){
	double s = generateNormalSample(mean, std);
	if (s < min) s = min;
	if (s > max) s = max;
	return s;
}

double RRT::generateUniformSample(){
	return rand() / double(RAND_MAX);
}

double RRT::generateUniformSample(double a, double b){
	return (b - a)*generateUniformSample() + a;
}


//There is no actual simulation in this code,
//except a pretty cool operator overloading!
void RRT::buildUniform(double* initialState){
	root = new node(d);
	root->set(initialState);
	root->setRoot();
	for (int i = 0; i < k; i++){
		//create a new sample
		node* q_sample = new node(d);
		q_sample->randomize(min, max);
		node* q_near = (getNearestNode(q_sample))[0];
		//find the best trajectory using integration or simulation, constructing the new node q_new
		double delta = 0.5;
		node q_new = *q_near + ((*q_sample - *q_near) * delta);
		q_new.dump();
		q_near->addChildren(&q_new);
		q_new.setParent(q_near);
	}
}

void RRT::build(double* initialState){
	root = new node(d);
	root->set(initialState);
	root->setRoot();

	for (int i = 0; i < k; i++){
		node* q_sample = new node(d);
		q_sample->randomize(min, max);
		node* q_near = (getNearestNode(q_sample))[0];

		double* state_near = q_near->get();
		double* ic = new double[d];
		for (int j = 0; j < d; j++) ic[j] = state_near[j];
		double delta = 0.5;
		//variation or input to the system
		//todo: should be defined in the main or system, not here
		vector<double> param;
		for (int i = 0; i < var; i++){
			param.push_back(generateUniformSample(0.29, 0.31));
		}
		vector<string> settings;
		vector<double> result = system->simulate(ic, param, settings, 0, delta);

		double* state = new double[d];
		for (int i = 0; i < d; i++){
			state[i] = result[i + 1];
		}
		delete ic;
		cout << "state*=" << state[0] << " " << state[1] << endl;
		node* q_new = new node(d);
		q_new->set(state);

		q_near->addChildren(q_new);
		q_new->setParent(q_near);
	}
}


//Define the minimum and maximum value for each dimension
void RRT::setBound(int i, double _min, double _max){
	min[i] = _min;
	max[i] = _max;
}

void RRT::setVariationBound(int i, double _min, double _max){
	variationMin[i] = _min;
	variationMax[i] = _max;
}

/*
node RRT::getNode(int i){
return nodes[i].second;
}

int RRT::getNodeParentId(int i){
return nodes[ nodes[i].first ].first;
}

//returns the parent of node i in the tree, we hold the parent# in nodes[i].first
node RRT::getParentNode(int i){
cout << "get parent for node " << i << " = " << nodes[i].first << endl ;
return nodes[ nodes[i].first ].second;
}

int RRT::sampleNumber(){
return nodes.size();
}
*/
int RRT::getDimension(){
	return d;
}

void RRT::setSystem(System* s){
	system = s;
}

double RRT::getMin(int i){
	return min[i];
}

double RRT::getMax(int i){
	return max[i];
}

//Saving the RRT into a file
//Stores each-node-id, parent-node-id, node-data, input
void RRT::save(string fileName){
	//cout << "Saving the RRT" << endl;
	ofstream file;
	file.open(fileName.c_str());
	file << "rrt" << endl;
	file << d << endl;		//number of dimensions
	file << var << endl;
	file << k << endl;		//number of nodes

	//saving the bounds on each dimensions
	for (int i = 0; i < d; i++){
		file << min[i] << endl;
		file << max[i] << endl;
	}

	//start from the root and recursively print every node
	//root->save(file);

	//linear printing of every nodes
	//nodes are being printed in this format: id parentid data timestamp

	//cout << "Nodes.size=" << nodes.size() << endl;
	for (int i = 0; i < nodes.size(); i++){
		//cout << "Saving node " << i << endl;
		file << i << " ";

		if (nodes[i]->isRoot()){
			file << "-1" << " ";
		}
		else{
			for (int j = 0; j < nodes.size(); j++){
				if (nodes[i]->getParent()->getID() == nodes[j]->getID()){
					file << j << " ";
					break;
				}
			}
		}

		for (int j = 0; j < nodes[i]->getDimension(); j++){
			file << nodes[i]->get(j) << " ";
		}

		for (int j = 0; j < var; j++){
			file << nodes[i]->getInput(j) << " ";
		}

		file << endl;
	}

	file.close();
}

void RRT::loadInput(string fileName){
	string line;
	ifstream file(fileName.c_str());
	int id;
	int bit;
	double jitter;
	double transition;
	if (file.is_open()){
		for (int i = 0; i < k; i++){
		file >> id;
		file >> bit;
		file >> jitter;
		file >> transition;
		bits.push_back(bit);
		}
	}
}

//Loads the RRT from the text file (default extension is the *.rrt)
void RRT::load(string fileName){
	
	if (config->checkParameter("edu.uiuc.crhc.core.options.eyediagram.bitstream", "1")){
		loadInput(config->get("edu.uiuc.crhc.core.options.eyediagram.inputfile"));
		eye->setBits(bits);
	}

	string line;
	ifstream file(fileName.c_str());

	if (file.is_open()){
		file >> name;	//The first line should read rrt
		cout << "Loading RRT " << name << endl;
		file >> d;
		file >> var;
		file >> k;
		min = new double[d];
		max = new double[d];
		cout << "d,k=" << d << " " << k << endl;
		for (int i = 0; i < d; i++){
			file >> min[i];
			file >> max[i];
			cout << min[i] << " " << max[i] << endl;
		}

		for (int i = 0; i < k; i++){
			cout << i << endl;
			double* data = new double[d];
			int id = -2;
			int parent_id = -2;
			file >> id;
			file >> parent_id;
			cout << "id :" << id << " , pid=" << parent_id << " ";
			cout << "data: "; 
			for (int j = 0; j < d; j++){
				file >> data[j];
				cout << data[j] << ", ";
			}
			cout << endl;
			vector<double> param;
			if (config->checkParameter("edu.uiuc.crhc.core.options.input.format", "2")){
				for (int j = 0; j < var; j++){
					double x;
					file >> x;
					param.push_back(x);
				}
			}
			
			node* newNode = new node(d, id, data);
			cout << newNode->toString() << endl;
			newNode->setInputVector(param);
			//If this node is a root node (i.e. the parent_id is -1), sets this as root, otherwise
			//this node has a parent. Find the parent and add this as children. 
			//I'm reading/writing nodes in a monotonic manner, so I can directly access nodes via nodes[id]
			if (parent_id == -1){
				newNode->setRoot();
				root = newNode;
			}else{
				newNode->setParent(nodes[parent_id]);
				nodes[parent_id]->addChildren(newNode);
			}
			if (i==201 || i==401){
				newNode->setJitter(1);
			}
			else{
				newNode->setJitter(0);
			}
			newNode->setIndex(id);	//This is unnecessary, because currently nodes are constructed with their id
			nodes.push_back(newNode);
			if (config->checkParameter("edu.uiuc.crhc.core.options.eyediagram", "1"))
				eye->push(newNode);
		}
		file.close();
	}
}


string RRT::getName(){
	return name;
}

void RRT::setdt(double d){
	dt = d;
}

double RRT::getdt(){
	return dt;
}

node* RRT::getRoot(){
	return root;
}

vector<node*> RRT::recursiveNearestNodeSearch(node* q_sample){
	vector<node*> results;
	results.push_back(root->getNearestNode(q_sample, min, max, false).first);
	return results;
}

vector<node*> RRT::getLatestNode(){
	vector<node*> results;
	//waste this sample on forward progress in time. 
	//we choose the sample with highest time stamp, regardless of how close that sample is to the q_sample
	//Can be efficiently implemented using priority queue
	double lastTime = -1;
	node* lastNode;
	int last;
	for (int i = 0; i < nodes.size(); i++){
		if (nodes[i]->getTime() > lastTime){
			last = i;
			lastTime = nodes[i]->getTime();
			lastNode = nodes[i];
		}
	}
	results.push_back(lastNode);
	return results;
}

vector<node*> RRT::NearestNodeInProjectiveSpace(node* q_sample){
	vector<node*> results;
	double closestDistance = 99999;
	int closestNodeIndex = -1;
	for (int i = 0; i < nodes.size(); i++){
		//Compute distance from an abstraction, not an entire model, 
		//I don't really care about every variable in the circuits, some are more important than the others (weighted distance model???).
		/*int projectiveDimension = 1; config->getParameter("edu.uiuc.csl.core.search.space.subset", &projectiveDimension);
		double distance = 0;
		for (int j = 0; j < projectiveDimension; j++){
			distance += (nodes[i]->get(j) - q_sample->get(j))*(nodes[i]->get(j) - q_sample->get(j));
		}
		distance = sqrt(distance);
		if (distance < closestDistance){
			closestDistance = distance;
			closestNodeIndex = i;
		}*/
		double iv = 2;
		double d = abs(nodes[i]->get(iv) - q_sample->get(iv));
		if (d<closestDistance){
			closestDistance = d;
			closestNodeIndex = i;
		}
	}
	results.push_back(nodes[closestNodeIndex]);
	return results;
}

vector<node*> RRT::randomizedNearestNodeSearchTimeDistance(node* q_sample){
	vector<node*> results;
	double minDistance = 99999;
	int nearest = -1;

	int searchSize = 2*log(nodes.size())+2;
	for (int i = 0; i < searchSize; i++){
		int s = generateUniformSample(0, nodes.size()-1);
		double d = q_sample->timed_distance(nodes[s], min, max);
		if (d < minDistance){
			minDistance = d;
			nearest = s;
		}
	}
	results.push_back(nodes[nearest]);
	return results;
}

vector<node*> RRT::NearestNodeUsingTimeDistance(node* q_sample){
	vector<node*> results;
	double minDistance = 99999;
	int nearest = -1;
	for (int i = 0; i < nodes.size(); i++){
		double d = q_sample->timed_distance(nodes[i], min, max);

		if (d < minDistance){
			minDistance = d;
			nearest = i;
		}
		//double eTol = 1e-6; config->getParameter("edu.uiuc.csl.core.search.tol", &eTol);
		//if (d <= eTol){
		//	if (nodes[i]->getID() != q_sample->getID()){
		//		results.push_back(nodes[i]);
		//	}
		//}
	}
	results.push_back(nodes[nearest]);
	return results;
}


vector<node*> RRT::randomizedNearestNodeSearchEuclideanDistance(node* q_sample){
	vector<node*> results;
	double minDistance = 99999;
	int nearest = -1;

	int searchSize = log(nodes.size())+2;
	for (int i = 0; i < searchSize; i++){
		int s = generateUniformSample(0, nodes.size());
		double d = q_sample->distance(nodes[s], min, max);
		if (d < minDistance){
			minDistance = d;
			nearest = i;
		}
	}
	results.push_back(nodes[nearest]);
	return results;
}

vector<node*> RRT::NearestNodeUsingEuclideanDistance(node* q_sample){
	vector<node*> results;
	double minDistance = 99999;
	int nearest = -1;
	for (int i = 0; i < nodes.size(); i++){
		double d = q_sample->distance(nodes[i], min, max);
		if (d < minDistance){
			minDistance = d;
			nearest = i;
		}
		//double eTol = 1e-6; config->getParameter("edu.uiuc.csl.core.search.tol", &eTol);
		//if (d <= eTol){
		//	if (nodes[i]->getID() != q_sample->getID()){
		//		results.push_back(nodes[i]);
		//	}
		//}
	}
	results.push_back(nodes[nearest]);
	return results;
}

vector<node*> RRT::getNearestNode(node* q_sample){
	vector<node*> results;

	if (config->checkParameter("edu.uiuc.csl.core.search.mode", "random")){
		int n = generateUniformSample(0, nodes.size());
		results.push_back(nodes[n]);
		return results;
	}

	if (config->checkParameter("edu.uiuc.csl.core.search.mode", "optimize")){
		int n = generateUniformSample(0, eye->getWindowSize());
		node* v = eye->getNode(n);
		results.push_back(v);
		return results;
	}

	if (config->checkParameter("edu.uiuc.csl.core.search.time-progress", "1")){
		if (q_sample->getTime() != -1){
			//time-progress is enabled
			double p = generateUniformSample(0, 1);
			double progressFactor; config->getParameter("edu.uiuc.csl.core.sampling.progressFactor", &progressFactor);
			if (p < progressFactor){
				return getLatestNode();
			}
		}
	}

	//follow the normal search
	//should we search the total state space or just a projection?
	if (config->checkParameter("edu.uiuc.csl.core.search.space", "total")){
		//which distance metric to be used for the search?
		if (config->checkParameter("edu.uiuc.csl.core.search.distance", "timed")){
			//timed-distance metric
			if (config->checkParameter("edu.uiuc.csl.core.search.randomized", "1")){
				return randomizedNearestNodeSearchTimeDistance(q_sample);
			}
			else{
				return NearestNodeUsingTimeDistance(q_sample);
			}
		}
		else{
			//eucledian metric
			if (config->checkParameter("edu.uiuc.csl.core.search.randomized", "1")){
				return randomizedNearestNodeSearchEuclideanDistance(q_sample);
			}
			else{
				return NearestNodeUsingEuclideanDistance(q_sample);
			}
		}
	}
	else if (config->checkParameter("edu.uiuc.csl.core.search.space", "subset")){
		return NearestNodeInProjectiveSpace(q_sample);
	}
	else{
		cout << "Uknown search algorithm selected" << endl;
	}
	return results;
}

void RRT::addMonitor(Monitor* m){
	cout << "Adding a new monitor to the RRT , total monitors=" << monitors.size() << endl;
	monitors.push_back(m);
}

void RRT::setConfig(Configuration* c){
	config = c;
}

EyeDiagram* RRT::getEyeDiagram(){
	return eye;
}

void RRT::setIterations(int _k){
	k = _k;
}

vector<node*> RRT::getTest(node* v){
	vector<node*> results;
	node* x = v;
	while (!x->isRoot()){
		results.push_back(x);
		x = x->getParent();
	}

	return results;
}

string RRT::drawTest(vector<node*> path, int color){
	if (path.size() == 0) return "";
	stringstream str;

	
	for (int i = 0; i < path.size()-1; i++){
		/*double t = path[i]->getTime();
		int cycles = ceilf(t / period) - 1;
		if (cycles == -1)cycles = 0;
		double t1 = t - period*cycles;

		t = path[i+1]->getTime();
		cycles = ceilf(t / period) - 1;
		if (cycles == -1)cycles = 0;
		double t2 = t - period*cycles;
		*/


		double iToX = path[i]->getTime();
		double iFromX = path[i + 1]->getTime();
		double iToY = path[i + 1]->getInput(0); //path[i+1]->get(voltage);
		double iFromY = path[i]->getInput(0); //path[i]->get(voltage);
		cout << "[" << iFromX << "," << iFromY << "] --> [" << iToX<< "," << iToY << "]"<<endl;
		if (path[i]->getInput(0) == 0.9){
			str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"grey\" lw 1 \n";
		}
		else{
			if (color == 1){
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"red\" lw 1 \n";
			}
			else{
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"blue\" lw 1 \n";
			}
			
		}
		
	}
	return str.str();
}