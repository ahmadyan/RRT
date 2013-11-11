#include "node.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
using namespace std;

node::node(int _n){
	id = generateID();
	root=false; 
	n=_n;
	data = new double[n];
	for(int i=0;i<n;i++) data[i]=0;
}

//copy constructor
node::node(const node& Node){
	id = generateID() ;
	root = Node.root ;
	n = Node.n ;
	data = new double[n];
	for(int i=0;i<n;i++){
		data[i] = Node.data[i];
	}
}

//This piece of code causes unknown problems with visual studio 2012.
node::~node(){
	//for(int i=0;i<children.size();i++){
	// delete children[i];
	//}
	//	delete data ;
}


int node::objectCount = 0 ;

int node::generateID(){
	cout << "Generating a new Object :" << objectCount+1 << endl;
	return objectCount++;
}

int node::getID(){
	return id;
}

double node::unifRand(){
	return rand() / double(RAND_MAX);
}
double node::unifRand(double a, double b){
	return (b-a)*unifRand() + a;
}

void node::set(double* _x){
	//todo: check for size x
	data=_x;
}

double* node::get(){
	return data ;
}
double  node::get(int i){
	return data[i];
}

int node::getDimension(){
	return n;
}

double node::getTime(){
	return data[n-1];
}

//assigns a uniform random value to the new node,
//this methods is usually called for generatign a new sample in the RRT, so
//the RRT can find the nearest node in the tree to this newly generated node.
//random number is between (min[i], max[i]) range.
void node::randomize(double* min, double* max){
	for(int i=0;i<n;i++){
		data[i] = unifRand(min[i], max[i]) ;// (max[i]-min[i])*( rand()/double(RAND_MAX) ) + min[i] ;
	}
}

//this is a bit offsetted to push the rrt through time
// i = current sample #
// k = maximum sample #
// min, max = minimum and maximum value for each dim
void node::timed_randomize(int i, int k, double* min, double* max){
	for(int c=0;c<n-1;c++){
		data[c] = unifRand(min[c], max[c]); //(max[c]-min[c])*( rand()/double(RAND_MAX) ) + min[c] ;
	}

	cout << max[n-1] << endl;
	if(i>2000){//1000-3000 uniform
		data[n-1]=unifRand(0,max[n-1]);
	}else if(i>1000){	//1000-2000  //standard
		double offset = ( (double)i/(double)k ) *unifRand(min[n-1], max[n-1]) ;
		data[n-1] = unifRand(  offset , max[n-1]);
	}else if (i>100){	//0-1000 //rapid-growth
		data[n-1] = (double)i/(double)k * max[n-1];
	}else{
		data[n-1] = max[n-1];
	}
	//cout << "trand = " << i << " " << offset << " "  << data[n-1] << endl ;

}


//compute the norm between current node and given node a
double node::distance(node* a, double* max, double* min){
	double d = 0;
	for(int i=0;i<n;i++){
		d +=  ( ( data[i]-a->data[i] ) * ( data[i] - a->data[i] ) / (max[i]-min[i]));
	}
	d = sqrt(d);
	return d;
}

//compute the norm between current node and given node a
//d = || x-x' ||_p + a.dt
double node::timed_distance(node* a, double* max, double* min){
	//if( data[n-1] < a.data[n-1] ){
	//    return 99999999;  //unfourtunately we cannot go back in time, yet.
	//}
	double alpha = 1 ;
	double beta  = 0.01 ;
	double d = 0;
	for(int i=0;i<n-1;i++){
		d +=  ( ( data[i]-a->data[i] ) * ( data[i] - a->data[i] ) / ((max[i]-min[i])*(max[i]-min[i])));
	}
	d = beta * sqrt(d);
	d += alpha * fabs((data[n-1] - a->data[n-1])/(max[n-1]-min[n-1]));
	return d;
}

//TODO: issue a warning when the vector size does not match in any of overloaded operators
node& node::operator=(const node &rhs){
	if(this==&rhs) return *this; 
	delete data ;   //deallocate previous copy of data, (we can's reuse it, since the n may be different between two node)
	n = rhs.n;
	data = new double[rhs.n];
	for(int i=0;i<n;i++){
		data[i] = rhs.data[i];
	}
	return *this;
}

node& node::operator+=(const node &rhs){
	for(int i=0;i<n;i++){
		data[i] += rhs.data[i] ;
	}
	return *this;
}

node& node::operator-=(const node &rhs){
	for(int i=0;i<n;i++){
		data[i] -= rhs.data[i] ;
	}
	return *this;
}

const node node::operator+(const node &rhs) const{
	return node(*this) += rhs ;
}

const node node::operator+(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator+(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator-(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator-(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator*(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] *= rhs ;
	}
	return result;
}

const node node::operator*(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] *= rhs ;
	}
	return result;
}

//result of multiplying two vector is a number.
const double node::operator*(const node &rhs) const{
	double result = 0 ;
	for(int i=0;i<n;i++){
		result += (this->data[i] * rhs.data[i]) ;
	}
	return result;
}

const node node::operator-(const node &rhs) const{
	return node(*this) -= rhs ;
}

string node::toString(){
	stringstream nstr ;
	for(int i=0;i<n;i++){
		nstr << data[i] << " " ;
	}
	return nstr.str();
}
void node::dump(){
	for(int i=0;i<n;i++){
		cout << data[i] << " " ;
	}
	cout << endl ;
}

void node::addChildren(node* k){
	children.push_back(k);
}

vector<node*> node::getChildren(){
	return children;
}

node* node::getChild(int i){
	return children[i];
}

int node::getSize(){
	return children.size();
}

void node::setParent(node* p){
	parent = p ;
}

node* node::getParent(){
	return parent;
}

bool node::isRoot(){
	return root ;
}

void node::setRoot(){
	root=true;
}

//recursively returns the nearestNode between this node and all of it's children w.r.t. q_sample
pair<node*, double> node::getNearestNode(node* q_sample, double* min, double* max, bool timed){
	node* nearestNode = this ;
	double distance = ( timed ? q_sample->timed_distance(this, min, max) : q_sample->distance(this, min, max) );
	for(int i=0;i<children.size(); i++){
		pair<node*, double> p = q_sample->getNearestNode(children[i], min, max, timed); 
		if(p.second < distance ){
			distance = p.second ;
			nearestNode = p.first ;
		}
	}
	return make_pair(nearestNode, distance);
}

void node::save(ofstream& of){
	cout << "Printing node " << getID() << endl ;
	if(root){
		of << getID() << " -1 " <<  toString() << endl ;
	}else{
		of << getID() << " " << getParent()->getID() << " " <<  toString() << endl ;
	}
	
	for(int i=0;i<children.size();i++){
		cout << "Child " << i << endl ;
		getChild(i)->save(of);
	}
}