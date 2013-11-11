#pragma once
#include <string>
#include <sstream>
#include <vector>

using namespace std;
// node class contains the tupe of points in n-1 dimension plus time, 
// each tuple is defined as <x_1, x_2, ..., x_{n-1}, t>
// the distance between two node is defined as ||x-x'||_p + \alpha \delta t


class node{
	node* parent ;
	vector<node*> children ;
	int n ; // n=node-dimension
	double* data;  // this array holds the state-time space, n-1 double for holding x and last double for holding time
	bool root ;
	int id;
    static int objectCount ;
public:
	node(int n);
	node(const node&);
	~node();
	void set(double* _data);
	double* get();
	double  get(int i);
	double  getTime();
	int getDimension();
	void randomize(double* min, double* max);
	void timed_randomize(int i, int k, double* min, double* max);
	double distance(node* a, double* max, double* min);
	double timed_distance(node* a, double* max, double* min);

	double unifRand();
	double unifRand(double a, double b);
	void dump();
	string toString();
	node& operator=(const node &rhs);
	node& operator+=(const node &rhs);
	node& operator-=(const node &rhs);
	const node operator+(const node &rhs) const;
	const node operator+(const int &rhs) const;    
	const node operator+(const double &rhs) const;
	const node operator-(const node &rhs) const; 
	const node operator-(const int &rhs) const; 
	const node operator-(const double &rhs) const; 
	const node operator*(const int &rhs) const;
	const node operator*(const double &rhs) const;
	const double operator*(const node &rhs) const;

	void addChildren(node*);
	void setParent(node*); 
	vector<node*> getChildren();
	node* getChild(int i);
	int getSize();
	node* getParent();
	bool isRoot();
	void setRoot();
	pair<node*, double> getNearestNode(node*, double*, double*, bool);

	static int generateID();
	int getID();
	void save(ofstream&);
};


