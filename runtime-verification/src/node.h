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
	vector<node*> cast;			//the casts are at the same state space, but in different time
	int n ; // n=node-dimension
	double* data;  // this array holds the state-time space, n-1 double for holding x and last double for holding time
	bool root ;
	int id;
    static int objectCount ;
	int index; //this node index in the nodes array
public:
	node(int n);
	node(const node&);
	node::node(int _n, int _id, double* _data);
	~node();
	void set(double* _data);
	void set(int, double);
	double* get();
	double  get(int i);
	double  getTime();
	int getDimension();
	void randomize(double* min, double* max);
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
	void addCast(node*);
	void setParent(node*); 
	vector<node*> getChildren();
	vector<node*> getCast();
	node* getChild(int i);
	node* getCast(int);
	int getSize();
	node* getParent();
	bool isRoot();
	void setRoot();
	pair<node*, double> getNearestNode(node*, double*, double*, bool);

	static int generateID();
	int getID();
	int getIndex();
	void setIndex(int);
	void save(ofstream&);
};


