#pragma once

#include <vector>
#include <string>
#include "point.h"
#include "node.h"
#include "config.h"
#include "circuit.h"
#include "kdtree.h"

using namespace std;



class Forest{
	int d;	//number of dimensions
	int k;	//total iterations
	int forward_rrt_count;	
	int backward_rrt_count;
	int connector_rrt_count;
    vector<Node*> nodes ;
    struct kdtree *kd;
    struct kdres  *set;

    Configuration* config;
    Circuit* circuit;
    void collision(Node* q_new, Node* q_near);
public:
	Forest(Configuration* config, Circuit* circuit);
	~Forest();
    
    Node* find_nearest_point(geometry::Point* ) ;
    pair<geometry::Point*, Type> find_optimum_trajectory(geometry::Point* , Node* );
    Node* addNode(geometry::Point* , Node* parent, Type t);
    Node* addRoot(geometry::Point*, Type);
    Node* addRoot(Node* new_node, Type type);
    int getDimension();
    
    void load(string fileName);
    void save(string fileName);
    string toString();
    string draw();
    string draw(int d0, int d1);
    string draw(int index);
    int getSize();
    Node* getNode(int i);
};