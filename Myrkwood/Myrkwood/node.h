#pragma once

#include "point.h"
#include "circuit.h"
#include <set>
enum Type {forward, backward, coupling, invalid};


class Node{
	geometry::Point* point;
    Node* parent;
    bool root;
    Type nodeType ;
    Type treeType ;
    double t;   //time-augmentation for time-variant systems
public:
	Node(geometry::Point*, bool);
	~Node();
    bool isRoot() const;
    double distance(Node* a, Circuit* circuit);
    double distance(geometry::Point* p, Circuit* circuit);
    double* getState();
    Node* getParent();
    void setParent(Node* p);
    Type getNodeType();
    Type getTreeType();
    void setNodeType(Type);
    void setTreeType(Type);
    string toString();
    double getData(int);
    geometry::Point* getPoint();
    void setTime(double);
    double getTime();
};