//
//  AbstractRRT.h
//  core
//
//  Created by Adel Ahmadyan on 2/25/13.
//  Copyright (c) 2013 University of Illinois at Urbana-Champaign. All rights reserved.
//
#pragma once
#include <iostream>
#include <string.h>

using namespace std;

class AbstractRRT{
public:
    virtual void setIntegrator()=0;
    virtual void addRoot()=0;
    virtual void addNode(Node*)=0;
    virtual Node* getNearestNode(Node*)=0;
    virtual string draw()=0;
    virtual int getSize()=0;
    virtual int getDimension()=0;
};