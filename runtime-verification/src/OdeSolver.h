#pragma once
#include "System.h"
#include <iostream>
using namespace std;

class OdeSolver : public System {
public:
    OdeSolver();
    ~OdeSolver();
    void simulate();
};
