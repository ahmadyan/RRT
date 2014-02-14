#include <iostream>

using namespace std;

int main(){
    int k=7; //number of inverter in the ring
    string prefix="edu.uiuc.csl.";
    //generating random vdd
    double vdd=0.9;
    double dvdd=0.1;
    for(int i=0;i<k;i++){
        cout << "#auto-generated parameter: vdd#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=" << vdd-dvdd << endl ;
        cout << prefix << "system.param.max["<< i<<"]=" << vdd+dvdd << endl ;
    }
    cout << endl ;

    //generating random vgnd
    double gnd=0;
    double dgnd=0.1;
    
    for(int i=1*k; i<2*k;i++){
        cout << "#auto-generated parameter: vgnd#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=" << gnd-dgnd << endl ;
        cout << prefix << "system.param.max["<< i<<"]=" << gnd+dgnd << endl ;
    }
    cout << endl ;

    //generating random vnoise
    double vnoise=0;
    double dvnoise=0.05;
    for(int i=2*k; i<3*k;i++){
        cout << "#auto-generated parameter: vnoise#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=" << vnoise-dvnoise<< endl ;
        cout << prefix << "system.param.max["<< i<<"]=" << vnoise+dvnoise<< endl ;
    }

    cout << endl ;

    //generating random vsubp
    double vsubp=0.9;
    double dvsubp=0.1;
    for(int i=3*k; i<4*k;i++){
        cout << "#auto-generated parameter: vsubp#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]="<<vsubp-dvsubp << endl ;
        cout << prefix << "system.param.max["<< i<<"]="<< vsubp+dvsubp << endl ;
    }
    cout << endl ;

    //generating random vsubn
    double vsubn=0;
    double dvsubn=0.1;
    for(int i=4*k; i<5*k;i++){
        cout << "#auto-generated parameter: vsubn#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=" << vsubn-dvsubn<< endl ;
        cout << prefix << "system.param.max["<< i<<"]="<< vsubn+dvsubn<< endl ;
    }

    cout << endl ;
    int var=7;
    for(int i=1;i<=var;i++){
        cout << "#auto-generated variable: var" << i << endl;
        cout << prefix << "system.var.ic["<< i<<"]=0"<< endl;
        cout << prefix << "system.var.min["<< i<<"]=-0.2"<< endl;
        cout << prefix << "system.var.max["<< i<<"]=1.1"<< endl;
        cout << prefix << "system.var.name["<< i<<"]=v"<< i << endl;
    }
}