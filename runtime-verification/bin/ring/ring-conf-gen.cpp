#include <iostream>

using namespace std;

int main(){
    int k=7; //number of inverter in the ring
    string prefix="edu.uiuc.csl.";
    //generating random vdd
    for(int i=0;i<k;i++){
        cout << "#auto-generated parameter: vdd#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=0.7" << endl ;
        cout << prefix << "system.param.max["<< i<<"]=1.1" << endl ;
    }
    cout << endl ;
    //generating random vgnd
    for(int i=1*k; i<2*k;i++){
        cout << "#auto-generated parameter: vgnd#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=-0.2" << endl ;
        cout << prefix << "system.param.max["<< i<<"]=0.2" << endl ;
    }
    cout << endl ;
    //generating random vnoise
    for(int i=2*k; i<3*k;i++){
        cout << "#auto-generated parameter: vnoise#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=-0.01" << endl ;
        cout << prefix << "system.param.max["<< i<<"]=0.01" << endl ;
    }

    cout << endl ;
    //generating random vsubp
    for(int i=3*k; i<4*k;i++){
        cout << "#auto-generated parameter: vsubp#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=0.7" << endl ;
        cout << prefix << "system.param.max["<< i<<"]=1.1" << endl ;
    }
    cout << endl ;

    //generating random vsubn
    for(int i=4*k; i<5*k;i++){
        cout << "#auto-generated parameter: vsubn#" << i << endl;
        cout << prefix << "system.param.min["<< i<<"]=-0.2" << endl ;
        cout << prefix << "system.param.max["<< i<<"]=0.2" << endl ;
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