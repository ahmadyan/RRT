#include "System.h"
#define USE_HSPICE	true

System::System(){}
System::~System(){}

void System::setSystem(int _d,int (*f)(double t, const double y[], double f[], void *params), 
                    int (*j) (double t, const double y[], double *dfdy, double dfdt[], void *params) ){
    d=_d;
    function = f ;
    jacobian = j ;
}

//void System::dumpSystemInfo(){
//}

//this is a simple integrator that would do a step integration
//double System::simulate(double* initialState, double* param, double dt){
//	if(USE_HSPICE)
//		return runSPICE(initialState, param, dt);
//	else
//		return integareODE(initialState, param, dt);
//}

double System::integareODE(double* initialState, double* param, double dt){
/*     double mu = param;
    gsl_odeiv2_system dynamic = {function, jacobian, d , &mu};  // {vector-element, vector-of-derivative, system-dimension, arbitrary-parameters}
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&dynamic, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t=0;
    int status = gsl_odeiv2_driver_apply (d, &t, dt, initialState);
    if (status != GSL_SUCCESS){
        cout << "error, return value=" <<  status << endl ;
        return 0;
    }
    cout <<  t  << " " << initialState[0] << " " <<  initialState[1] << endl ;
    gsl_odeiv2_driver_free (d);
    return t;*/
	return 0 ;
}

//double System::runSPICE(double* initialState, double* param, double dt){
	/*spice hspice;
	vector<double> result = hspice.simulate(type, initialState, param, dt);
	cout << "result.size()=" << result.size() << endl;
	for (int i = 0; i < result.size(); i++){
		cout << i << " " << result[i] << endl;
	}
	for(int i=0;i<result.size()-1;i++){
		initialState[i]=result[i+1];
		cout << i << " " << result[i+1];
	}
	return result[0];*/
//}


//this is a compelete integrator, does not have any particular use, except for debugging the previous function
//double System::simulate(){
/*    double mu = 10;
    gsl_odeiv2_system dynamic = {function, jacobian, 2, &mu};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&dynamic, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 100.0;
    double y[2] = { 1.0, 0.0 };
    for (int i = 1; i <= 100; i++){
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }
    gsl_odeiv2_driver_free (d);
    return t1 ;*/
//	return 0 ;
//}

void System::setDimension(int _d){
    d=_d;
}

int System::getDimension(){
    return d;
}
#define	unluckyThirteen	13

bool System::is_only_ascii_whitespace(const std::string& str){
	//bool isOnlyWhiteSpace = false;

	//return str.find_first_not_of (' ') == str.npos)
	std::string::const_iterator  it = str.begin();
	do {
		if (it == str.end()) return true;
	} while (*it >= 0 && *it <= 0x7f && std::isspace(*(it++)));
	// one of these conditions will be optimized away by the compiler,
	// which one depends on whether char is signed or not
	return false;
}

double System::unit(char u){
	switch (u){
	case 'f':
	case 'F':
		return 1e-15;
		break;
	case 'p':
	case 'P':
		return 1e-12;
		break;
	case 'n':
	case 'N':
		return 1e-9;
		break;
	case 'u':
	case 'U':
		return 1e-6;
		break;
	case 'm':
	case 'M':
		return 1e-3;
		break;
	case 'k':
	case 'K':
		return 1e3;
		break;

	default:
		return 1;
	}
}

//The input is something like this:
//5.0000000000n   226.6820055575m   71.3468319456m  -71.3468319456m
// The first number is transient simulation time, the rest are variables that need to be converted to double
vector<double> System::parse(string s){
	vector<double> result;
	string word;
	string str = s;
	if (s.c_str()[s.length() - 1] == unluckyThirteen){	//guess what this does!
		str = s.substr(0, s.length() - 1);
	}
	stringstream stream(str);
	while (getline(stream, word, ' ')){
		if (!is_only_ascii_whitespace(word)){
			double d;
			stringstream ss(word);
			ss >> d;
			char c = word.c_str()[word.length() - 1];
			d = d*unit(c);
			result.push_back(d);
		}
	}
	return result;
}
