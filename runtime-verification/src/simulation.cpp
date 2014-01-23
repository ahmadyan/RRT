#include "simulation.h"
#include "pll.h"
#define	unluckyThirteen	13

bool Simulation::is_only_ascii_whitespace(const std::string& str){
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

double Simulation::unit(char u){
	switch(u){
	case 'f':
	case 'F':
		return 1e-15 ;
		break;
	case 'p':
	case 'P':
		return 1e-12 ;
		break;
	case 'n':
	case 'N':
		return 1e-9 ;
		break;
	case 'u':
	case 'U':
		return 1e-6 ;
		break;
	case 'm':
	case 'M':
		return 1e-3 ;
		break;
	case 'k':
	case 'K':
		return 1e3 ;
		break;

	default:
		return 1 ;
	}
}

//The input is something like this:
//5.0000000000n   226.6820055575m   71.3468319456m  -71.3468319456m
// The first number is transient simulation time, the rest are variables that need to be converted to double
vector<double> Simulation::parse(string s){
	vector<double> result ;
	string word ;
	string str = s ;
	if(s.c_str()[s.length()-1] == unluckyThirteen){	//guess what this does!
		str = s.substr(0, s.length()-1);
	}
	stringstream stream(str);
	while( getline(stream, word, ' ') ){
		if( !is_only_ascii_whitespace(word) ){
			double d ;
			stringstream ss(word);
			ss >> d ;
			char c = word.c_str()[word.length()-1] ;
			d = d*unit(c);
			result.push_back(d);
		}
	}
	return result ;
}
