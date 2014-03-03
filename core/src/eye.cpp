#include "eye.h"

//The constructor for the EyeDiagram classes, 
//this constructor will get the initial parameter from the config file, then allocate the eyelids and sets the initial values for max/min
EyeDiagram::EyeDiagram(Configuration* c){
	config = c;

	config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
	config->getParameter("edu.uiuc.csl.core.simulation.dt", &sampleRate);
	config->getParameter("edu.uiuc.csl.core.simulation.window", &window);
	config->getParameter("edu.uiuc.crhc.core.options.eyediagram.var", &voltage);
	period = 1 / freq;

	size = window / sampleRate;

	palpebraSuperior = vector< vector<node*> >(size);

	maxSuperior = new double[size];
	minSuperior = new double[size];
	maxSuperiorIndex = new int[size];
	minSuperiorIndex = new int[size];

	palpebraInferior = vector< vector<node*> >(size);
	maxInferior = new double[size];
	minInferior = new double[size];
	maxInferiorIndex = new int[size];
	minInferiorIndex = new int[size];

	for (int i = 0; i < size; i++){
		maxSuperior[i] = -9909;
		minSuperior[i] = 9909;
		maxSuperiorIndex[i] = -1;
		minSuperiorIndex[i] = -1;
		maxInferior[i] = -9909;
		minInferior[i] = 9909;
		maxInferiorIndex[i] = -1;
		minInferiorIndex[i] = -1;
	}
}

EyeDiagram::~EyeDiagram(){
	delete maxSuperior;
	delete maxInferior;
	delete minSuperior;
	delete minInferior;
	delete maxSuperiorIndex;
	delete minSuperiorIndex;
	delete maxInferiorIndex;
	delete minInferiorIndex;
}

//Perform the initial integration of the eyediagram after applying the 00110 pattern
void EyeDiagram::sum(){
	double palpebraSuperiorInsideIntegral = 0;
	double palpebraSuperiorOutsideIntegral = 0;
	double palpebraInferiorInsideIntegral = 0;
	double palpebraInferiorOutsideIntegral = 0;

	//todo: these should be a vector, not a fixed value
	double NadirSuperiorHigh = 1.1;
	double NadirSuperiorLow = -0.2;
	double NadirInferiorHigh = 1.1;
	double NadirInferiorLow = -0.2;

	for (int i = 0; i < size; i++){
		palpebraSuperiorInsideIntegral += (minSuperior[i] - NadirSuperiorLow )* sampleRate;
		palpebraSuperiorOutsideIntegral += (NadirSuperiorHigh - maxSuperior[i] ) * sampleRate;
		palpebraInferiorInsideIntegral += (NadirInferiorHigh - maxInferior[i])* sampleRate;
		palpebraInferiorOutsideIntegral += (minInferior[i] - NadirInferiorLow) * sampleRate;
	}

	cout << "Result" << endl << endl;
	cout << "palpebraSuperiorInsideIntegral=" << palpebraSuperiorInsideIntegral << endl;
	cout << "palpebraSuperiorOutsideIntegral=" << palpebraSuperiorOutsideIntegral << endl;
	cout << "palpebraInferiorInsideIntegral=" << palpebraInferiorInsideIntegral << endl;
	cout << "palpebraInferiorOutsideIntegral=" << palpebraInferiorOutsideIntegral << endl;
	cout << "========" << endl << endl;
}

void EyeDiagram::push(node* v){
	push(v, tboot);
}

void EyeDiagram::push(node* v, Transition tran){
	double t = v->getTime();
	Transition transition = tnull;
	
	if (tran == tboot){
	//for MC result only
		double t = v->getTime();

		if (t <= 2e-10){
			Transition boot[5] = { t01, t11};
			int cycles = floor(t / period);
			transition = boot[cycles];
		}
		else{
			Transition boot[10] = { t10, t00, t01, t11, t11, t01, t00, t10, t11, t11 };
			t = t - 2e-10;
			int cycles = floor(t / period);
			cycles = cycles % 10;
			transition = boot[cycles];
		}

		/*
		Transition boot[5] = { t11, t11, t10, t00, t01 };
		//Transition boot[5] = { t01, t11, t10, t00, t01 };
		double t = v->getTime();
		int cycles = ceilf(t / period) - 1;
		cycles = cycles % 5;
		transition = boot[cycles];
		
		if (t < 1e-10){
			transition = t01;
		}
		else{
			transition = boot[cycles];
		}
		
		if (t<=8e-10 && t > 7e-10){
			transition = t01;
			cout << "Transition = " << transition << endl;
		}
		else if (t <= 9e-10 &&t > 8e-10){
			transition = t00;
			cout << "Transition = " << transition << endl;
		}
		else if (t <= 10e-10 &&t > 9e-10){
			transition = t10;
		}
		else if (t <= 11e-10 &&t > 10e-10){
			transition = t11;
		}


		*/

		cout << t << " " << transition << endl;
	}
	else{
		transition = tran;
	}
	
	int fullWindow = 2* (period / sampleRate);

	int i = int( round(t / sampleRate) ) % fullWindow;
	int j = int( round((t + period) / sampleRate)) % fullWindow;
	double volt = v->get(voltage);

	switch (transition){
	case t00:
		if (i < size){
			(palpebraInferior[i]).push_back(v);
			if (volt > maxInferior[i]){
				maxInferior[i] = volt;
				maxInferiorIndex[i] = palpebraInferior[i].size() - 1;
			}
			if (volt < minInferior[i]){
				minInferior[i] = volt;
				minInferiorIndex[i] = palpebraInferior[i].size() - 1;
			}
		}
		if (j < size){
			(palpebraInferior[j]).push_back(v);
			if (volt > maxInferior[j]){
				maxInferior[j] = volt;
				maxInferiorIndex[j] = palpebraInferior[j].size() - 1;
			}
			if (volt < minInferior[j]){
				minInferior[j] = volt;
				minInferiorIndex[j] = palpebraInferior[j].size() - 1;
			}
		}
		break;
	case t01:
		if (i < size){
			(palpebraSuperior[i]).push_back(v);
			if (volt > maxSuperior[i]){
				maxSuperior[i] = volt;
				maxSuperiorIndex[i] = palpebraSuperior[i].size() - 1;
			}
			if (volt < minSuperior[i]){
				minSuperior[i] = volt;
				minSuperiorIndex[i] = palpebraSuperior[i].size() - 1;
			}
		}
		if (j < size){
			(palpebraInferior[j]).push_back(v);
			if (volt > maxInferior[j]){
				maxInferior[j] = volt;
				maxInferiorIndex[j] = palpebraInferior[j].size() - 1;
			}
			if (volt < minInferior[j]){
				minInferior[j] = volt;
				minInferiorIndex[j] = palpebraInferior[j].size() - 1;
			}
		}
		break;
	case t10:
		if (i < size){
			(palpebraInferior[i]).push_back(v);
			if (volt > maxInferior[i]){
				maxInferior[i] = volt;
				maxInferiorIndex[i] = palpebraInferior[i].size() - 1;
			}
			if (volt < minInferior[j]){
				minInferior[i] = volt;
				minInferiorIndex[i] = palpebraInferior[i].size() - 1;
			}
		}
		if (j < size){
			(palpebraSuperior[j]).push_back(v);
			if (volt > maxSuperior[j]){
				maxSuperior[j] = volt;
				maxSuperiorIndex[j] = palpebraSuperior[j].size() - 1;
			}
			if (volt < minSuperior[j]){
				minSuperior[j] = volt;
				minSuperiorIndex[j] = palpebraSuperior[j].size() - 1;
			}
		}
		break;
	case t11:
		if (i < size){
			(palpebraSuperior[i]).push_back(v);
			if (volt > maxSuperior[i]){
				maxSuperior[i] = volt;
				maxSuperiorIndex[i] = palpebraSuperior[i].size()-1;
			}
			if (volt < minSuperior[i]){
				minSuperior[i] = volt;
				minSuperiorIndex[i] = palpebraSuperior[i].size() - 1;
			}
		}
		if (j < size){
			(palpebraSuperior[j]).push_back(v);
			if (volt > maxSuperior[j]){
				maxSuperior[j] = volt;
				maxSuperiorIndex[j] = palpebraSuperior[j].size() - 1;
			}
			if (volt < minSuperior[j]){
				minSuperior[j] = volt;
				minSuperiorIndex[j] = palpebraSuperior[j].size() - 1;
			}
		}
		break;
	default:
		cout << "Transition type not specified" << endl;
	}
}

void EyeDiagram::test(){
	//recompute the maxSuperior
	for (int i = 0; i < size; i++){
		double max = -100000;
		int index = -1;
		for (int j = 0; j < palpebraSuperior[i].size(); j++){
			if (palpebraSuperior[i][j]->get(voltage) > max){
				max = palpebraSuperior[i][j]->get(voltage);
				index = j;
			}
		}
		cout << "Max @" << i << " = " << max << " @" << index << "     == " << maxSuperior[i] << " @" << maxSuperiorIndex[i] << endl;
	}


	for (int i = 0; i < size; i++){
		if (palpebraSuperior[i].size() == 0){
			cout << "palpebraSuperior[" << i << "] is empty." << endl;
		}
		if (palpebraInferior[i].size() == 0){
			cout << "palpebraInferior[" << i << "] is empty." << endl;
		}

		if (minSuperiorIndex[i] == -1){
			cout << "minSuperiorIndex[" << i << "] is empty." << endl;
		}

		if (maxSuperiorIndex[i] == -1){
			cout << "maxSuperiorIndex[" << i << "] is empty." << endl;
		}
		if (minSuperiorIndex[i] == -1){
			cout << "minSuperiorIndex[" << i << "] is empty." << endl;
		}

		if (maxSuperiorIndex[i] == -1){
			cout << "maxSuperiorIndex[" << i << "] is empty." << endl;
		}

	}

}

string EyeDiagram::toString(){
	stringstream str;

	double vmin = -0.1, vmax = 1, tmin = 0, tmax = window;
	config->getParameter("edu.uiuc.csl.system.var.min", voltage, &vmin);
	config->getParameter("edu.uiuc.csl.system.var.max", voltage, &vmax);

	str << "plot [ " << tmin << ":" << tmax << "][" << vmin << ":" << vmax << "] 0 with linespoints lt \"white\" pt 0.01";
	str << " title \"" << " " << "\"  \n";
	
	//draw palpebra Superior
	for (int i = 0; i < size; i++){
		for (int j = 0; j < palpebraSuperior[i].size(); j++){
			if (!palpebraSuperior[i][j]->isRoot()){
				//double iToX = palpebraSuperior[i][j]->getTime();
				double iToX = i*sampleRate;
				double iToY = palpebraSuperior[i][j]->get(voltage);
				//double iFromX = palpebraSuperior[i][j]->getParent()->getTime();
				double iFromX = (i - 1)*sampleRate;
				double iFromY = palpebraSuperior[i][j]->getParent()->get(voltage);
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"blue\" lw 2 \n";
			}
		}
	}

	//draw palpebra inferior
	for (int i = 0; i < size; i++){	
		for (int j = 0; j < palpebraInferior[i].size(); j++){
			if (!palpebraInferior[i][j]->isRoot()){
				//double iToX = palpebraSuperior[i][j]->getTime();
				double iToX = i*sampleRate;
				double iToY = palpebraInferior[i][j]->get(voltage);
				//double iFromX = palpebraSuperior[i][j]->getParent()->getTime();
				double iFromX = (i - 1)*sampleRate;
				double iFromY = palpebraInferior[i][j]->getParent()->get(voltage);
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"red\" lw 2 \n";
			}
		}
	}
	
	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minInferior[i-1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"yellow\" lw 2 \n";
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxInferior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"black\" lw 2 \n";
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"green\" lw 2 \n";
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"purple\" lw 2 \n";
	}
	return str.str();
}

node* EyeDiagram::getNode(int i){
	return palpebraSuperior[i][minSuperiorIndex[i]];
}

int EyeDiagram::getWindowSize(){
	return size;
}

//todo
//autodetermining the window size