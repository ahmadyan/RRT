#include "eye.h"

//The constructor for the EyeDiagram classes, 
//this constructor will get the initial parameter from the config file, then allocate the eyelids and sets the initial values for max/min
EyeDiagram::EyeDiagram(Configuration* c){
	config = c;

	double unusuallySmallNumber = -9909;
	double unusuallyBigNumber = +9909;
	int nonExistentIndex = -1;

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
		maxSuperior[i] = unusuallySmallNumber;
		minSuperior[i] = unusuallyBigNumber;
		maxSuperiorIndex[i] = nonExistentIndex;
		minSuperiorIndex[i] = nonExistentIndex;
		maxInferior[i] = unusuallySmallNumber;
		minInferior[i] = unusuallyBigNumber;
		maxInferiorIndex[i] = nonExistentIndex;
		minInferiorIndex[i] = nonExistentIndex;
	}
		
	nadir = new double[8];

	nadir[0] = -0.2;	//area inside higher eyelid
	nadir[1] = 1.2;	//area inside higher eyelid
	nadir[2] = 4e-11;	//jitter-1
	nadir[3] = 8e-11;	//jitter-2
	nadir[4] = 4e-11;	//jitter-2
	nadir[5] = 8e-11;	//jitter-2
	nadir[6] = 1.2;	//area outside higher eyelid
	nadir[7] = -0.2;	//area outside lower eyelid


	vSize=20;	//for lebesge integrals dv
	dv = (0.8 - 0.2) / vSize;
	
	leftSuperior = new double[vSize];
	rightSuperior = new double[vSize];
	leftSuperiorIndex = new int[vSize];
	rightSuperiorIndex = new int[vSize];
	leftInferior = new double[vSize];
	rightInferior = new double[vSize];
	leftInferiorIndex = new int[vSize];
	rightInferiorIndex = new int[vSize];

	for (int i = 0; i < vSize; i++){
		leftSuperior[i] = unusuallySmallNumber;
		leftSuperiorIndex[i] = nonExistentIndex;
		
		rightSuperior[i] = unusuallyBigNumber;
		rightSuperiorIndex[i] = nonExistentIndex;
		
		leftInferior[i] = unusuallySmallNumber;
		leftInferiorIndex[i] = nonExistentIndex;

		rightInferior[i] = unusuallyBigNumber;
		rightInferiorIndex[i] = nonExistentIndex;
		//real programmer's don't comment!
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

	delete leftSuperior;
	delete rightSuperior;
	delete leftSuperiorIndex;
	delete rightSuperiorIndex;
	delete leftInferior;
	delete rightInferior;
	delete leftInferiorIndex;
	delete rightInferiorIndex;

	delete nadir;
}

//Perform the initial integration of the eyediagram after applying the 00110 pattern
void EyeDiagram::sum(){
	double* g = new double[8]();
	double dt = sampleRate;
	for (int i = 0; i < size; i++){
		g[0] += (minSuperior[i] - nadir[0])* dt;
		g[1] += (nadir[1] - maxSuperior[i]) * dt;
		g[6] += (nadir[6] - maxInferior[i])* dt;
		g[7] += (minInferior[i] - nadir[7]) * dt;
	}

	//Computing the lebesgue integrals
	for (int i = 0; i < vSize; i++){
		g[2] += (nadir[2] - leftSuperior[i])*dv;
		g[3] += (rightSuperior[i] - nadir[3])*dv;
		g[4] += (nadir[4]- leftInferior[i])*dv;
		g[5] += (rightInferior[i]-nadir[5])*dv;
	}
	cout << "Result" << endl << endl;
	for (int i = 0; i < 8; i++){
		cout << "g_" << i << " = " << g[i] << endl;
	}
	delete g;
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
	double volt = v->get(voltage);
	
	int cycles = ceilf(t / period) - 1;
	double t1 = t - period*cycles;
	double t2 = t1 + period;

	int i = int( round(t / sampleRate) ) % fullWindow;
	int j = int( round((t + period) / sampleRate)) % fullWindow;
	int k = floor((volt - 0.2) / dv);

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
			if (volt >= 0.2 && volt <= 0.8){
				if (t1 > leftSuperior[k]){
					leftSuperior[k] = t1;	
					leftSuperiorIndex[k] = palpebraSuperior[i].size() - 1;
				}
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
			if (volt >= 0.2 && volt <= 0.8){
				if (t2 < rightInferior[k]){
					rightInferior[k] = t2;
					rightInferiorIndex[k] = palpebraInferior[j].size() - 1;
				}
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
			if (volt >= 0.2 && volt <= 0.8){
				if (t1 > leftInferior[k]){
					leftInferior[k] = t1; 
					leftInferiorIndex[k] = palpebraInferior[i].size() - 1;
				}
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
			if (volt >= 0.2 && volt <= 0.8){
				if (t2 < rightSuperior[k]){
					rightSuperior[k] = t2;
					rightSuperiorIndex[k] = palpebraSuperior[j].size() - 1;
				}
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
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead lc rgb \"grey\" lw 2 \n";
			}
		}
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minInferior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"blue\" lw 2 \n";
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxInferior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"blue\" lw 2 \n";
	}




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

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"red\" lw 2 \n";	
	}

	for (int i = 1; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"red\" lw 2 \n";
	}


	

	return str.str();
}

enum EYE_FUNCTIONALS {G1, G2, G3, G4, G5, G7, G8};
node* EyeDiagram::getNode(int i){
	//todo: fix this
	return palpebraSuperior[i][minSuperiorIndex[i]];
}

int EyeDiagram::getWindowSize(){
	return size;
}

//todo
//autodetermining the window size