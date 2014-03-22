#include "eye.h"
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>

#include "node.h"
//The constructor for the EyeDiagram classes, 
//this constructor will get the initial parameter from the config file, then allocate the eyelids and sets the initial values for max/min
EyeDiagram::EyeDiagram(Configuration* c){
	config = c;
	bitsAreSet = 0;
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


	vSize = 10;	//for lebesge integrals dv
	dv = (0.8 - 0.2) / vSize;

	leftSuperior = new double[vSize];
	rightSuperior = new double[vSize];
	leftSuperiorIndex = new int[vSize];
	rightSuperiorIndex = new int[vSize];
	leftInferior = new double[vSize];
	rightInferior = new double[vSize];
	leftInferiorIndex = new int[vSize];
	rightInferiorIndex = new int[vSize];
	leftSuperiorNodes = new node*[vSize];
	rightSuperiorNodes = new node*[vSize];
	leftInferiorNodes = new node*[vSize];
	rightInferiorNodes = new node*[vSize];

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

	maxOne = unusuallySmallNumber;
	minOne = unusuallyBigNumber;
	maxZero = unusuallySmallNumber;
	minZero = unusuallyBigNumber;
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
	for (int i = 1; i < size; i++){
		g[0] += (minSuperior[i] - nadir[0])* dt;
		g[1] += (nadir[1] - maxSuperior[i]) * dt;
		g[6] += (nadir[6] - maxInferior[i])* dt;
		g[7] += (minInferior[i] - nadir[7]) * dt;
	}

	//Computing the lebesgue integrals
	for (int i = 1; i < vSize; i++){
		g[2] += (nadir[2] - leftSuperior[i])*dv;
		g[3] += (rightSuperior[i] - nadir[3])*dv;
		g[4] += (nadir[4] - leftInferior[i])*dv;
		g[5] += (rightInferior[i] - nadir[5])*dv;
	}
	//cout << "Result" << endl << endl;
	//for (int i = 0; i < 8; i++){
	//	cout << "g_" << i << " = " << g[i] << endl;
	//}


	//cout << "Classic measurements:" << endl;
	double noiseMargin = minSuperior[size / 2] - maxInferior[size / 2];
	double jitterMargin1 = rightInferior[vSize / 2] - leftInferior[vSize / 2];
	double jitterMargin2 = rightInferior[vSize / 2] - leftInferior[vSize / 2];
	double zeroCrossingRight = leftSuperior[vSize / 2] - rightInferior[vSize / 2] + period;
	double zeroCrossingLeft = leftInferior[vSize / 2] - rightSuperior[vSize / 2] + period;
	/*
	double maxsup = -1;
	for (int i = 0; i < size; i++){
	if ((maxSuperior[i] - minSuperior[i])>maxsup){
	maxsup = maxSuperior[i] - minSuperior[i];
	}
	if ((maxInf))
	}*/

	//cout << "Noise margin= " << noiseMargin << endl;
	//cout << "jitterMargin(1) = " << jitterMargin1 << endl;
	//cout << "jitterMargin(2) = " << jitterMargin2 << endl;
	//cout << "zeroCrossingRight = " << zeroCrossingRight << endl;
	//cout << "zeroCrossingLeft = " << zeroCrossingLeft << endl;

	if (config->checkParameter("edu.uiuc.crhc.core.options.eyediagram.dump", "1")){
		//print the current stat for the eye diagram
		ofstream file;
		file.open(config->get("edu.uiuc.crhc.core.options.eyediagram.dump.filename"), std::ofstream::app);

		for (int i = 0; i < 8; i++){
			file << g[i] << " ";
		}

		file << noiseMargin << " ";
		file << jitterMargin1 << " ";
		file << jitterMargin2 << " ";
		file << zeroCrossingRight << " ";
		file << zeroCrossingLeft;
		file << endl;

		file.close();
	}
	delete g;
}

void EyeDiagram::push(node* v){
	push(v, tboot);
	if (config->checkParameter("edu.uiuc.crhc.core.options.eyediagram.dump", "1"))
		sum();
}

void EyeDiagram::push(node* v, Transition tran){
	double t = v->getTime();
	Transition transition = getTransition(v);

	int fullWindow = 2 * (period / sampleRate);
	double volt = v->get(voltage);

	int cycles = ceilf(t / period) - 1;
	if (cycles == -1)cycles = 0;
	double t1 = t - period*cycles;
	double t2 = t1 + period;
	double x = floorf( t1 / sampleRate ); 

	int periodWindow = round(period / sampleRate);
	int i = round(t1 / sampleRate);
	int j = i + periodWindow;

	cout << periodWindow << " " << i << " " << j << endl;
	//int i = int(floorf(t / sampleRate)) % fullWindow;
	//int j = int(floorf(t + period) / sampleRate) % fullWindow;



	int k = floor((volt - 0.2) / dv);

	switch (transition){
	case t00:
		if (volt > maxZero) maxZero = volt;
		if (volt < minZero) minZero = volt;
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
		if (v->getJitter() == 1){
			jitterFrontierSet01.push_back(v);
		}
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
					leftSuperiorNodes[k] = v;
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
					rightInferiorNodes[k] = v;
				}
			}
		}
		break;
	case t10:
		if (v->getJitter() == 1){
			jitterFrontierSet10.push_back(v);
		}
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
			if (volt >= 0.2 && volt <= 0.8){
				if (t1 > leftInferior[k]){
					leftInferior[k] = t1;
					leftInferiorIndex[k] = palpebraInferior[i].size() - 1;
					leftInferiorNodes[k] = v;
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
					rightSuperiorNodes[k] = v;
				}
			}
		}
		break;
	case t11:
		if (volt>maxOne) maxOne = volt;
		if (volt < minOne) minOne = volt;
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
	for (int i = 1; i < vSize; i++){
		if (leftSuperiorIndex[i] == -1){
			leftSuperiorIndex[i] = leftSuperiorIndex[i - 1];
			leftSuperior[i] = leftSuperior[i - 1];
			leftSuperiorNodes[i] = leftSuperiorNodes[i - 1];
		}
		if (rightSuperiorIndex[i] == -1){
			rightSuperiorIndex[i] = rightSuperiorIndex[i - 1];
			rightSuperior[i] = rightSuperior[i - 1];
			rightSuperiorNodes[i] = rightSuperiorNodes[i - 1];
		}
		if (leftInferiorIndex[i] == -1){
			leftInferiorIndex[i] = leftInferiorIndex[i - 1];
			leftInferior[i] = leftInferior[i - 1];
			leftInferiorNodes[i] = leftInferiorNodes[i - 1];
		}
		if (rightInferiorIndex[i] == -1){
			rightInferiorIndex[i] = rightInferiorIndex[i - 1];
			leftInferiorNodes[i] = leftInferiorNodes[i - 1];
			rightInferior[i] = rightInferior[i - 1];
		}
	}

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


	for (int i = 0; i < vSize; i++){
		if (leftSuperiorIndex[i] == -1) cout << "leftSuperiorIndex " << i << " is empty " << endl;
		if (rightSuperiorIndex[i] == -1) cout << "rightSuperiorIndex!" << i << " is empty " << endl;
		if (leftInferiorIndex[i] == -1) cout << "leftInferiorIndex!" << i << " is empty " << endl;
		if (rightInferiorIndex[i] == -1) cout << "rightInferiorIndex!" << i << " is empty " << endl;
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
				str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead lc rgb \"blue\" lw 2 \n";
			}
		}
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

	for (int i = 2; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minInferior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"red\" lw 2 \n";
	}

	for (int i = 2; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxInferior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxInferior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"black\" lw 2 \n";
	}


	for (int i = 2; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = minSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = minSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"green\" lw 2 \n";
	}

	for (int i = 2; i < size; i++){
		double iToX = i*sampleRate;
		double iToY = maxSuperior[i];
		double iFromX = (i - 1)*sampleRate;
		double iFromY = maxSuperior[i - 1];
		str << " set arrow from " << iFromX << "," << iFromY << "   to     " << iToX << "," << iToY << "  nohead  lc rgb \"yellow\" lw 2 \n";
	}
	return str.str();
}

//	Input determines which objective
//	0: jitter (0->1)
//	1: jitter (1->0)
//	2: inside 1
//	3: inside 0
//	4: outside 1
//	5: outside 0
//	6: any of the lebesgue
node* EyeDiagram::getNode(int i){
	int x = -1, y = -1;
	switch (i){
	case 0:
		//int q = rand() % jitterFrontierSet01.size();
		//return jitterFrontierSet01[q];
		return jitterFrontierSet01[jitterFrontierSet01.size() - 1];
	case 1:
		//int q = rand() % jitterFrontierSet10.size();
		//return jitterFrontierSet10[q];
		return jitterFrontierSet10[jitterFrontierSet10.size() - 1];
	case 2:
		x = rand() % palpebraSuperior.size();
		return palpebraSuperior[x][minSuperiorIndex[x]];
		break;
	case 3:
		x = rand() % palpebraInferior.size();
		return palpebraInferior[x][maxInferiorIndex[x]];
		break;
	case 4:
		x = rand() % palpebraSuperior.size();
		return palpebraSuperior[x][maxSuperiorIndex[x]];
		break;
	case 5:
		x = rand() % palpebraSuperior.size();
		return palpebraInferior[x][minInferiorIndex[x]];
		break;
	case 6:
		x = rand() % 4;
		y = rand() % vSize;

		if (x == 0){
			return leftInferiorNodes[y];
		}
		else if (x == 1){
			return rightInferiorNodes[y];
		}
		else if (x == 2){
			return leftSuperiorNodes[y];
		}
		else{
			return rightSuperiorNodes[y];
		}
		break;
	}
}

int EyeDiagram::getWindowSize(){
	return size;
}

void EyeDiagram::setBits(vector<int> b){
	for (int i = 0; i < b.size(); i++){
		bits.push_back(b[i]);
	}
	bitsAreSet = 1;
}

//This is hard-coded, fix it.
Transition EyeDiagram::getTransition(node* v){
	Transition transition;
	if (bitsAreSet == 1){
		double t = v->getTime();
		int cycles = ceilf(t / period) - 1;
		cout << cycles << endl;
		if (t <= period){
			return t01;	//forst bit
		}

		//bits has the input, in inverter we are going to invert it.
		int currentBit = 1-bits[cycles];
		int previousBit = 1-bits[cycles - 1];

		if (previousBit == 0 && currentBit == 0){
			return t00;
		}
		else if (previousBit == 0 && currentBit == 1){
			return t01;
		}
		else if (previousBit == 1 && currentBit == 0){
			return t10;
		}
		else if (previousBit == 1 && currentBit == 1){
			return t11;
		}
		else{
			cout << "Uknown transition " << currentBit << " " << previousBit << endl;
			return tnull;
		}

	}
	else{
		double t = v->getTime();

		if (t < 2e-10){
			Transition boot[5] = { t01, t11 };
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
	}
}

//todo
//autodetermining the window size