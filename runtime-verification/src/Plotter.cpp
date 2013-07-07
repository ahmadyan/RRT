#include "Plotter.h"
#include "TimedRRT.h"
#include <queue>
using namespace std;


//#include "gnuplot.h" //Gnuplot class handles POSIX-Pipe-communikation with Gnuplot

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)
#include <conio.h>   //for getch(), needed in wait_for_key()
#include <windows.h> //for Sleep()
void sleep(int i) { Sleep(i*1000); }
#endif


#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

//void wait_for_key(); // Programm halts until keypress

using std::cout;
using std::endl;


Plotter::Plotter(string path){
	gnuPlotPath = path ;
<<<<<<< HEAD
	#ifdef _WIN32
    gnuplotPipe = _popen(gnuPlotPath.c_str(),"w");
	#else
    gnuplotPipe = popen(gnuPlotPath.c_str(),"w");
	#endif
	emptyPlot("", 0, 1, 0, 1);
	closed=false;
=======
>>>>>>> origin/master
}

Plotter::~Plotter(){
    if(!closed) close();
}

<<<<<<< HEAD
void Plotter::close(){
    string buffer = "replot\n";
    fprintf(gnuplotPipe, buffer.c_str());
    fflush(gnuplotPipe);
    
    waitForKey();
    
#ifdef _WIN32
    _pclose(gnuplotPipe);
#else
    pclose(gnuplotPipe);
#endif
	closed=true;
}

string drawLine(const double iFromX, const double iFromY, const double iToX, const double iToY){
	stringstream cmdstr;
	cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
	printf(cmdstr.str().c_str());
	return cmdstr.str();
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel, string zlabel){
    stringstream cmdstr;
    cmdstr << "splot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "][" << zmin << ":" << zmax <<"] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    cmdstr << "set xlabel \"$" << xlabel << "$\" \n";
    cmdstr << "set ylabel \"$" << ylabel << "$\" \n";
    cmdstr << "set zlabel \"$"<< zlabel<< "$\" \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::emptyPlot(string title, double xmin, double xmax, double ymin, double ymax ){
    stringstream cmdstr;
    cmdstr << "plot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << title << "\"  \n";
    
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
=======
string drawLine(const double iFromX, const double iFromY,
				const double iToX, const double iToY){
					stringstream cmdstr;
					cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
					printf(cmdstr.str().c_str());
					return cmdstr.str();
>>>>>>> origin/master
}



void Plotter::drawLine(double iFromX, double iFromY, double iFromZ, double iToX, double iToY, double iToZ){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ <<  "     to     " << iToX << "," << iToY  << "," << iToZ << "        nohead  lc rgb \"blue\" lw 2 \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::drawLine(double iFromX, double iFromY, double iToX, double iToY){
    stringstream cmdstr;
    cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
    cout << cmdstr.str() << endl ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}


void Plotter::drawArray(vector< vector<double> > trace, int index1, int index2){
    if(trace.size()>0){
            //scaling
        
        double xmin=9999, xmax=-99999, ymin=9999, ymax=-9999;
        for(int i=0;i<trace.size();i++){
            if(trace.at(i).at(index1) > xmax) xmax=trace.at(i).at(index1);
            if(trace.at(i).at(index1) < xmin) xmin=trace.at(i).at(index1);
            if(trace.at(i).at(index2) > ymax) ymax=trace.at(i).at(index2);
            if(trace.at(i).at(index2) < ymin) ymin=trace.at(i).at(index2);

        }
        stringstream cmdstr;
        cmdstr << "plot [ " << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
        fprintf(gnuplotPipe, cmdstr.str().c_str());
        fflush(gnuplotPipe);

    for(int i=0;i<trace.size()-1;i++){
        vector<double> point1 = trace[i] ;
        vector<double> point2 = trace[i+1] ;
        drawLine( point1[index1], point1[index2], point2[index1], point2[index2] ) ;
    }
        
        
        stringstream cmdstr2;
        cmdstr2 << " replot \n" ;
        fprintf(gnuplotPipe, cmdstr2.str().c_str());
        fflush(gnuplotPipe);
    }
    else{ cout << "empty trace" << endl ;
    }
    
    
}

/*void Plotter::drawTrace(Trace* trace, int index1, int index2){
    drawArray(trace->getData(), index1, index2);
}*/


/*void Plotter::drawTrace(Trace* trace, double t_max, int index){
    
    double min=9999, max=-99999;
    for(int i=0;i<trace->getSampleSize();i++){
        if(trace->getSample(i)[index] > max) max=trace->getSample(i)[index];
        if(trace->getSample(i)[index] < min) min=trace->getSample(i)[index];
    }
    stringstream cmdstr;
    cmdstr << "plot [ " << 0 << ":" << t_max << "][" << min << ":" << max << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << " " << "\"  \n";
    //cmdstr << "set ylabel " << "y(" << index << " )"  << " \n";
    //cmdstr << "set xlabel " << "time" << " \n";
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
    
    double time_step = t_max / trace->getSampleSize();
    for(int i=0;i<trace->getData().size()-1;i++){
        vector<double> point1 = trace->getSample(i) ;
        vector<double> point2 = trace->getSample(i+1) ;
        drawLine( i*time_step, point1[index], (i+1)*time_step, point2[index] ) ;
    }
    stringstream cmdstr2;
    cmdstr2 << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr2.str().c_str());
    fflush(gnuplotPipe);
}*/



/*string drawLine(const double iFromX, const double iFromY, const double iFromZ, const double iToX, const double iToY, const double iToZ){
	stringstream cmdstr;
	cmdstr << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ <<  "     to     " << iToX << "," << iToY  << "," << iToZ << "        nohead  lc rgb \"blue\" lw 2 \n" ;
	//cout << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ << " to " << iToX << "," << iToY  << "," << iToZ << " nohead  lc rgb \"blue\" lw 2 " ;
	return cmdstr.str();
}*/

void Plotter::plotRRT(string name, string title, string output, RRT rrt, string xlabel, string ylabel, string zlabel){
	if(rrt.getDimension()>3){
		cout << " Error - This is 2D/3D drawing method." << endl ;
		return ;
	}

	bool rrt3d = false ;
	if( rrt.getDimension()==3) rrt3d = true;
	ofstream out("plot.txt");
	FILE *gnuplotPipe = _popen(gnuPlotPath.c_str(),"w");
	string buffer ;
	fflush(gnuplotPipe);
	double xmin = rrt.getMin(0);
	double xmax = rrt.getMax(0);
	double ymin = rrt.getMin(1);
	double ymax = rrt.getMax(1);
	double zmin = 0 ;
	double zmax = 0 ;
	if (rrt3d) {
		zmin = rrt.getMin(2);
		zmax = rrt.getMax(2);
		emptyPlot(title, xmin, xmax, ymin, ymax, zmin, zmax, 45, 45, xlabel, ylabel, zlabel);
	}else{
		emptyPlot(title, xmin, xmax, ymin, ymax);
	}
	fflush(gnuplotPipe);

	cout << "ZMAX=" <<  zmax << endl ;

	queue<node*> q ; 
	q.push(rrt.getRoot());
	while(!q.empty()){
		node* n = q.front();
		q.pop();

		if(!n->isRoot()){
			double* begin = n->get();
			double* end   = n->getParent()->get();
			int xscale = 1, yscale=1, zscale=1 ;
			double x0 = xscale*begin[0] ;
			double y0 = yscale*begin[1] ;
			double x1 = xscale*end[0] ;
			double y1 = yscale*end[1] ;
			double z0 = 0 ;
			double z1 = 0 ;
			if(rrt3d){
				z0 = zscale*begin[2];
				z1 = zscale*end[2] ;
				drawLine(x0, y0, z0, x1, y1, z1);
				out << " set arrow from " << x0 << "," << y0 << "," << z0 << " to " << x1 << "," << y1  << "," << z1 << " nohead  lc rgb \"blue\" lw 2 " << endl;
			}else{
				drawLine(x0, y0, x1, y1);
			}
			fflush(gnuplotPipe);
		}

		for(int i=0;i<n->getSize();i++){
			q.push(n->getChild(i));
		}
	}

	waitForKey();
	_pclose(gnuplotPipe);
	out.close();
}

void Plotter::plotTrace(RRT rrt, int v1, int v2, int tdim, double simTime, double dt){
	ofstream out("plot.txt");
	FILE *gnuplotPipe = _popen(gnuPlotPath.c_str(),"w");
	string buffer ;
	fflush(gnuplotPipe);
	//double simTime = ((TimedRRT)(rrt)).getSimTime();
	int n = (int)(simTime / dt) ;
	double* min = new double[n+1];
	double* max = new double[n+1];
	for(int i=0;i<n;i++){
		min[i]=+1;
		max[i]=-1;
	}
	emptyPlot("test", -1e-7, 100e-6, -0.5, 0.5  );
	
	queue<node*> q ; 
	q.push(rrt.getRoot());
	while(!q.empty()){
		node* n = q.front();
		q.pop();

		if(!n->isRoot()){
			double* begin = n->get();
			double* end   = n->getParent()->get();

			int xscale = 1, yscale=1, zscale=1 ;
			double x1 = xscale*(begin[v1]-begin[v2]) ;
			double t1 = yscale*begin[tdim] ;
			double x0 = xscale*(end[v1]-end[v2]) ;
			double t0 = yscale*end[tdim] ;


			int it0 = t0/dt; int it1 = t1/dt;
			if(x1>max[it1]) max[it1]=x1 ;
			if(x1<min[it1]) min[it1]=x1 ;
			if(x0>max[it0]) max[it0]=x0 ;
			if(x0<min[it0]) min[it0]=x0 ;
			cout << "**" << endl ;
//			cout << t0 << " " << t1 << " " << x0 << " " << x1 << endl ;
<<<<<<< HEAD
			drawLine(t0, x0, t1, x1);
			fflush(gnuplotPipe);

=======
			buffer = drawLine(t0, x0, t1, x1);

			fprintf(gnuplotPipe, buffer.c_str());
			fflush(gnuplotPipe);

>>>>>>> origin/master
		
		}
		for(int i=0;i<n->getSize();i++){
			q.push(n->getChild(i));
		}

	}
<<<<<<< HEAD
	buffer = "replot\n";
	fprintf(gnuplotPipe, buffer.c_str());
	fflush(gnuplotPipe);
	waitForKey();
=======
		buffer = "replot\n";
			fprintf(gnuplotPipe, buffer.c_str());
			fflush(gnuplotPipe);
	wait_for_key();
>>>>>>> origin/master

	_pclose(gnuplotPipe);
	out.close();

}

// draws the given signal w.r.t. time to gnuplot
void Plotter::drawTrace(vector<double> signal, double t){
	int n = signal.size() ;
	double dt = t / signal.size();
	//find the min/max
	double min=9999, max=-99999;
    for(int i=0;i<n;i++){
        if(signal[i] > max) max=signal[i];
        if(signal[i] < min) min=signal[i];
    }

	stringstream cmdstr;
    cmdstr << "plot [ " << 0 << ":" << t << "][" << min << ":" << max << "] 0 with linespoints lt \"white\" pt 0.01" ;
    cmdstr  << " title \"" << " " << "\"  \n";
    //cmdstr << "set ylabel " << "y(" << index << " )"  << " \n";
    //cmdstr << "set xlabel " << "time" << " \n";
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);

	for(int i=0;i<n-1;i++){
		drawLine(i*dt, signal[i], (i+1)*dt, signal[i+1]);
	}
	 
	stringstream cmdstr2;
    cmdstr2 << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr2.str().c_str());
    fflush(gnuplotPipe);
	waitForKey();
}


void Plotter::waitForKey(){
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__TOS_WIN__)  // every keypress registered, also arrow keys
    cout << endl << "Press any key to continue..." << endl;
    
    FlushConsoleInputBuffer(GetStdHandle(STD_INPUT_HANDLE));
    _getch();
#elif defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
    cout << endl << "Press ENTER to continue..." << endl;
    
    std::cin.clear();
    std::cin.ignore(std::cin.rdbuf()->in_avail());
    std::cin.get();
#endif
    return;
}

void Plotter::saveToPdf(string path){
    stringstream cmdstr;
    cmdstr << " set term post \n" ;
    cmdstr << " set output \"" << path << "\"\n" ;
    cmdstr << " replot \n" ;
    fprintf(gnuplotPipe, cmdstr.str().c_str());
    fflush(gnuplotPipe);
}

void Plotter::execute(string str){
    fprintf(gnuplotPipe, str.c_str());
    fflush(gnuplotPipe);
    cout << str << endl ;
    string buffer = "replot\n";
    fprintf(gnuplotPipe, buffer.c_str());
    fflush(gnuplotPipe);
    
}
