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
}

Plotter::~Plotter(){}

string drawLine(const double iFromX, const double iFromY,
				const double iToX, const double iToY){
					stringstream cmdstr;
					cmdstr << " set arrow from " << iFromX << "," << iFromY << " to " << iToX << "," << iToY << " nohead  lc rgb \"blue\" lw 2 \n" ;
					printf(cmdstr.str().c_str());
					return cmdstr.str();
}

string emptyPlot(string title, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, int angleX, int angleY, string xlabel, string ylabel, string zlabel){
	stringstream cmdstr;
	cmdstr << "splot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "][" << zmin << ":" << zmax <<"] 0 with linespoints lt \"white\" pt 0.01" ;
	cmdstr  << " title \"" << title << "\"  \n";
	cmdstr << "set xlabel \"$" << xlabel << "$\" \n";
	cmdstr << "set ylabel \"$" << ylabel << "$\" \n";
	cmdstr << "set zlabel \"$"<< zlabel<< "$\" \n";
	return cmdstr.str();
}

string emptyPlot(string title, double xmin, double xmax, double ymin, double ymax ){
	stringstream cmdstr;
	cmdstr << "plot [" << xmin << ":" << xmax << "][" << ymin << ":" << ymax << "] 0 with linespoints lt \"white\" pt 0.01" ;
	cmdstr  << " title \"" << title << "\"  \n";
	return cmdstr.str();
}


string drawLine(const double iFromX, const double iFromY, const double iFromZ, const double iToX, const double iToY, const double iToZ){
	stringstream cmdstr;
	cmdstr << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ <<  "     to     " << iToX << "," << iToY  << "," << iToZ << "        nohead  lc rgb \"blue\" lw 2 \n" ;
	//cout << " set arrow from " << iFromX << "," << iFromY << "," << iFromZ << " to " << iToX << "," << iToY  << "," << iToZ << " nohead  lc rgb \"blue\" lw 2 " ;
	return cmdstr.str();
}

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
		buffer = emptyPlot(title, xmin, xmax, ymin, ymax, zmin, zmax, 45, 45, xlabel, ylabel, zlabel);
		fprintf(gnuplotPipe, buffer.c_str());
	}else{
		buffer = emptyPlot(title, xmin, xmax, ymin, ymax);
		fprintf(gnuplotPipe, buffer.c_str());
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
				buffer = drawLine(x0, y0, z0, x1, y1, z1);
				out << " set arrow from " << x0 << "," << y0 << "," << z0 << " to " << x1 << "," << y1  << "," << z1 << " nohead  lc rgb \"blue\" lw 2 " << endl;
			}else{
				buffer = drawLine(x0, y0, x1, y1);
			}
			fprintf(gnuplotPipe, buffer.c_str());
			fflush(gnuplotPipe);
		}

		for(int i=0;i<n->getSize();i++){
			q.push(n->getChild(i));
		}
	}

	wait_for_key();
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
	buffer = emptyPlot("test", -1e-7, 100e-6, -0.5, 0.5  );
	fprintf(gnuplotPipe, buffer.c_str());
	fflush(gnuplotPipe);


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
			buffer = drawLine(t0, x0, t1, x1);

			fprintf(gnuplotPipe, buffer.c_str());
			fflush(gnuplotPipe);

		
		}
		for(int i=0;i<n->getSize();i++){
			q.push(n->getChild(i));
		}

	}
		buffer = "replot\n";
			fprintf(gnuplotPipe, buffer.c_str());
			fflush(gnuplotPipe);
	wait_for_key();

	_pclose(gnuplotPipe);
	out.close();

}
void Plotter::wait_for_key(){
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
